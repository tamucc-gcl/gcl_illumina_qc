// workflows/optimize_denovo.nf
// ============================================================================
// De novo cutoff/similarity OPTIMIZATION subworkflow (replaces the map-back sweep)
// ============================================================================
// SKELETON / CHUNK 2: channel topology with STUB processes (echo placeholders).
// This exists to validate the three-branch wiring on a cheap -resume run BEFORE
// any expensive process internals are written. Real signal processes are filled
// in chunks 3-6.
//
// Three branches off the repaired reads (passed in as cleaned_reads):
//   1. PRODUCTION  — intact individuals; handled in main.nf, NOT here.
//   2. ASSEMBLY    — intact subset -> per-sample uniq.seqs -> candidate references.
//   3. CONCORDANCE — N highest-depth individuals split 50/50 (pseudo-reps),
//                    used ONLY for the stage-2 concordance signal. Never published.
//
// The assembly subset and concordance individuals MAY OVERLAP (default).
//
// Emits (contract that denovo_assembly.nf -> main.nf depends on):
//   reference, assembly_stats   (always)
//   plus diagnostics/optimization artifacts for the report (chunk 7).

nextflow.enable.dsl = 2

include { extract_unique_seqs }  from '../modules/extract_unique_seqs.nf'
include { assembly_diagnostics } from '../modules/assembly_diagnostics.nf'
include { split_pseudo_rep }     from '../modules/split_pseudo_rep.nf'
include { filter_unique_seqs_candidate } from '../modules/filter_unique_seqs_candidate.nf'
include { assemble_rainbow_candidate }   from '../modules/assemble_rainbow_candidate.nf'
include { compute_cheap_signals }       from '../modules/compute_cheap_signals.nf'
include { fit_nb_mixture }               from '../modules/fit_nb_mixture.nf'
include { provisional_rank }             from '../modules/provisional_rank.nf'

// ---------------------------------------------------------------------------
// STUB PROCESSES — echo placeholders. Replaced in later chunks.
// Each prints what it WOULD do and emits correctly-typed placeholder files so
// the channel graph is exercised end to end.
// ---------------------------------------------------------------------------

// (chunk 3a) stub_build_candidate REMOVED — real candidates now come from
// filter_unique_seqs_candidate + assemble_rainbow_candidate.

// (chunk 3b) stub_cheap_signals and stub_provisional_rank REMOVED — real signals
// now come from compute_cheap_signals.nf + fit_nb_mixture.nf + provisional_rank.nf.

// Stand-in for: stage-2 bcftools SNP recovery + pseudo-rep concordance per survivor.
// Real version = chunk 5. Takes a survivor candidate + the pseudo-rep reads.
process stub_stage2 {
    label 'basic'
    tag "${meta.id}"
    input:
        tuple val(meta), path(candidate), path(pseudo_reads)
    output:
        path("stage2_${meta.id}.tsv"), emit: stage2_row
    script:
    """
    echo "[STUB] stage-2 bcftools for survivor ${meta.id} using \$(ls *.fq.gz 2>/dev/null | wc -l) pseudo-rep fastqs"
    # id n_snps concordance
    printf "%s\\t0\\t0\\n" "${meta.id}" > stage2_${meta.id}.tsv
    """
}

// Stand-in for: weight-free rank aggregation -> winning meta.id + summary/plot.
// Real version = chunk 6.
process stub_aggregate {
    label 'basic'
    tag "aggregate"
    input:
        path(provisional)
        path(stage2_rows)
    output:
        path("best_id.value"),       emit: best_id
        path("optimize_summary.tsv"),emit: summary
        path("optimize_plot.png"),   emit: plot
    script:
    """
    echo "[STUB] aggregate provisional + stage-2 rows"
    cat ${provisional} > optimize_summary.tsv
    # STUB winner = first survivor row's id from stage-2
    cat stage2_*.tsv | head -n1 | cut -f1 > best_id.value
    echo "[STUB] best_id:"; cat best_id.value
    echo "stub plot" > optimize_plot.png
    """
}

// Stand-in for: SELECT the winning candidate as the production reference and
// publish it. Real version (chunk 6) — under Position B, candidates are already
// full-sample assemblies, so finalize does NOT re-assemble. It picks the winner's
// existing denovo_reference.fa by meta.id and publishes it to a clean, stable
// location (the single publishDir for the production reference; candidate
// publishDirs will be removed once finalize owns publishing).
process stub_finalize_reference {
    label 'basic'
    tag "finalize:${meta.id}"

    // Position B: this is the SINGLE publish point for the production reference.
    publishDir "${params.outdir}/denovo_assembly", mode: params.publish_dir_mode

    input:
        tuple val(meta), path(winning_reference, stageAs: 'winner_input.fa')
    output:
        path("denovo_reference.fa"), emit: reference
        path("assembly_stats.txt"),  emit: stats
    script:
    """
    echo "Selected winning de novo reference: ${meta.id} (cutoff1=${meta.c1} cutoff2=${meta.c2} sim=${meta.sim})"
    # Winner is already a full-sample assembly — stage it under the canonical name.
    cp winner_input.fa denovo_reference.fa
    echo "Final reference contigs: \$(grep -c '^>' denovo_reference.fa)" > assembly_stats.txt
    cat assembly_stats.txt
    """
}


workflow optimize_denovo {
    take:
        cleaned_reads          // tuple(sample_id, r1, r2) — repaired reads (intact)

    main:
        // ----- ASSEMBLY BRANCH: ALL samples -> per-sample uniq.seqs -----
        // Position B: candidate references are built from the FULL sample set, NOT
        // a subset. Each candidate is therefore a real, production-grade assembly;
        // the winner IS the production reference (finalize just selects/publishes
        // it — no re-assembly). This makes cutoff2 (n-individuals) mean the SAME
        // thing during optimization and in the final reference, which a subset
        // would distort. The assembly seed/subset are gone here; optimize_seed and
        // snp_sample_pct now govern the STAGE-2 SNP-recovery subset only (chunk 5).
        assembly_all = cleaned_reads

        extract_unique_seqs( assembly_all )
        all_uniq_seqs = extract_unique_seqs.out.uniq_seqs
            .map{ sid, f -> f }
            .collect()

        // Diagnostics: knee values (grid CEILINGS) + curves + NB-mixture (5b later)
        assembly_diagnostics( all_uniq_seqs )

        // ----- Build the (c1,c2,sim) grid from floor..knee (or pinned params) -----
        // Knee values are read from diagnostics .value files and act as CEILINGS.
        // Explicit --cutoff1/--cutoff2 PIN that dimension to a single value.
        // Independent floors: cutoff1_floor (per-individual coverage) and
        // cutoff2_floor (n individuals). cutoff2_floor defaults to 3 — the k2 corner
        // is biologically junky AND the slowest to assemble (CD-HIT cost).
        def c1_floor = (params.cutoff1_floor ?: 2) as int
        def c2_floor = (params.cutoff2_floor ?: 3) as int

        c1_vals = (params.cutoff1 != null)
            ? Channel.value( [ params.cutoff1 as int ] )
            : assembly_diagnostics.out.cutoff1_value
                .map { f -> def k = f.text.trim() as int; (c1_floor..Math.max(c1_floor, k)).toList() }
        c2_vals = (params.cutoff2 != null)
            ? Channel.value( [ params.cutoff2 as int ] )
            : assembly_diagnostics.out.cutoff2_value
                .map { f -> def k = f.text.trim() as int; (c2_floor..Math.max(c2_floor, k)).toList() }

        // Cutoff combos = c1 x c2 (filtering is similarity-INDEPENDENT, so filter
        // runs ONCE per (c1,c2) and is later crossed with similarities).
        // Wrap each list with .map{ [it] } so combine keeps it as ONE slot rather
        // than spreading the list elements across tuple positions.
        cutoff_combos = c1_vals.map{ [it] }
            .combine( c2_vals.map{ [it] } )
            .flatMap { c1list, c2list ->
                def out = []
                for (c1 in c1list) for (c2 in c2list) {
                    out << [ c1: c1 as int, c2: c2 as int, id_cutoff: "c${c1}_k${c2}".toString() ]
                }
                out
            }

        // Pair each cutoff meta with the collected uniq.seqs (wrap list -> one slot)
        filter_in = cutoff_combos
            .combine( all_uniq_seqs.map{ [it] } )
            .map { meta, files -> tuple(meta, files) }

        filter_unique_seqs_candidate( filter_in )

        // Cross each filtered (c1,c2) result with the similarity grid -> full candidate set.
        sims = Channel.fromList( params.optimize_cluster_similarity )

        assemble_in = filter_unique_seqs_candidate.out.filtered
            .combine( sims )
            .map { meta, uniq_fasta, totaluniqseq, fstats, sim ->
                def m = meta + [ sim: sim, id: "${meta.id_cutoff}_s${sim}".toString() ]
                tuple(m, uniq_fasta, totaluniqseq)
            }

        assemble_rainbow_candidate(
            assemble_in,
            params.div_f, params.div_K, params.merge_r, params.final_similarity
        )

        candidates = assemble_rainbow_candidate.out.reference   // tuple(meta, denovo_reference.fa)

        // ----- STAGE 1: CHEAP signals per candidate -> provisional rank -----
        // 1a inflection inputs + 1b redundancy + 2 anchor input, per candidate.
        // redundancy self-cluster identity: tighter than any grid sim so it surfaces
        // residual near-duplicate contigs (default 0.98).
        def redun_identity = (params.optimize_redundancy_identity ?: 0.98) as double
        compute_cheap_signals( candidates, redun_identity )
        cheap_rows = compute_cheap_signals.out.cheap_row.collect()

        // 5b: global NB-mixture cutoff1 from the pooled coverage distribution
        // (assembly_diagnostics now emits coverage_freq.txt; knee value is fallback).
        fit_nb_mixture(
            assembly_diagnostics.out.coverage_freq,
            assembly_diagnostics.out.cutoff1_value
        )

        // 2 anchor: expected RAD locus count from enzyme/genome/size-selection params.
        // Cut-site model: each enzyme cuts ~once per 4^site_len bp; double-digest
        // fragments bounded by the two cut frequencies; fraction in the size window
        // approximated by window_width / mean_fragment_length. If any input missing,
        // expected = "NA" and signal 2 is dropped in the rank script.
        def expected_loci = "NA"
        if (params.genome_size_est && params.size_select_min && params.size_select_max) {
            def G   = params.genome_size_est as double
            def s1  = (params.enzyme1_site_len ?: 6) as int
            def s2  = (params.enzyme2_site_len ?: 4) as int
            // expected cut sites for each enzyme across the genome (both strands ~ /2 cancels in ratio)
            def cuts1 = G / Math.pow(4, s1)
            def cuts2 = G / Math.pow(4, s2)
            // double-digest fragments flanked by one of each enzyme: density ~ 2*cuts1*cuts2/G
            // (Poisson-spacing approximation for AB/BA fragments along the genome)
            def dd_fragments = 2.0 * cuts1 * cuts2 / G
            def mean_frag = G / dd_fragments
            def win = (params.size_select_max as double) - (params.size_select_min as double)
            // fraction of an exponential fragment-length distribution within the window,
            // centered near the window: approximate retained fraction = win / mean_frag, capped at 1
            def frac = Math.min(1.0, win / mean_frag)
            def exp_loci = Math.round(dd_fragments * frac)
            expected_loci = exp_loci.toString()
            log.info "Biological anchor: expected ~${expected_loci} RAD loci (enzymes ${s1}/${s2}-cutter, genome ${G} bp, window ${params.size_select_min}-${params.size_select_max} bp)"
        } else {
            log.info "Biological anchor (signal 2) disabled — supply genome_size_est + size_select_min/max to enable"
        }

        def min_surv = (params.stage2_min_candidates ?: 3) as int
        provisional_rank( cheap_rows, fit_nb_mixture.out.nb_cutoff1, expected_loci, min_surv )
        survivors_ch = provisional_rank.out.survivors



        // ----- CONCORDANCE BRANCH: N highest-depth individuals, split 50/50 -----
        // STUB selection: deterministic first-N by SORTED sample id (chunk 5 will
        // replace with depth-based selection). Sort before take so the set is
        // stable across runs (same cache-determinism issue as the assembly subset).
        n_reps = params.n_pseudo_reps ?: 6
        concordance_indivs = cleaned_reads
            .toList()
            .flatMap { rows -> new ArrayList(rows).sort { a, b -> a[0] <=> b[0] }.take( n_reps as int ) }
        split_pseudo_rep( concordance_indivs, params.optimize_seed ?: 42 )

        // Flatten pseudo-rep halves into a collected bag of fastqs for stage-2.
        // Each emission: (parent, a_id, a_r1, a_r2, b_id, b_r1, b_r2)
        pseudo_fastqs = split_pseudo_rep.out.pseudo_reps
            .flatMap{ parent, aid, ar1, ar2, bid, br1, br2 -> [ar1, ar2, br1, br2] }
            .collect()

        // ----- STAGE 2: bcftools on SURVIVORS only -----
        // Join survivor ids back to their candidate fastas, attach pseudo-rep reads.
        survivor_ids = survivors_ch
            .splitText()
            .map{ it.trim() }
            .filter{ it }

        survivor_candidates = candidates
            .map{ meta, fa -> tuple(meta.id.toString(), meta, fa) }
            .join( survivor_ids.map{ id -> tuple(id.toString(), true) } )
            .map{ id, meta, fa, flag -> tuple(meta, fa) }

        stage2_in = survivor_candidates
            .combine( pseudo_fastqs.map{ [it] } )
            .map{ meta, fa, reads -> tuple(meta, fa, reads) }

        stub_stage2( stage2_in )
        stage2_rows = stub_stage2.out.stage2_row.collect()

        // ----- AGGREGATE: weight-free rank aggregation -> winner -----
        stub_aggregate( provisional_rank.out.provisional, stage2_rows )

        // ----- FINALIZE: SELECT the winning candidate as the production reference -----
        // Position B: candidates are already full-sample assemblies, so finalize just
        // selects the winner — it does NOT collect all candidates (they're all named
        // denovo_reference.fa, which collides) and does NOT re-assemble. Join the
        // winning id back to its single candidate and pass only that one through.
        best_id_ch = stub_aggregate.out.best_id
            .map { f -> f.text.trim() }

        winning_candidate = candidates
            .map { meta, fa -> tuple(meta.id.toString(), meta, fa) }
            .combine( best_id_ch )
            .filter { id, meta, fa, best -> id == best.toString() }
            .map { id, meta, fa, best -> tuple(meta, fa) }

        stub_finalize_reference( winning_candidate )

    emit:
        reference      = stub_finalize_reference.out.reference
        assembly_stats = stub_finalize_reference.out.stats
        // Optimization/diagnostics artifacts (consumed by report in chunk 7):
        cutoff1_plot   = assembly_diagnostics.out.cutoff1_plot
        cutoff2_plot   = assembly_diagnostics.out.cutoff2_plot
        diag_summary   = assembly_diagnostics.out.summary
        optimize_summary = stub_aggregate.out.summary
        optimize_plot    = stub_aggregate.out.plot
}
