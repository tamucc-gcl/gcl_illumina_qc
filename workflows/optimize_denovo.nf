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
include { compute_snp_signals }          from '../modules/compute_snp_signals.nf'
include { fit_nb_mixture }               from '../modules/fit_nb_mixture.nf'
include { rank_and_select }              from '../modules/rank_and_select.nf'
include { write_anchor_calc }            from '../modules/write_anchor_calc.nf'

// ---------------------------------------------------------------------------
// STUB PROCESSES — echo placeholders. Replaced in later chunks.
// Each prints what it WOULD do and emits correctly-typed placeholder files so
// the channel graph is exercised end to end.
// ---------------------------------------------------------------------------

// (chunk 3a) stub_build_candidate REMOVED — real candidates now come from
// filter_unique_seqs_candidate + assemble_rainbow_candidate.

// (chunk 3b) stub_cheap_signals and stub_provisional_rank REMOVED — real signals
// now come from compute_cheap_signals.nf + fit_nb_mixture.nf + rank_and_select.nf.
// (chunk 5b) stub_stage2 + stub_aggregate REMOVED — 1-pass design: compute_snp_signals
// runs on ALL candidates, rank_and_select does the single weight-free aggregation.

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
        // Grid is floor..ceiling on each cutoff (EXPLICIT ceilings — NOT the
        // auto-detect knee, which capped c1 at 5 / c2 at 4 and made the NB signal
        // degenerate by truncating the grid at the NB value itself). Extending the
        // ceiling above the expected optimum lets NB/quality signals bracket it in
        // the interior. The diagnostics knee is still computed (for the report and
        // as the NB fallback) but no longer bounds the optimization grid.
        def c1_floor   = (params.cutoff1_floor ?: 2) as int
        def c2_floor   = (params.cutoff2_floor ?: 3) as int
        def c1_ceiling = (params.cutoff1_ceiling ?: 8) as int
        def c2_ceiling = (params.cutoff2_ceiling ?: 6) as int

        c1_vals = (params.cutoff1 != null)
            ? Channel.value( [ params.cutoff1 as int ] )
            : Channel.value( (c1_floor..Math.max(c1_floor, c1_ceiling)).toList() )
        c2_vals = (params.cutoff2 != null)
            ? Channel.value( [ params.cutoff2 as int ] )
            : Channel.value( (c2_floor..Math.max(c2_floor, c2_ceiling)).toList() )

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
        // Signal set (chunk 3c):
        //   NB cutoff1 (5b)         : ranks the c1 axis (coverage-cutoff authority)
        //   coverage CV + Gini (2)  : rank assembly QUALITY along c2/similarity
        //   anchor (2, optional)    : ranks contig count vs expected RAD loci
        //   inflection (1a)         : computed but REPORT-ONLY (excluded from agg;
        //                             it duplicates NB on the c1 axis)
        // Dropped: redundancy (inherently ~0 for CD-HIT'd refs), paired-locus
        // fraction (all contigs paired by construction), length-dist (read-length
        // governed -> flat). See BUILD_STATUS chunk 3c rationale.

        // contig stats for inflection (report) + n_contigs for anchor
        compute_cheap_signals( candidates )
        cheap_rows = compute_cheap_signals.out.cheap_row.collect()

        // (5b) coverage-uniformity signal REMOVED — every coverage statistic we
        // tried (idxstats CV/Gini, per-base frac_hi/Gini/CV) stayed correlated with
        // assembly size on real data (frac_hi corr n_contigs = -0.80), i.e. it read
        // pruning-aggressiveness not an independent quality axis. The genotype
        // signals (concordance + r80) ARE size-orthogonal, so they carry quality now.

        // 5b: global NB-mixture cutoff1 from the pooled coverage distribution
        // (assembly_diagnostics now emits coverage_freq.txt; knee value is fallback).
        fit_nb_mixture(
            assembly_diagnostics.out.coverage_freq,
            assembly_diagnostics.out.cutoff1_value
        )

        // 2 anchor: expected ddRAD locus count. The CALCULATION now lives entirely
        // in write_anchor_calc (awk), which is the single source of truth and emits
        // both expected_loci.value (consumed by the rank step) and a published
        // human-readable breakdown. Here we only decide enabled vs disabled and pass
        // raw params through. IMPORTANT: insert window = genomic DNA between cut
        // sites, NOT the fragment-analyzer trace (which includes adapters).
        def anchor_enabled = (params.genome_size_est && params.size_select_min && params.size_select_max) as boolean
        if (anchor_enabled) {
            log.info "Biological anchor ENABLED: computing expected ddRAD loci (genome ${params.genome_size_est} bp, enzymes ${params.enzyme1_site_len}/${params.enzyme2_site_len}-cutter, insert ${params.size_select_min}-${params.size_select_max} bp)"
        } else {
            log.info "Biological anchor (signal 2) DISABLED — supply genome_size_est + size_select_min/max (INSERT window, not trace) to enable"
        }

        write_anchor_calc(
            anchor_enabled,
            params.genome_size_est ?: '',
            params.enzyme1_site_len ?: 6,
            params.enzyme2_site_len ?: 4,
            params.size_select_min ?: '',
            params.size_select_max ?: '',
            params.enzyme_pair
        )

        // expected_loci flows to the rank step as a value read from the process file
        // ("NA" when disabled; rank R drops the anchor signal on "NA").
        expected_loci_ch = write_anchor_calc.out.expected_loci.map { f -> f.text.trim() }

        // ----- PSEUDO-REPLICATES: N highest-depth individuals, split 50/50 -----
        // Used ONLY for the concordance signal (PASS A in compute_snp_signals).
        // STUB selection: first-N by SORTED sample id (depth-based selection is a
        // later refinement). Sort-before-take for cache stability.
        n_reps = params.n_pseudo_reps ?: 6
        concordance_indivs = cleaned_reads
            .toList()
            .flatMap { rows -> new ArrayList(rows).sort { a, b -> a[0] <=> b[0] }.take( n_reps as int ) }
        split_pseudo_rep( concordance_indivs, params.optimize_seed ?: 42 )

        // Flatten pseudo-rep halves into a collected bag of fastqs (every candidate
        // gets the same bag). Each emission: (parent, a_id, a_r1, a_r2, b_id, b_r1, b_r2)
        pseudo_fastqs = split_pseudo_rep.out.pseudo_reps
            .flatMap{ parent, aid, ar1, ar2, bid, br1, br2 -> [ar1, ar2, br1, br2] }
            .collect()

        // ----- SNP-subset: seeded snp_sample_pct of samples for r80 + density -----
        // sorted-before-take for cache stability; overlap with pseudo-reps is fine.
        def snp_pct = (params.snp_sample_pct ?: 25) as int
        snp_subset_reads = cleaned_reads
            .toList()
            .flatMap { rows ->
                def sorted = new ArrayList(rows).sort { a, b -> a[0] <=> b[0] }
                int n = Math.max(2, (int) Math.ceil(sorted.size() * snp_pct / 100.0))
                sorted.take(n)
            }
            .flatMap { sid, r1, r2 -> [r1, r2] }
            .collect()

        // ----- 1-PASS SNP SIGNALS on ALL candidates (no gate) -----
        // Each candidate gets the pseudo-rep bag (PASS A concordance) + the snp
        // subset bag (PASS B r80/density). Wrap each collected bag as one slot.
        snp_in = candidates
            .combine( pseudo_fastqs.map{ [it] } )
            .combine( snp_subset_reads.map{ [it] } )
            .map { meta, fa, preads, sreads -> tuple(meta, fa, preads, sreads) }

        compute_snp_signals( snp_in, (params.r80_threshold ?: 0.8) )
        snp_rows = compute_snp_signals.out.snp_row.collect()

        // ----- SINGLE rank aggregation (no provisional/final split) -----
        // Signals: NB cutoff1 proximity + anchor proximity + concordance (primary)
        // + r80 loci + SNPs-per-locus. Inflection report-only. Emits final_rank.tsv
        // + best_id.value (the winner). Weight-free rank aggregation.
        rank_and_select( cheap_rows, snp_rows, fit_nb_mixture.out.nb_cutoff1, expected_loci_ch )

        // ----- FINALIZE: select winning candidate as the production reference -----
        best_id_ch = rank_and_select.out.best_id.map { f -> f.text.trim() }

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
        optimize_summary = rank_and_select.out.final_rank
        optimize_plot    = rank_and_select.out.plot
}
