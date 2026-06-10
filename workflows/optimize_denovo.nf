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

// ---------------------------------------------------------------------------
// STUB PROCESSES — echo placeholders. Replaced in later chunks.
// Each prints what it WOULD do and emits correctly-typed placeholder files so
// the channel graph is exercised end to end.
// ---------------------------------------------------------------------------

// (chunk 3a) stub_build_candidate REMOVED — real candidates now come from
// filter_unique_seqs_candidate + assemble_rainbow_candidate.

// Stand-in for: cheap signals (1a inflection, 1b redundancy, 5b NB-mixture, 2 anchor)
// per candidate -> one provisional-score row. Real version = chunks 3.
process stub_cheap_signals {
    label 'basic'
    tag "${meta.id}"
    input:
        tuple val(meta), path(candidate)
    output:
        path("cheap_${meta.id}.tsv"), emit: cheap_row
    script:
    """
    echo "[STUB] cheap signals for ${meta.id}"
    # id c1 c2 sim n_contigs redundancy inflection nbfit anchor_dev
    printf "%s\\t%s\\t%s\\t%s\\t0\\t0\\t0\\t0\\tNA\\n" "${meta.id}" "${meta.c1}" "${meta.c2}" "${meta.sim}" > cheap_${meta.id}.tsv
    """
}

// Stand-in for: provisional rank + gap-based top-N selection. Real version = chunk 4.
process stub_provisional_rank {
    label 'basic'
    tag "provisional_rank"
    input:
        path(cheap_rows)
    output:
        path("survivors.txt"),        emit: survivors   // one meta.id per line
        path("provisional_rank.tsv"), emit: provisional
    script:
    """
    echo "[STUB] provisional rank over \$(ls cheap_*.tsv 2>/dev/null | wc -l) candidates"
    cat cheap_*.tsv > provisional_rank.tsv
    # STUB: 'select' the first up-to-3 candidate ids as survivors
    cut -f1 provisional_rank.tsv | head -n 3 > survivors.txt
    echo "[STUB] survivors:"; cat survivors.txt
    """
}

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

        // ----- STAGE 1: cheap signals per candidate -> provisional rank -----
        stub_cheap_signals( candidates )
        cheap_rows = stub_cheap_signals.out.cheap_row.collect()
        stub_provisional_rank( cheap_rows )
        survivors_ch = stub_provisional_rank.out.survivors


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
        stub_aggregate( stub_provisional_rank.out.provisional, stage2_rows )

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
