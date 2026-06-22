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
    echo "Selected winning de novo reference: ${meta.id} (cutoff1=${meta.c1} cutoff2=${meta.c2} init_sim=${meta.isim} div_f=${meta.divf} merge_r=${meta.mr} final_sim=${meta.fsim})"
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

        // Diagnostics: knee curves (report) + coverage_freq for the NB fit
        assembly_diagnostics( all_uniq_seqs )

        // NB-mixture cutoff1 must be computed BEFORE grid construction, because the
        // cutoff1 axis defaults to the NB crossover (params.cutoff1 == null). Knee
        // value is the graceful fallback inside fit_nb_mixture.
        fit_nb_mixture(
            assembly_diagnostics.out.coverage_freq,
            assembly_diagnostics.out.cutoff1_value
        )
        // nb value as an int channel (single value, broadcast)
        nb_c1_ch = fit_nb_mixture.out.nb_cutoff1.map { f -> f.text.trim() as int }

        // N = number of individuals entering assembly. Counted from the per-sample
        // uniq.seqs outputs (one per sample), NOT by re-tapping cleaned_reads. Used to
        // resolve the N-RELATIVE cutoff2 percentages into integer individual-counts.
        n_samples_ch = extract_unique_seqs.out.uniq_seqs.count()

        // ===== GRID MODEL =====================================================
        // Six axes, each a SCALAR (fixed) or a LIST (swept). Normalize each to a
        // list; the candidate grid is the Cartesian product. The cutoff axes are a
        // monotone stringency dial (no interior optimum), so by default they are
        // PINNED (c1=NB crossover, c2=3) and the assembly axes (final_sim, div_f,
        // optionally merge_r/init_sim) are SWEPT — that is where the STACKS-M-style
        // collapse-vs-split tradeoff lives. Selection uses the r80-vs-n_contigs
        // ELBOW (diminishing returns), not any single input axis, so this works
        // regardless of which axes the user chose to sweep.
        // Normalize an axis to a list. Accepts: a real Groovy list (from a config
        // or -params-file), a scalar number, or a STRING from the CLI. Nextflow does
        // NOT parse "--cutoff1 [3,5,7]" into a list — it arrives as the string
        // "[3,5,7]". So we also split comma/bracket/space strings here, making
        // --cutoff1 "3,5,7", --cutoff1 "[3,5,7]", and --cutoff1 5 all work.
        def asList = { v ->
            if (v == null) return [null]
            if (v instanceof List) return v
            if (v instanceof CharSequence) {
                def s = v.toString().trim().replaceAll(/^\[|\]$/, '').trim()
                if (s == '') return [null]
                if (s.contains(',') || s.contains(' '))
                    return s.split(/[,\s]+/).findAll { it }.collect { it.trim() }
                return [s]   // single value as string; coerced below
            }
            return [v]
        }

        // c1: null sentinel means "use the NB crossover" (resolved per-value below).
        def c1_axis   = asList(params.cutoff1)
        // c2: resolved at runtime. If params.cutoff2 (absolute) is set it OVERRIDES;
        // otherwise cutoff2 = round(N * pct/100) for each pct in params.cutoff2_pct,
        // floored at params.cutoff2_min and de-duplicated. (Resolved in the flatMap
        // below, where the runtime N value is available.)
        def c2_abs_override = (params.cutoff2 != null) ? asList(params.cutoff2).collect { it as int } : null
        def c2_pct_axis     = asList(params.cutoff2_pct).collect { it as double }
        def c2_min          = (params.cutoff2_min ?: 2) as int
        def isim_axis = asList(params.cluster_similarity).collect { it as double }
        def fsim_axis = asList(params.final_similarity).collect { it as double }
        def divf_axis = asList(params.div_f).collect { it as double }
        def mr_axis   = asList(params.merge_r).collect { it as int }

        // Resolve the c1 axis against the runtime NB value AND the c2 axis against the
        // runtime N, then build the full grid. Combine the two single-value channels
        // so both runtime numbers are in scope when the Cartesian product is built.
        grid_metas = nb_c1_ch.combine( n_samples_ch ).flatMap { nb, n_samples ->
            def c1_resolved = c1_axis.collect { it == null ? nb : (it as int) }.unique()
            // N-relative cutoff2 (or absolute override)
            def c2_axis = (c2_abs_override != null)
                ? c2_abs_override.unique()
                : c2_pct_axis.collect { pct -> Math.max(c2_min, Math.round(n_samples * pct / 100.0) as int) }.unique()
            log.info "cutoff2 resolved (N=${n_samples}): ${ c2_abs_override != null ? "absolute ${c2_axis}" : "${c2_pct_axis}% -> ${c2_axis} individuals" }"
            def combos = []
            for (c1 in c1_resolved)
              for (c2 in c2_axis)
                for (isim in isim_axis)
                  for (divf in divf_axis)
                    for (mr in mr_axis)
                      for (fsim in fsim_axis) {
                        def id = "c${c1}_k${c2}_is${isim}_f${divf}_r${mr}_fs${fsim}".toString()
                        combos << [ c1:c1, c2:c2, isim:isim, divf:divf, mr:mr, fsim:fsim,
                                    id_cutoff:"c${c1}_k${c2}".toString(), id:id ]
                      }
            // Grid-size guardrail
            int gmax = (params.max_grid_candidates ?: 64) as int
            if (combos.size() > gmax)
                log.warn "De novo grid = ${combos.size()} candidates (> max_grid_candidates=${gmax}). Each candidate is a full assembly + SNP pass; this may be very expensive. Pin more axes to scalars to shrink the grid."
            else
                log.info "De novo grid = ${combos.size()} candidates (c1=${c1_resolved} c2=${c2_axis} init_sim=${isim_axis} div_f=${divf_axis} merge_r=${mr_axis} final_sim=${fsim_axis})"
            combos
        }

        // ----- FILTER: depends only on (c1,c2); run ONCE per distinct cutoff pair -----
        // Collapse the grid to its distinct (c1,c2) pairs for filtering, then
        // re-expand by joining each filtered result back to every full meta that
        // shares its (c1,c2). id_cutoff is the join key.
        cutoff_pairs = grid_metas
            .map { m -> tuple(m.id_cutoff, [c1:m.c1, c2:m.c2]) }
            .unique { it[0] }

        filter_in = cutoff_pairs
            .combine( all_uniq_seqs.map{ [it] } )
            .map { idc, cm, files -> tuple(cm + [id_cutoff: idc], files) }

        filter_unique_seqs_candidate( filter_in )

        // Join filtered (c1,c2) outputs back to the full grid metas by id_cutoff.
        filtered_by_cutoff = filter_unique_seqs_candidate.out.filtered
            .map { meta, uniq_fasta, totaluniqseq, fstats ->
                tuple(meta.id_cutoff, uniq_fasta, totaluniqseq)
            }

        // Re-expand: each full grid meta joins to its (c1,c2)'s single filtered
        // result. combine(by:0) is the RIGHT operator here (NOT join): the left
        // side has DUPLICATE id_cutoff keys (many metas share a cutoff pair), which
        // join does not support, while the right side has exactly ONE item per key.
        // combine(by:0) does the within-key product = each meta x its 1 filtered row.
        // Ordering non-determinism is harmless: 1 right item per key means each meta
        // gets exactly its filtered fasta, and candidates are keyed by meta.id after.
        assemble_in = grid_metas
            .map { m -> tuple(m.id_cutoff, m) }
            .combine( filtered_by_cutoff, by: 0 )
            .map { idc, m, uniq_fasta, totaluniqseq -> tuple(m, uniq_fasta, totaluniqseq) }

        // assemble_rainbow_candidate now reads ALL assembly params from meta
        // (init_sim, div_f, merge_r, final_sim); only div_K stays a global param.
        assemble_rainbow_candidate( assemble_in, params.div_K )

        candidates = assemble_rainbow_candidate.out.reference   // tuple(meta, denovo_reference.fa)

        // ----- contig stats per candidate (n_contigs for the r80-vs-size curve + anchor) -----
        compute_cheap_signals( candidates )
        cheap_rows = compute_cheap_signals.out.cheap_row.collect()

        // (5b) coverage-uniformity signal REMOVED — every coverage statistic tried
        // stayed correlated with assembly size on real data. Selection is now the
        // r80-vs-n_contigs elbow (see rank_and_select.R).

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
        // Default 75% — r80 stability scales with sample count (the 80%-shared set
        // is better estimated with more samples).
        def snp_pct = (params.snp_sample_pct ?: 75) as int
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
        // Each candidate gets the pseudo-rep bag (PASS A concordance, REPORTED) + the
        // snp subset bag (PASS B r80/density). Wrap each collected bag as one slot.
        snp_in = candidates
            .combine( pseudo_fastqs.map{ [it] } )
            .combine( snp_subset_reads.map{ [it] } )
            .map { meta, fa, preads, sreads -> tuple(meta, fa, preads, sreads) }

        compute_snp_signals( snp_in, (params.r80_threshold ?: 0.8) )
        snp_rows = compute_snp_signals.out.snp_row.collect()

        // ----- SELECT via r80-vs-n_contigs ELBOW -----
        // PRIMARY selection = elbow of the r80(n_contigs) curve: the smallest
        // assembly past which marginal r80 gain per added contig falls below
        // r80_elbow_min_slope (diminishing returns -> only junk/paralogs added).
        // This is parameter-agnostic (works for any swept axes). concordance/anchor/
        // NB are REPORTED for context but do NOT drive selection (concordance is
        // size-monotone; the elbow on r80 is the real, size-aware optimum).
        rank_and_select(
            cheap_rows,
            snp_rows,
            fit_nb_mixture.out.nb_cutoff1,
            expected_loci_ch,
            (params.r80_elbow_min_slope ?: 30)
        )

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
