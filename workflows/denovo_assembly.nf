// workflows/denovo_assembly.nf
// de novo assemble genome from ddrad: https://ddocent.com/assembly/
//
// Branches on params.do_optimize:
//   false -> single assembly at resolved cutoffs (explicit params, else knee)
//   true  -> optimize_denovo subworkflow (multi-signal rank-aggregation selector)

nextflow.enable.dsl = 2

include { extract_unique_seqs }  from '../modules/extract_unique_seqs.nf'
include { assembly_diagnostics } from '../modules/assembly_diagnostics.nf'
include { fit_nb_mixture }       from '../modules/fit_nb_mixture.nf'
include { filter_unique_seqs }   from '../modules/filter_unique_seqs.nf'
include { assemble_rainbow }     from '../modules/assemble_rainbow.nf'
include { optimize_denovo }      from './optimize_denovo.nf'

workflow denovo_assembly {
    take:
        cleaned_reads

    main:
        // Axes can now be scalars OR lists (grid model). The single-assembly path
        // needs ONE value per axis, so coerce: list -> first element, scalar -> itself.
        // First value of an axis (single-assembly path needs ONE value). Handles
        // real lists, scalars, and CLI list-strings ("[3,5]" / "3,5") consistently
        // with the optimize-path normalizer.
        def firstOf = { v ->
            if (v == null) return null
            if (v instanceof List) return v.isEmpty() ? null : v[0]
            if (v instanceof CharSequence) {
                def s = v.toString().trim().replaceAll(/^\[|\]$/, '').trim()
                if (s == '') return null
                return (s.split(/[,\s]+/) as List).find { it } 
            }
            return v
        }

        if (params.do_optimize) {
            log.info "De novo OPTIMIZATION enabled: grid over c1=${params.cutoff1 ?: 'NB'} c2=${params.cutoff2} init_sim=${params.cluster_similarity} final_sim=${params.final_similarity} div_f=${params.div_f} merge_r=${params.merge_r}; r80 SNP subset ${params.snp_sample_pct}%, ${params.n_pseudo_reps} pseudo-replicates"

            optimize_denovo( cleaned_reads )

            reference_ch        = optimize_denovo.out.reference
            assembly_stats_ch   = optimize_denovo.out.assembly_stats
            cutoff1_plot_ch     = optimize_denovo.out.cutoff1_plot
            cutoff2_plot_ch     = optimize_denovo.out.cutoff2_plot
            diag_summary_ch     = optimize_denovo.out.diag_summary
            optimize_summary_ch = optimize_denovo.out.optimize_summary
            optimize_plot_ch    = optimize_denovo.out.optimize_plot
            filter_stats_ch     = Channel.empty()   // no single filter step in optimize mode

        } else {
            // ---- Single run: ONE value per axis (first element if a list) ----
            extract_unique_seqs( cleaned_reads )
            all_uniq_seqs = extract_unique_seqs.out.uniq_seqs
                .map{ sid, f -> f }
                .collect()

            assembly_diagnostics( all_uniq_seqs )

            // NB-mixture cutoff1 (same principled value used in optimize mode), with
            // the diagnostics knee as graceful fallback inside fit_nb_mixture.
            fit_nb_mixture(
                assembly_diagnostics.out.coverage_freq,
                assembly_diagnostics.out.cutoff1_value
            )

            def c1_fixed   = firstOf(params.cutoff1)
            def isim_fixed = firstOf(params.cluster_similarity)
            def divf_fixed = firstOf(params.div_f)
            def mr_fixed   = firstOf(params.merge_r)
            def fsim_fixed = firstOf(params.final_similarity)

            // cutoff1: explicit value if given (first element if a list); else the
            // NB crossover (consistent with optimize mode — NOT the geometric knee).
            cutoff1_ch = (c1_fixed != null)
                ? Channel.value( c1_fixed as int )
                : fit_nb_mixture.out.nb_cutoff1.map{ f -> f.text.trim() as int }

            // cutoff2: N-RELATIVE, same logic as optimize mode but resolved to ONE
            // value. Precedence: explicit absolute params.cutoff2 (first if a list)
            //   -> else round(N * firstOf(cutoff2_pct)/100) floored at cutoff2_min
            //   -> else (no pct supplied) the diagnostics knee.
            def c2_abs    = firstOf(params.cutoff2)
            def c2_pct    = firstOf(params.cutoff2_pct)
            def c2_min    = (params.cutoff2_min ?: 2) as int
            n_samples_ch  = extract_unique_seqs.out.uniq_seqs.count()
            cutoff2_ch =
                (c2_abs != null)
                    ? Channel.value( (c2_abs as Double) as int )
                    : (c2_pct != null)
                        ? n_samples_ch.map{ n -> Math.max(c2_min, Math.round(n * (c2_pct as double) / 100.0) as int) }
                        : assembly_diagnostics.out.cutoff2_value.map{ f -> f.text.trim() as int }

            filter_unique_seqs( all_uniq_seqs, cutoff1_ch, cutoff2_ch )

            assemble_rainbow(
                filter_unique_seqs.out.filtered_fasta,
                filter_unique_seqs.out.totaluniqseq,
                isim_fixed,
                divf_fixed, params.div_K, mr_fixed, fsim_fixed
            )

            reference_ch        = assemble_rainbow.out.reference
            assembly_stats_ch   = assemble_rainbow.out.stats
            cutoff1_plot_ch     = assembly_diagnostics.out.cutoff1_plot
            cutoff2_plot_ch     = assembly_diagnostics.out.cutoff2_plot
            diag_summary_ch     = assembly_diagnostics.out.summary
            optimize_summary_ch = Channel.empty()
            optimize_plot_ch    = Channel.empty()
            filter_stats_ch     = filter_unique_seqs.out.stats
        }

    emit:
        reference        = reference_ch
        assembly_stats   = assembly_stats_ch
        filter_stats     = filter_stats_ch
        cutoff1_plot     = cutoff1_plot_ch
        cutoff2_plot     = cutoff2_plot_ch
        diag_summary     = diag_summary_ch
        optimize_summary = optimize_summary_ch
        optimize_plot    = optimize_plot_ch
}
