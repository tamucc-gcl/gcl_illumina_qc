// workflows/denovo_assembly.nf
// de novo assemble genome from ddrad: https://ddocent.com/assembly/
//
// Branches on params.do_optimize:
//   false -> single assembly at resolved cutoffs (explicit params, else knee)
//   true  -> optimize_denovo subworkflow (multi-signal rank-aggregation selector)

nextflow.enable.dsl = 2

include { extract_unique_seqs }  from '../modules/extract_unique_seqs.nf'
include { assembly_diagnostics } from '../modules/assembly_diagnostics.nf'
include { filter_unique_seqs }   from '../modules/filter_unique_seqs.nf'
include { assemble_rainbow }     from '../modules/assemble_rainbow.nf'
include { optimize_denovo }      from './optimize_denovo.nf'

workflow denovo_assembly {
    take:
        cleaned_reads

    main:
        if (params.do_optimize) {
            log.info "De novo OPTIMIZATION enabled: grid c1 ${params.cutoff1_floor}..knee x c2 ${params.cutoff2_floor}..knee x cluster_similarity ${params.optimize_cluster_similarity}; assembly subset ${params.optimize_sample_pct}%, ${params.n_pseudo_reps} pseudo-replicates"

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
            // ---- Single run: resolved cutoffs at params.cluster_similarity ----
            extract_unique_seqs( cleaned_reads )
            all_uniq_seqs = extract_unique_seqs.out.uniq_seqs
                .map{ sid, f -> f }
                .collect()

            assembly_diagnostics( all_uniq_seqs )

            cutoff1_ch = (params.cutoff1 != null)
                ? Channel.value( params.cutoff1 as int )
                : assembly_diagnostics.out.cutoff1_value.map{ f -> f.text.trim() as int }
            cutoff2_ch = (params.cutoff2 != null)
                ? Channel.value( params.cutoff2 as int )
                : assembly_diagnostics.out.cutoff2_value.map{ f -> f.text.trim() as int }

            filter_unique_seqs( all_uniq_seqs, cutoff1_ch, cutoff2_ch )

            assemble_rainbow(
                filter_unique_seqs.out.filtered_fasta,
                filter_unique_seqs.out.totaluniqseq,
                params.cluster_similarity,
                params.div_f, params.div_K, params.merge_r, params.final_similarity
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
