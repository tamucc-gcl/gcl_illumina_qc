// workflows/denovo_assembly.nf
// de novo assemble genome from ddrad: https://ddocent.com/assembly/

nextflow.enable.dsl = 2

include { extract_unique_seqs }  from '../modules/extract_unique_seqs.nf'
include { assembly_diagnostics } from '../modules/assembly_diagnostics.nf'
include { filter_unique_seqs }   from '../modules/filter_unique_seqs.nf'
include { assemble_rainbow }     from '../modules/assemble_rainbow.nf'

// Joint cutoff + cluster_similarity sweep (only when params.do_sweep is true)
include { assembly_sweep } from './assembly_sweep.nf'

workflow denovo_assembly {
    take:
        cleaned_reads

    main:
        // Step 1: per-sample unique sequences
        extract_unique_seqs( cleaned_reads )
        all_uniq_seqs = extract_unique_seqs.out.uniq_seqs
            .map{ sid, file -> file }
            .collect()

        // Step 2: diagnostics — always run (plots + knee values).
        // In sweep mode the knees are the grid CEILINGS; otherwise the chosen values.
        assembly_diagnostics( all_uniq_seqs )

        if (params.do_sweep) {
            log.info "Sweeping cutoffs (floor ${params.cutoff_sweep_floor}..knee) x cluster_similarity ${params.sweep_cluster_similarity}; sample subset = ${params.sweep_sample_pct}%"

            assembly_sweep(
                all_uniq_seqs,
                assembly_diagnostics.out.cutoff1_value,
                assembly_diagnostics.out.cutoff2_value,
                cleaned_reads
            )

            reference_ch      = assembly_sweep.out.reference
            assembly_stats_ch = assembly_sweep.out.assembly_stats
            sweep_summary_ch  = assembly_sweep.out.sweep_summary
            sweep_plot_ch     = assembly_sweep.out.sweep_plot

        } else {
            // Single run: resolved cutoffs (explicit param, else knee) at params.cluster_similarity
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

            reference_ch      = assemble_rainbow.out.reference
            assembly_stats_ch = assemble_rainbow.out.stats
            sweep_summary_ch  = Channel.empty()
            sweep_plot_ch     = Channel.empty()
        }

    emit:
        reference      = reference_ch
        assembly_stats = assembly_stats_ch
        filter_stats   = (params.do_sweep ? Channel.empty() : filter_unique_seqs.out.stats)
        cutoff1_plot   = assembly_diagnostics.out.cutoff1_plot
        cutoff2_plot   = assembly_diagnostics.out.cutoff2_plot
        diag_summary   = assembly_diagnostics.out.summary
        sweep_summary  = sweep_summary_ch
        sweep_plot     = sweep_plot_ch
}
