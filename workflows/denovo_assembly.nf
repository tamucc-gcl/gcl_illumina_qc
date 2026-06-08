// workflows/denovo_assembly.nf
// de novo assemble genome from ddrad: https://ddocent.com/assembly/

nextflow.enable.dsl = 2

// Import necessary modules
include { extract_unique_seqs }  from '../modules/extract_unique_seqs.nf'
include { assembly_diagnostics } from '../modules/assembly_diagnostics.nf'
include { filter_unique_seqs }   from '../modules/filter_unique_seqs.nf'
include { assemble_rainbow }     from '../modules/assemble_rainbow.nf'

// Sweep subworkflow (only used when params.do_cluster_sweep is true)
include { cluster_similarity_sweep } from './cluster_similarity_sweep.nf'

workflow denovo_assembly {
    take:
        cleaned_reads

    main:
        // Step 1: Extract unique sequences for each sample
        extract_unique_seqs( cleaned_reads )

        // Step 2: Collect all unique sequence files
        all_uniq_seqs = extract_unique_seqs.out.uniq_seqs
            .map{ sid, file -> file }
            .collect()

        // Step 2a: Diagnostics — always run, auto-selects cutoffs from the data
        assembly_diagnostics( all_uniq_seqs )

        // Step 2b: Resolve effective cutoffs (non-null param overrides auto)
        cutoff1_ch = (params.cutoff1 != null)
            ? Channel.value( params.cutoff1 as int )
            : assembly_diagnostics.out.cutoff1_value.map{ f -> f.text.trim() as int }

        cutoff2_ch = (params.cutoff2 != null)
            ? Channel.value( params.cutoff2 as int )
            : assembly_diagnostics.out.cutoff2_value.map{ f -> f.text.trim() as int }

        // Step 3: Filter unique sequences with resolved cutoffs
        filter_unique_seqs(
            all_uniq_seqs,
            cutoff1_ch,
            cutoff2_ch
        )

        // Step 4: Assembly — sweep cluster_similarity, or single run
        if (params.do_cluster_sweep) {
            log.info "Sweeping cluster_similarity over: ${params.sweep_cluster_similarity}"

            cluster_similarity_sweep(
                filter_unique_seqs.out.filtered_fasta.first(),
                filter_unique_seqs.out.totaluniqseq.first(),
                cleaned_reads
            )

            reference_ch      = cluster_similarity_sweep.out.reference
            assembly_stats_ch = cluster_similarity_sweep.out.assembly_stats
            sweep_summary_ch  = cluster_similarity_sweep.out.sweep_summary
            sweep_plot_ch     = cluster_similarity_sweep.out.sweep_plot

        } else {
            assemble_rainbow(
                filter_unique_seqs.out.filtered_fasta,
                filter_unique_seqs.out.totaluniqseq,
                params.cluster_similarity,
                params.div_f,
                params.div_K,
                params.merge_r,
                params.final_similarity
            )

            reference_ch      = assemble_rainbow.out.reference
            assembly_stats_ch = assemble_rainbow.out.stats
            sweep_summary_ch  = Channel.empty()
            sweep_plot_ch     = Channel.empty()
        }

    emit:
        reference      = reference_ch
        assembly_stats = assembly_stats_ch
        filter_stats   = filter_unique_seqs.out.stats
        // diagnostics outputs (for the final report)
        cutoff1_plot   = assembly_diagnostics.out.cutoff1_plot
        cutoff2_plot   = assembly_diagnostics.out.cutoff2_plot
        diag_summary   = assembly_diagnostics.out.summary
        // sweep outputs (empty channels when sweep disabled)
        sweep_summary  = sweep_summary_ch
        sweep_plot     = sweep_plot_ch
}
