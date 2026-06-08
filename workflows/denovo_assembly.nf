// workflows/denovo_assembly.nf
// de novo assemble genome from ddrad: https://ddocent.com/assembly/

nextflow.enable.dsl = 2

// Import necessary modules
include { extract_unique_seqs }  from '../modules/extract_unique_seqs.nf'
include { assembly_diagnostics } from '../modules/assembly_diagnostics.nf'
include { filter_unique_seqs }   from '../modules/filter_unique_seqs.nf'
include { assemble_rainbow }     from '../modules/assemble_rainbow.nf'

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

        // Step 2a: Diagnostics — always run (produces report plots), and
        //          auto-selects cutoffs from the data.
        assembly_diagnostics( all_uniq_seqs )

        // Step 2b: Resolve effective cutoffs.
        //          A non-null params value overrides the auto-detected one.
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

        // Step 4: Perform Rainbow assembly
        assemble_rainbow(
            filter_unique_seqs.out.filtered_fasta,
            filter_unique_seqs.out.totaluniqseq,
            params.cluster_similarity,
            params.div_f,
            params.div_K,
            params.merge_r,
            params.final_similarity
        )

    emit:
        reference      = assemble_rainbow.out.reference
        assembly_stats = assemble_rainbow.out.stats
        filter_stats   = filter_unique_seqs.out.stats
        // diagnostics outputs (for the final report)
        cutoff1_plot   = assembly_diagnostics.out.cutoff1_plot
        cutoff2_plot   = assembly_diagnostics.out.cutoff2_plot
        diag_summary   = assembly_diagnostics.out.summary
}
