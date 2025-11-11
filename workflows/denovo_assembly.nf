// workflows/denovo_assembly.nf
// de novo assemble genome from ddrad: https://ddocent.com/assembly/

nextflow.enable.dsl = 2

// Import necessary modules
include { extract_unique_seqs } from '../modules/extract_unique_seqs.nf'
include { filter_unique_seqs } from '../modules/filter_unique_seqs.nf'
include { assemble_rainbow } from '../modules/assemble_rainbow.nf'

workflow denovo_assembly {
    take:
        cleaned_reads
    
    main:
        // Step 1: Extract unique sequences for each sample
        extract_unique_seqs( cleaned_reads )
        
        // Step 2: Collect all unique sequence files and filter
        all_uniq_seqs = extract_unique_seqs.out.uniq_seqs
            .map{ sid, file -> file }
            .collect()
        
        filter_unique_seqs( 
            all_uniq_seqs,
            params.cutoff1,
            params.cutoff2
        )
        
        // Step 3: Perform Rainbow assembly
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
        reference = assemble_rainbow.out.reference
        assembly_stats = assemble_rainbow.out.stats
        filter_stats = filter_unique_seqs.out.stats
    
}