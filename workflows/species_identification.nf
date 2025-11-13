// workflows/species_identification.nf
// Species identification subworkflow using mitochondrial gene extraction and BLAST

nextflow.enable.dsl = 2

// Import necessary modules
include { get_mito_genes } from '../modules/get_mito_genes.nf'
include { blast_mito_genes } from '../modules/blast_mito_genes.nf'
include { summarize_species_id } from '../modules/summarize_species_id.nf'

workflow species_identification {
    take:
        reads_ch  // Channel of tuple(sample_id, read1, read2)
    
    main:
        // Extract mitochondrial genes from reads
        get_mito_genes(
            reads_ch,
            Channel.value(file(params.mito_reference)),
            Channel.value(params.genetic_code)
        )
        
        // BLAST mitochondrial genes
        blast_mito_genes(
            get_mito_genes.out,
            Channel.value(params.blast_db),
            Channel.value(params.taxonomy_db)
        )
        
        // Collect all BLAST results for summary
        all_blast_results = blast_mito_genes.out.blast_results
            .map{ sid, blast_file -> blast_file }
            .collect()
        
        // Summarize species identification across all samples
        summarize_species_id(all_blast_results)
    
    emit:
        // Individual sample results
        blast_results = blast_mito_genes.out.blast_results
        
        // Summary outputs from summarize_species_id
        combined_blast = summarize_species_id.out.combined_blast
        raw_pie_chart = summarize_species_id.out.raw_pie_chart
        summary_pie_chart = summarize_species_id.out.summary_pie_chart
        top_hits = summarize_species_id.out.top_hits
        posteriors = summarize_species_id.out.posteriors
}