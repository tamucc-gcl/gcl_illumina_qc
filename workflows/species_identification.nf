// workflows/species_identification.nf
// Species identification subworkflow using mitochondrial gene extraction and BLAST

nextflow.enable.dsl = 2

// Import necessary modules
include { get_mito_genes } from '../modules/get_mito_genes.nf'
include { blast_mito_genes } from '../modules/blast_mito_genes.nf'
include { summarize_species_id } from '../modules/summarize_species_id.nf'

workflow species_identification {
    take:
        fastp_trim5_reads  // Channel of tuples: [sample_id, read1, read2] from fastp_trim_5
    
    main:
        // Log the start of species identification
        log.info "Starting species identification workflow"
        
        // Step 1: Extract mitochondrial genes from each sample
        // Process each sample independently
        get_mito_genes(
            fastp_trim5_reads,
            file(params.mito_reference),
            params.genetic_code
        )
        
        
        // Step 2: BLAST extracted mitochondrial genes
        // Create channel for BLAST database
        blast_db = Channel.value(params.blast_db)
        
        blast_mito_genes(
            get_mito_genes.out,
            blast_db
        )
        
        
        // Step 3: Collect all results and generate summary
        // Collect all BLAST results
        all_blast_results = blast_mito_genes.out.blast_results
            .map{ sid, file -> file }
            .collect()

        /*
        // Collect all mitochondrial gene files
        all_mito_genes = get_mito_genes.out
            .map{ sid, file -> file }
            .collect()
        */
        
        // Generate comprehensive species identification summary
        summarize_species_id(
            all_blast_results
        )
        
        // Log completion
        summarize_species_id.out.report.view { 
            "Species identification completed. Report saved: $it" 
        }
    
    emit:
        mito_genes = get_mito_genes.out
        blast_results = blast_mito_genes.out.blast_results
        /*
        species_report = summarize_species_id.out.report
        species_consensus = summarize_species_id.out.consensus
        species_stats = summarize_species_id.out.stats
        */
}