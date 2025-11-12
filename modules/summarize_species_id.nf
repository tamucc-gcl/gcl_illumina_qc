// modules/summarize_species_id.nf
process summarize_species_id {
    label 'species_summary'
    tag "species_identification_summary"
    
    publishDir "${params.outdir}/species_id", mode: 'copy'
    
    input:
        path(blast_results)  // Collected blast result files from all samples
    
    output:
        path("combined_blast_results.tsv"), emit: combined_blast
    
    script:
    """
    # Combine all BLAST results into a single file with sample ID column
    echo "Combining BLAST results from all samples..."
    
    # Create header for combined file
    echo -e "sample_id\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxids\tsscinames\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies" > combined_blast_results.tsv

    # Process each BLAST result file (e.g., *.blast_with_taxa.header.tsv)
    # Input files are staged in working directory
    for blast_file in *.blast_results.txt; do
        if [ -f "\${blast_file}" ]; then
            # Extract sample name (everything before .blast)
            sample=\$(basename "\${blast_file}" | sed 's/\.blast.*//')

            # Skip header lines if present and prepend sample ID
            awk -v sample="\${sample}" -F'\t' 'BEGIN {OFS="\t"} NR>1 {print sample, \$0}' "\$blast_file" >> combined_blast_results.tsv

            echo "  Added results from sample: \$sample"
        fi
    done
    
    # Report summary
    total_hits=\$(tail -n +2 combined_blast_results.tsv | wc -l)
    n_samples=\$(tail -n +2 combined_blast_results.tsv | cut -f1 | sort -u | wc -l)
    
    echo ""
    echo "Summary:"
    echo "  Total BLAST hits: \$total_hits"
    echo "  Number of samples: \$n_samples"
    echo "  Output file: combined_blast_results.tsv"
    """
}