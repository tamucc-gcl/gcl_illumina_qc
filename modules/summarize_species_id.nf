// modules/summarize_species_id.nf
process summarize_species_id {
    label 'species_summary'
    tag "species_identification_summary"
    
    publishDir "${params.outdir}/species_id", mode: 'copy', pattern: "*.{tsv,png,csv}"
    publishDir "${params.outdir}/species_id/blast_posteriors", mode: 'copy', pattern: "blast_posteriors/*.csv"
    
    input:
        path(blast_results)  // Collected blast result files from all samples
    
    output:
        path("blast_results.tsv"), emit: combined_blast
        path("blast_raw_pie.png"), emit: raw_pie_chart, optional: true
        path("blast_summary_pie.png"), emit: summary_pie_chart, optional: true
        path("top_blast_hits.csv"), emit: top_hits, optional: true
        path("blast_posteriors/*.csv"), emit: posteriors, optional: true
    
    script:
    """
    # Combine all BLAST results into a single file with sample ID column
    echo "Combining BLAST results from all samples..."
    
    # Create header for combined file
    echo -e "sample_id\\tqseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tstaxids\\tsscinames\\tkingdom\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies" > blast_results.tsv

    # Process each BLAST result file
    for blast_file in *.blast_results.txt; do
        if [ -f "\${blast_file}" ]; then
            # Extract sample name (everything before .blast_results.txt)
            sample=\$(basename "\${blast_file}" .blast_results.txt)

            # Skip header lines if present and prepend sample ID
            awk -v sample="\${sample}" -F'\\t' 'BEGIN {OFS="\\t"} NR>1 {print sample, \$0}' "\${blast_file}" >> blast_results.tsv

            echo "  Added results from sample: \${sample}"
        fi
    done
    
    # Check if any files were processed
    if [ ! -s blast_results.tsv ] || [ \$(wc -l < blast_results.tsv) -eq 1 ]; then
        echo "Warning: No BLAST results were found or processed"
        echo "Files in directory:"
        ls -la
        # Create empty output files to prevent pipeline failure
        touch blast_raw_pie.png blast_summary_pie.png top_blast_hits.csv
        mkdir -p blast_posteriors
        touch blast_posteriors/empty.csv
        exit 0
    fi
    
    # Report summary
    total_hits=\$(tail -n +2 blast_results.tsv | wc -l)
    n_samples=\$(tail -n +2 blast_results.tsv | cut -f1 | sort -u | wc -l)
    
    echo ""
    echo "Summary:"
    echo "  Total BLAST hits: \$total_hits"
    echo "  Number of samples: \$n_samples"
    echo "  Output file: blast_results.tsv"
    
    # Run R script for analysis
    echo ""
    echo "Running R analysis for species identification..."
    
    # Copy R script from module directory to working directory
    cp ${projectDir}/modules/blast_summary.R .
    
    # Run the R script
    Rscript blast_summary.R
    
    # Check if R script produced outputs
    if [ ! -f "top_blast_hits.csv" ]; then
        echo "Warning: R script did not produce expected outputs"
        # Create minimal outputs to prevent pipeline failure
        touch blast_raw_pie.png blast_summary_pie.png top_blast_hits.csv
        mkdir -p blast_posteriors
        touch blast_posteriors/empty.csv
    fi
    
    echo "Species identification analysis complete"
    """
}