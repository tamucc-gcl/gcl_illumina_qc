// modules/summarize_species_id.nf
process summarize_species_id {
    label 'species_summary'
    tag "species_identification_summary"
    
    publishDir "${params.outdir}/species_id", mode: 'copy'
    
    input:
        path(blast_results)  // Collected blast result files from all samples
    
    output:
        path("combined_blast_results.txt"), emit: combined_blast
        path("species_identification_report.txt"), emit: report
        path("species_consensus.txt"), emit: consensus
        path("species_identification_stats.txt"), emit: stats
        path("species_identification_plot.png"), emit: plot optional true
        path("species_by_sample.txt"), emit: by_sample optional true
    
    script:
    """
    # Combine all BLAST results into a single file with sample ID column
    echo "Combining BLAST results from all samples..."
    
    # Create header for combined file
    echo -e "sample_id\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxids\tsscinames" > combined_blast_results.txt
    
    # Process each BLAST result file
    for blast_file in *.blast_results.txt; do
        if [ -f "\$blast_file" ]; then
            # Extract sample name from filename (remove .blast_results.txt extension)
            sample=\$(basename "\$blast_file" .blast_results.txt)
            
            # Add sample ID as first column to each line and append to combined file
            awk -v sample="\$sample" 'BEGIN {OFS="\\t"} {print sample, \$0}' "\$blast_file" >> combined_blast_results.txt
            
            echo "  Added results from sample: \$sample"
        fi
    done
    
    # Count total hits
    total_hits=\$(tail -n +2 combined_blast_results.txt | wc -l)
    echo "Total BLAST hits across all samples: \$total_hits"
    
    # Get unique sample count
    n_samples=\$(tail -n +2 combined_blast_results.txt | cut -f1 | sort -u | wc -l)
    echo "Number of samples with BLAST results: \$n_samples"
    
    # R script for species identification analysis
    cat > analyze_species_id.R <<'REOF'
    #!/usr/bin/env Rscript
    
    # Load libraries
    suppressPackageStartupMessages({
        library(tidyverse)
    })
    
    # Read combined BLAST results
    cat("Reading combined BLAST results...\\n")
    blast_data <- read_delim("combined_blast_results.txt", 
                             delim = "\\t", 
                             show_col_types = FALSE,
                             col_names = c("sample_id", "qseqid", "sseqid", "pident", 
                                         "length", "mismatch", "gapopen", "qstart", 
                                         "qend", "sstart", "send", "evalue", 
                                         "bitscore", "staxids", "sscinames"),
                             skip = 1)
    
    cat(sprintf("Loaded %d BLAST hits from %d samples\\n", 
                nrow(blast_data), 
                n_distinct(blast_data\$sample_id)))
    
    # STUB IMPLEMENTATION - Replace with your actual R analysis
    # This is where you'll add your species identification logic
    
    # Example analysis structure:
    # 1. Filter for best hits per query per sample
    best_hits <- blast_data %>%
        group_by(sample_id, qseqid) %>%
        slice_max(bitscore, n = 1, with_ties = FALSE) %>%
        ungroup()
    
    # 2. Get top species per sample
    top_species_per_sample <- best_hits %>%
        group_by(sample_id) %>%
        count(sscinames, sort = TRUE) %>%
        slice_max(n, n = 3) %>%
        ungroup()
    
    # 3. Overall species consensus
    species_consensus <- best_hits %>%
        count(sscinames, sort = TRUE) %>%
        mutate(percent = n / sum(n) * 100)
    
    # Write outputs (placeholder for now)
    write_lines("Species identification analysis completed (R stub)", 
                "species_identification_report.txt")
    
    write_lines("Consensus species: To be implemented", 
                "species_consensus.txt")
    
    # Write sample-by-sample results
    top_species_per_sample %>%
        write_delim("species_by_sample.txt", delim = "\\t")
    
    # Statistics
    stats_text <- sprintf(
        "Total samples: %d\\nTotal BLAST hits: %d\\nUnique species identified: %d",
        n_distinct(blast_data\$sample_id),
        nrow(blast_data),
        n_distinct(blast_data\$sscinames)
    )
    write_lines(stats_text, "species_identification_stats.txt")
    
    # Placeholder for plot generation
    # library(ggplot2)
    # p <- ggplot(...) + ...
    # ggsave("species_identification_plot.png", p, width = 10, height = 8, dpi = 300)
    
    cat("R analysis complete\\n")
REOF
    
    # Run R script (comment out if you want to run manually for debugging)
    # Rscript analyze_species_id.R
    
    # For now, create stub outputs if R script hasn't run
    if [ ! -f "species_identification_report.txt" ]; then
        echo "==== Species Identification Summary Report ====" > species_identification_report.txt
        echo "" >> species_identification_report.txt
        echo "Analysis Date: \$(date)" >> species_identification_report.txt
        echo "Combined BLAST results saved to: combined_blast_results.txt" >> species_identification_report.txt
        echo "Total BLAST hits: \$total_hits" >> species_identification_report.txt
        echo "Number of samples: \$n_samples" >> species_identification_report.txt
        echo "" >> species_identification_report.txt
        echo "R analysis script created: analyze_species_id.R" >> species_identification_report.txt
        echo "Uncomment the Rscript line to run the analysis" >> species_identification_report.txt
    fi
    
    if [ ! -f "species_consensus.txt" ]; then
        cat > species_consensus.txt <<-EOF
		Consensus Species Identification
		=================================
		Combined BLAST results: combined_blast_results.txt
		Total hits: \$total_hits
		Samples analyzed: \$n_samples
		
		Run the R script for detailed analysis:
		  Rscript analyze_species_id.R
		EOF
    fi
    
    if [ ! -f "species_identification_stats.txt" ]; then
        cat > species_identification_stats.txt <<-EOF
		Species Identification Statistics
		==================================
		Total Samples: \$n_samples
		Total BLAST Hits: \$total_hits
		
		Detailed statistics will be generated by R script
		EOF
    fi
    
    echo ""
    echo "Species identification summary completed"
    echo "Files generated:"
    echo "  - combined_blast_results.txt (input for R analysis)"
    echo "  - analyze_species_id.R (R script for analysis)"
    echo "  - species_identification_report.txt"
    echo "  - species_consensus.txt"
    echo "  - species_identification_stats.txt"
    """
}