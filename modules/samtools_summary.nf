// modules/samtools_summary.nf
process samtools_summary {
    label 'samtools_summary'
    tag "mapping_summary"
    
    publishDir "${params.outdir}/mapping_stats", mode: 'copy'
    
    input:
        path(stats_files)
    
    output:
        path("mapping_summary.txt")
        path("debug_stats.txt")
    
    script:
    """
    # Create header
    echo -e "Sample\tTotal_Reads\tMapped_Reads\tMapped_Paired\tProperly_Paired\tMapping_Rate" > mapping_summary.txt
    
    # Create debug file to see what we're parsing
    echo "=== DEBUG: Examining samtools stats format ===" > debug_stats.txt
    
    # Process each stats file
    for stats_file in *.stats; do
        if [ -f "\$stats_file" ]; then
            # Extract sample name from filename (remove .stats extension)
            sample=\$(basename "\$stats_file" .stats)
            
            echo "" >> debug_stats.txt
            echo "=== Processing \$sample ===" >> debug_stats.txt
            echo "Stats file: \$stats_file" >> debug_stats.txt
            echo "" >> debug_stats.txt
            
            # Show the relevant lines from the stats file for debugging
            echo "Raw total sequences line:" >> debug_stats.txt
            grep "^SN.*raw total sequences:" "\$stats_file" >> debug_stats.txt
            echo "" >> debug_stats.txt
            
            echo "Reads mapped line:" >> debug_stats.txt
            grep "^SN.*reads mapped:" "\$stats_file" | head -1 >> debug_stats.txt
            echo "" >> debug_stats.txt
            
            echo "Reads mapped and paired line:" >> debug_stats.txt
            grep "^SN.*reads mapped and paired:" "\$stats_file" >> debug_stats.txt
            echo "" >> debug_stats.txt
            
            echo "Reads properly paired line:" >> debug_stats.txt
            grep "^SN.*reads properly paired:" "\$stats_file" >> debug_stats.txt
            echo "" >> debug_stats.txt
            
            # Extract key statistics using corrected field extraction
            # samtools stats format: SN<tab>description:<tab>number
            total_reads=\$(grep "^SN.*raw total sequences:" "\$stats_file" | cut -f3)
            mapped_reads=\$(grep "^SN.*reads mapped:" "\$stats_file" | head -1 | cut -f3)
            mapped_paired=\$(grep "^SN.*reads mapped and paired:" "\$stats_file" | cut -f3)
            properly_paired=\$(grep "^SN.*reads properly paired:" "\$stats_file" | cut -f3)
            
            # Alternative extraction using awk with tab separator
            if [ -z "\$total_reads" ]; then
                total_reads=\$(grep "^SN.*raw total sequences:" "\$stats_file" | awk -F'\t' '{print \$3}')
            fi
            if [ -z "\$mapped_reads" ]; then
                mapped_reads=\$(grep "^SN.*reads mapped:" "\$stats_file" | head -1 | awk -F'\t' '{print \$3}')
            fi
            if [ -z "\$mapped_paired" ]; then
                mapped_paired=\$(grep "^SN.*reads mapped and paired:" "\$stats_file" | awk -F'\t' '{print \$3}')
            fi
            if [ -z "\$properly_paired" ]; then
                properly_paired=\$(grep "^SN.*reads properly paired:" "\$stats_file" | awk -F'\t' '{print \$3}')
            fi
            
            # Set defaults for empty values
            total_reads=\${total_reads:-0}
            mapped_reads=\${mapped_reads:-0}
            mapped_paired=\${mapped_paired:-0}
            properly_paired=\${properly_paired:-0}
            
            # Debug output
            echo "Extracted values:" >> debug_stats.txt
            echo "  total_reads: \$total_reads" >> debug_stats.txt
            echo "  mapped_reads: \$mapped_reads" >> debug_stats.txt
            echo "  mapped_paired: \$mapped_paired" >> debug_stats.txt
            echo "  properly_paired: \$properly_paired" >> debug_stats.txt
            echo "" >> debug_stats.txt
            
            # Calculate mapping rate as percentage
            if [ "\$total_reads" -gt 0 ] 2>/dev/null; then
                mapping_rate=\$(awk "BEGIN {printf \"%.2f\", \$mapped_reads * 100 / \$total_reads}")
            else
                mapping_rate="0.00"
            fi
            
            # Add to summary file
            echo -e "\$sample\t\$total_reads\t\$mapped_reads\t\$mapped_paired\t\$properly_paired\t\$mapping_rate" >> mapping_summary.txt
            
            echo "Processed \$sample: \$mapped_reads/\$total_reads mapped (\$mapping_rate%)"
        fi
    done
    
    echo ""
    echo "=== Final Mapping Summary ==="
    cat mapping_summary.txt
    echo ""
    echo "Debug information written to debug_stats.txt"
    """
}