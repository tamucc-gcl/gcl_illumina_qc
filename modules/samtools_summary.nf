// modules/samtools_summary.nf
process samtools_summary {
    label 'samtools_summary'
    tag "mapping_summary"
    
    publishDir "${params.outdir}/mapping_stats", mode: 'copy'
    
    input:
        path(stats_files)
    
    output:
        path("mapping_summary.txt")
    
    script:
    """
    # Create header
    echo -e "Sample\tTotal_Reads\tMapped_Reads\tMapped_Paired\tProperly_Paired\tMapping_Rate" > mapping_summary.txt
    
    # Process each stats file
    for stats_file in *.stats; do
        if [ -f "\$stats_file" ]; then
            # Extract sample name from filename (remove .stats extension)
            sample=\$(basename "\$stats_file" .stats)
            
            # Extract key statistics using grep and awk
            total_reads=\$(grep "^SN.*raw total sequences:" "\$stats_file" | awk '{print \$NF}')
            mapped_reads=\$(grep "^SN.*reads mapped:" "\$stats_file" | awk '{print \$NF}')
            mapped_paired=\$(grep "^SN.*reads mapped and paired:" "\$stats_file" | awk '{print \$NF}')
            properly_paired=\$(grep "^SN.*reads properly paired:" "\$stats_file" | awk '{print \$NF}')
            
            # Calculate mapping rate as percentage using awk
            if [ "\$total_reads" -gt 0 ]; then
                mapping_rate=\$(awk "BEGIN {printf \"%.2f\", \$mapped_reads * 100 / \$total_reads}")
            else
                mapping_rate="0.00"
            fi
            
            # Add to summary file
            echo -e "\$sample\t\$total_reads\t\$mapped_reads\t\$mapped_paired\t\$properly_paired\t\$mapping_rate" >> mapping_summary.txt
            
            echo "Processed \$sample: \$mapped_reads/\$total_reads mapped (\$mapping_rate%)"
        fi
    done
    
    echo "=== Final Mapping Summary ==="
    cat mapping_summary.txt
    """
}