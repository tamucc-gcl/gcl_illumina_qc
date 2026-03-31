// modules/samtools_summary.nf
process samtools_summary {
    label 'samtools_summary'
    tag "mapping_summary"

    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode, pattern: "insert_size_violin.png"
    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode, pattern: "soft_clipping_violin.png"
    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode, pattern: "alignment_score_violin.png"
    
    input:
        path(stats_files)
        path(soft_clipping_files)
        path(alignment_score_files)
    
    output:
        path("mapping_summary.txt")
        path("insert_size_violin.png")
        path("soft_clipping_violin.png")
        path("alignment_score_violin.png")
    
    script:
    """
    # Create header
    echo -e "Sample\tTotal_Reads\tMapped_Reads\tMapped_Paired\tProperly_Paired\tMapping_Rate" > mapping_summary.txt
    
    # Process each stats file
    for stats_file in *.stats; do
        if [ -f "\$stats_file" ]; then
            sample=\$(basename "\$stats_file" .stats)
            
            total_reads=\$(grep "^SN.*raw total sequences:" "\$stats_file" | cut -f3)
            mapped_reads=\$(grep "^SN.*reads mapped:" "\$stats_file" | head -1 | cut -f3)
            mapped_paired=\$(grep "^SN.*reads mapped and paired:" "\$stats_file" | cut -f3)
            properly_paired=\$(grep "^SN.*reads properly paired:" "\$stats_file" | cut -f3)
            
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
            
            total_reads=\${total_reads:-0}
            mapped_reads=\${mapped_reads:-0}
            mapped_paired=\${mapped_paired:-0}
            properly_paired=\${properly_paired:-0}
            
            if [ "\$total_reads" -gt 0 ] 2>/dev/null; then
                mapping_rate=\$(echo "scale=2; \$mapped_reads * 100 / \$total_reads" | bc -l)
            else
                mapping_rate="0.00"
            fi
            
            echo -e "\$sample\t\$total_reads\t\$mapped_reads\t\$mapped_paired\t\$properly_paired\t\$mapping_rate" >> mapping_summary.txt
            
            echo "Processed \$sample: \$mapped_reads/\$total_reads mapped (\$mapping_rate%)"
        fi
    done
    
    echo ""
    echo "=== Final Mapping Summary ==="
    cat mapping_summary.txt

    # Make insert size violin plot using R script
    Rscript ${projectDir}/r_scripts/insert_violin.R

    # Make soft clipping and alignment score plots using R script
    Rscript ${projectDir}/r_scripts/bam_diagnostics.R
    """
}