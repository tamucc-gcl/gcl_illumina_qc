// modules/multiqc.nf
process multiqc {
    label 'multiqc'
    tag "${step_name}"
    
    // Only publish the HTML report to publishDir
    publishDir "${params.outdir}/multiqc_reports", mode: 'copy', pattern: "*.html"

    input:
        path(input_files)
        val(step_name)

    output:
        path("multiqc_${step_name}.html")
        path("multiqc_${step_name}_general_stats.txt")

    script:
    """
    multiqc \\
        --title "MultiQC Report - ${step_name}" \\
        --filename multiqc_${step_name}.html \\
        --export \\
        --data-dir \\
        .
    
    # The general stats file is created in multiqc_data directory
    # Copy it to the main directory with a descriptive name
    if [ -f multiqc_data/multiqc_general_stats.txt ]; then
        cp multiqc_data/multiqc_general_stats.txt multiqc_${step_name}_general_stats.txt
    else
        # Create empty file if no stats were generated
        touch multiqc_${step_name}_general_stats.txt
    fi
    """
}