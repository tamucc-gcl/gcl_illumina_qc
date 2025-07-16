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
    echo "=== MultiQC Debug for ${step_name} ==="
    echo "Input files:"
    ls -la
    
    echo "=== Running MultiQC ==="
    multiqc \\
        --title "MultiQC Report - ${step_name}" \\
        --filename multiqc_${step_name}.html \\
        --export \\
        --data-dir \\
        .
    
    echo "=== After MultiQC ==="
    echo "All files in current directory:"
    ls -la
    
    echo "=== Looking for multiqc_data directory ==="
    if [ -d multiqc_data ]; then
        echo "multiqc_data directory exists, contents:"
        ls -la multiqc_data/
    else
        echo "multiqc_data directory does not exist"
    fi
    
    echo "=== Looking for any general stats files ==="
    find . -name "*general_stats*" -type f
    
    echo "=== Looking for any multiqc_data files ==="
    find . -name "multiqc_data*" -type f
    
    echo "=== Trying to find and copy general stats file ==="
    if [ -f multiqc_data/multiqc_general_stats.txt ]; then
        echo "Found multiqc_data/multiqc_general_stats.txt"
        cp multiqc_data/multiqc_general_stats.txt multiqc_${step_name}_general_stats.txt
        echo "Copied to multiqc_${step_name}_general_stats.txt"
        echo "File size: \$(wc -l multiqc_${step_name}_general_stats.txt)"
    elif [ -f multiqc_general_stats.txt ]; then
        echo "Found multiqc_general_stats.txt in current directory"
        cp multiqc_general_stats.txt multiqc_${step_name}_general_stats.txt
        echo "Copied to multiqc_${step_name}_general_stats.txt"
        echo "File size: \$(wc -l multiqc_${step_name}_general_stats.txt)"
    else
        echo "ERROR: Could not find any general stats file"
        echo "Creating empty file"
        touch multiqc_${step_name}_general_stats.txt
    fi
    
    echo "=== Final directory contents ==="
    ls -la
    
    echo "=== Final stats file check ==="
    if [ -f multiqc_${step_name}_general_stats.txt ]; then
        echo "Stats file exists with size: \$(wc -l multiqc_${step_name}_general_stats.txt)"
        echo "First few lines:"
        head -5 multiqc_${step_name}_general_stats.txt
    else
        echo "Stats file does not exist!"
    fi
    """
}