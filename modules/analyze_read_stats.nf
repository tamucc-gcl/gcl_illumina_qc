// modules/analyze_read_stats.nf
process analyze_read_stats {
    label 'r_analysis'
    tag "read_stats_exploration"
    
    publishDir "${params.outdir}/read_analysis", mode: 'copy'
    
    input:
        path(stats_files)
    
    output:
        path("file_listing.txt")
        path("file_info.txt")
        path("*.txt")  // Pass through all input files
    
    script:
    """
    echo "=== MultiQC Stats Files Found ===" > file_listing.txt
    echo "" >> file_listing.txt
    
    # List all files with their full paths
    for file in *.txt; do
        if [[ -f "\$file" ]]; then
            echo "File: \$file" >> file_listing.txt
            echo "Full path: \$(realpath \$file)" >> file_listing.txt
            echo "Size: \$(stat -c%s \$file) bytes" >> file_listing.txt
            echo "" >> file_listing.txt
        fi
    done
    
    echo "=== File Contents Preview ===" > file_info.txt
    echo "" >> file_info.txt
    
    # Show first few lines of each file
    for file in *.txt; do
        if [[ -f "\$file" && "\$file" != "file_listing.txt" && "\$file" != "file_info.txt" ]]; then
            echo "=== \$file ===" >> file_info.txt
            echo "First 10 lines:" >> file_info.txt
            head -10 "\$file" >> file_info.txt
            echo "" >> file_info.txt
            echo "Column headers:" >> file_info.txt
            head -1 "\$file" | tr '\t' '\n' | nl >> file_info.txt
            echo "" >> file_info.txt
            echo "---" >> file_info.txt
            echo "" >> file_info.txt
        fi
    done
    
    echo "Files ready for interactive R exploration:"
    echo "- file_listing.txt: List of all files with paths"
    echo "- file_info.txt: Preview of file contents and structure"
    echo "- Individual MultiQC stats files are also copied to output directory"
    """
}