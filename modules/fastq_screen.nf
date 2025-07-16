process fastq_screen {
    label 'fastq_screen'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(read), val(read_num)
        path(config_file)

    output:
        tuple val(sample_id), 
              path("${sample_id}_screen_${read_num}.fq.gz"), 
              path("${sample_id}_R${read_num}_screen.txt"),
              val(read_num)

    script:
    """
    # Count the number of databases in the config file
    # This counts lines that start with "DATABASE" (case insensitive)
    DB_COUNT=\$(grep -i "^DATABASE" ${config_file} | wc -l)
    
    # Generate filter string with the correct number of zeros
    FILTER_STRING=\$(printf "%0\${DB_COUNT}d" 0)
    
    echo "Found \$DB_COUNT databases in config file"
    echo "Using filter string: \$FILTER_STRING"
    
    echo "=== Before fastq_screen ==="
    ls -la
    
    fastq_screen \\
        --aligner bowtie2 \\
        --conf ${config_file} \\
        --threads ${task.cpus ?: 4} \\
        --tag \\
        --force \\
        --filter \$FILTER_STRING \\
        --subset 0 \\
        ${read}
    
    echo "=== After fastq_screen ==="
    ls -la
    
    echo "=== Looking for specific patterns ==="
    echo "Files matching *.tagged.fq.gz:"
    ls -la *.tagged.fq.gz 2>/dev/null || echo "No files found matching *.tagged.fq.gz"
    
    echo "Files matching *.tagged.fastq.gz:"
    ls -la *.tagged.fastq.gz 2>/dev/null || echo "No files found matching *.tagged.fastq.gz"
    
    echo "Files matching *_screen.txt:"
    ls -la *_screen.txt 2>/dev/null || echo "No files found matching *_screen.txt"
    
    echo "Files matching *.txt:"
    ls -la *.txt 2>/dev/null || echo "No files found matching *.txt"
    
    echo "All files in directory:"
    find . -type f -name "*" | sort
    
    echo "=== Attempting to rename files ==="
    
    # Try to find and rename the tagged file
    if ls *.tagged.fq.gz 1> /dev/null 2>&1; then
        echo "Found .tagged.fq.gz files"
        mv *.tagged.fq.gz ${sample_id}_screen_${read_num}.fq.gz
    elif ls *.tagged.fastq.gz 1> /dev/null 2>&1; then
        echo "Found .tagged.fastq.gz files"
        mv *.tagged.fastq.gz ${sample_id}_screen_${read_num}.fq.gz
    else
        echo "ERROR: Could not find expected tagged output file"
        echo "Looking for any files with 'tagged' in the name:"
        find . -name "*tagged*" -type f
        exit 1
    fi
    
    # Try to find and rename the screen report file
    if ls *_screen.txt 1> /dev/null 2>&1; then
        echo "Found _screen.txt files"
        mv *_screen.txt ${sample_id}_R${read_num}_screen.txt
    else
        echo "ERROR: Could not find expected screen report file"
        echo "Looking for any .txt files:"
        find . -name "*.txt" -type f
        exit 1
    fi
    
    echo "=== Final file check ==="
    ls -la
    echo "FastQ Screen step completed for ${sample_id} R${read_num}"
    """
}