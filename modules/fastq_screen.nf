process fastq_screen {
    label 'fastq_screen'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(read), val(read_num)
        path(config_file)

    output:
        tuple val(sample_id), 
              path("${sample_id}_fp1-clmp-fp2-fqscrn_r${read_num}.fq.gz"), 
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
    

    fastq_screen \\
        --aligner bowtie2 \\
        --conf ${config_file} \\
        --threads ${task.cpus ?: 4} \\
        --tag \\
        --force \\
        --filter \${FILTER_STRING} \\
        --subset 0 \\
        ${read}
    
    # Rename output files to match expected naming convention
    # Note: Adjust these mv commands based on actual fastq_screen output naming
    mv *.tagged.fq.gz ${sample_id}_fp1-clmp-fp2-fqscrn_r${read_num}.fq.gz
    mv *_screen.txt ${sample_id}_R${read_num}_screen.txt
    
    echo "FastQ Screen step completed for ${sample_id} R${read_num}"
    """
}