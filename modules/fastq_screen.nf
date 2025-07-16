process fastq_screen {
    label 'fastq_screen'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(read), val(read_num)

    output:
        tuple val(sample_id), 
              path("${sample_id}_screen_${read_num}.fq.gz"), 
              path("${sample_id}_R${read_num}_screen.txt"),
              val(read_num)

    script:
    """
    fastq_screen \\
        --aligner bowtie2 \\
        --conf ${params.decontam_conffile} \\
        --threads ${task.cpus ?: 4} \\
        --tag \\
        --force \\
        --filter 000000000000 \\
        --subset 0 \\
        ${read}
    
    # Rename output files to match expected naming convention
    # Note: Adjust these mv commands based on actual fastq_screen output naming
    mv *.tagged.fq.gz ${sample_id}_screen_${read_num}.fq.gz
    mv *_screen.txt ${sample_id}_R${read_num}_screen.txt
    
    echo "FastQ Screen step completed for ${sample_id} R${read_num}"
    """
}