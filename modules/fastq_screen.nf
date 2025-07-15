process fastq_screen {
    label 'fastq_screen'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        tuple val(sample_id), 
              path("${sample_id}_screen_1.fq.gz"), 
              path("${sample_id}_screen_2.fq.gz"),
              path("${sample_id}_R1_screen.txt"),
              path("${sample_id}_R2_screen.txt")

    script:
    """
    # Run fastq_screen on both reads
    fastq_screen --threads ${task.cpus ?: 4} --outdir . ${read1}
    fastq_screen --threads ${task.cpus ?: 4} --outdir . ${read2}
    
    # Rename output files to predictable names
    mv \$(basename ${read1} .fq.gz)_screen.txt ${sample_id}_R1_screen.txt
    mv \$(basename ${read2} .fq.gz)_screen.txt ${sample_id}_R2_screen.txt
    
    # Create symbolic links for the reads (fastq_screen doesn't modify the reads)
    ln -s ${read1} ${sample_id}_screen_1.fq.gz
    ln -s ${read2} ${sample_id}_screen_2.fq.gz
    """
}