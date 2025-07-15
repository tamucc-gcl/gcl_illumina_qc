process fastq_screen {

    label 'fastq_screen'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(read1), path(read2), path(json), path(html)

    output:
        tuple val(sample_id),
              path("${sample_id}_screen_1.fq.gz"),
              path("${sample_id}_screen_2.fq.gz"),
              path("${sample_id}_screen_1.txt"),
              path("${sample_id}_screen_2.txt")

    script:
    """
    fastq_screen --threads ${task.cpus ?: 4} ${read1}
    fastq_screen --threads ${task.cpus ?: 4} ${read2}

    cp ${read1} ${sample_id}_screen_1.fq.gz
    cp ${read2} ${sample_id}_screen_2.fq.gz
    """
}
