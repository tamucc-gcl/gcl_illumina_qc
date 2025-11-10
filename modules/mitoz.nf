process mitoz {
    label 'mitoz'
    tag   "${sample_id}"

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        tuple val(sample_id),
              path("${sample_id}_fp1-clmp-fp2.r1.fq.gz"),
              path("${sample_id}_fp1-clmp-fp2.r2.fq.gz"),
              path("${sample_id}_trim5_fastp.json"),
              path("${sample_id}_trim5_fastp.html")

    script:
    """

    #Run MitoZ
    mitoz all \
        --outprefix ${sample_id} \
        --thread_number ${task.cpus ?: 8} \
        --clade Chordata \
        --requiring_taxa 7711 \
        --skip_filter \
        --fq1 ${read1} \
        --fq2 ${read2}

    #Process and classify into either failure or success


    """
}