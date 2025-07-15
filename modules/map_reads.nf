process map_reads {

    label 'map_reads'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(read1), path(read2)
        path genome

    output:
        path("${sample_id}.bam")

    script:
    """
    bwa mem2 ${genome} ${read1} ${read2} | samtools view -Sb - > ${sample_id}.bam
    """
}
