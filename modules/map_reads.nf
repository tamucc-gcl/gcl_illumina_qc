process map_reads {
    label 'map_reads'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(read1), path(read2)
        tuple path(genome), val(genome_path), path(index_files)

    output:
        tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    """
    bwa-mem2 mem -t ${task.cpus ?: 8} \
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
        ${genome} \
        ${read1} \
        ${read2} | \
    samtools view -@ ${task.cpus ?: 8} -Sb - | \
    samtools sort -@ ${task.cpus ?: 8} -o ${sample_id}.sorted.bam -

    # Mark duplicates
    samtools markdup -@ ${task.cpus ?: 8} ${sample_id}.sorted.bam ${sample_id}.bam
    rm ${sample_id}.sorted.bam

    # Index the BAM file
    samtools index ${sample_id}.bam
    """
}