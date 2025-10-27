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
    # 1) Align -> BAM
    bwa-mem2 mem -t ${task.cpus ?: 8} \
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
        ${genome} ${read1} ${read2} \
    | samtools view -@ ${task.cpus ?: 8} -b -o ${sample_id}.unsorted.bam -

    # 2) Name-sort (required before fixmate)
    samtools sort -@ ${task.cpus ?: 8} -n -o ${sample_id}.nsrt.bam ${sample_id}.unsorted.bam

    # 3) Fixmate (adds required ms/MC tags for paired-end)
    samtools fixmate -@ ${task.cpus ?: 8} -m ${sample_id}.nsrt.bam ${sample_id}.fxmt.bam

    # 4) Coordinate-sort for markdup
    samtools sort -@ ${task.cpus ?: 8} -o ${sample_id}.csrt.bam ${sample_id}.fxmt.bam

    # 5) Mark duplicates (marks; use -r to remove)
    samtools markdup -@ ${task.cpus ?: 8} ${sample_id}.csrt.bam ${sample_id}.bam

    # 6) Index final BAM
    samtools index ${sample_id}.bam

    # Optional: clean up intermediates to save space
    rm -f ${sample_id}.unsorted.bam ${sample_id}.nsrt.bam ${sample_id}.fxmt.bam ${sample_id}.csrt.bam
    """
}