process map_reads {
    label 'map_reads'
    tag "$sample_id"

    publishDir "${params.outdir}/bam", mode: params.publish_dir_mode

    input:
        tuple val(sample_id), path(read1), path(read2)
        tuple path(genome), val(genome_path), path(index_files)
        val(sequencing_type)

    output:
        tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    // ddRAD CAVEAT: in double-digest RAD both read ends sit at fixed restriction
    // cut sites, so every read from a locus maps to IDENTICAL coordinates. PCR
    // duplicates are therefore indistinguishable from independent observations of
    // the locus, and `samtools markdup` would flag ~all-but-one read per locus as a
    // duplicate. Both genotypers we use then DROP those reads by default
    // (FreeBayes excludes duplicate-flagged reads unless --use-duplicate-reads;
    // ANGSD's default -remove_bads 1 discards flag>=256, and duplicate=1024),
    // collapsing per-locus coverage to ~1. dDocent never marks duplicates for RAD
    // for exactly this reason. So: SKIP markdup for ddRAD; run it normally for
    // randomly-sheared libraries (whole_genome), where duplicates are identifiable
    // and removal is correct.
    def is_ddrad = sequencing_type == 'ddrad'
    if (is_ddrad)
        """
        # 1) Align -> BAM
        bwa-mem2 mem -t ${task.cpus ?: 8} \
            -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
            ${genome} ${read1} ${read2} \
        | samtools view -@ ${task.cpus ?: 8} -b -o ${sample_id}.unsorted.bam -

        # 2) Coordinate-sort (no markdup for ddRAD — see note above)
        samtools sort -@ ${task.cpus ?: 8} -o ${sample_id}.bam ${sample_id}.unsorted.bam

        # 3) Index final BAM
        samtools index ${sample_id}.bam

        echo "ddRAD: duplicate marking SKIPPED (PCR duplicates are indistinguishable from real coverage at fixed cut sites)"

        rm -f ${sample_id}.unsorted.bam
        """
    else
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

        # 5) Mark duplicates (marks; use -r to remove). Correct for randomly-sheared
        #    libraries where duplicates share both ends and are identifiable.
        samtools markdup -@ ${task.cpus ?: 8} ${sample_id}.csrt.bam ${sample_id}.bam

        # 6) Index final BAM
        samtools index ${sample_id}.bam

        # Optional: clean up intermediates to save space
        rm -f ${sample_id}.unsorted.bam ${sample_id}.nsrt.bam ${sample_id}.fxmt.bam ${sample_id}.csrt.bam
        """
}
