process repair {
    label 'repair'
    tag "$sample_id"

    // publishDir "${params.outdir}/fqgz/repair", mode: params.publish_dir_mode

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        tuple val(sample_id),
              path("${sample_id}_fp1-clmp-fp2-fqscrn-rprd.r1.fq.gz"),
              path("${sample_id}_fp1-clmp-fp2-fqscrn-rprd.r2.fq.gz")

    script:
    // PIN the JVM heap explicitly. BBTools (repair.sh -> SplitPairsAndSingles, and
    // rename.sh) autodetects available RAM from the cgroup slice visible at launch,
    // NOT the SLURM allocation; when that slice is tiny it sets a tiny/garbage -Xmx
    // (e.g. -Xmx44m) and OOMs on larger samples. Derive the heap from task.memory
    // (set by the 'repair' label in nextflow.config), leaving ~2 GB headroom.
    def heap_gb = (task.memory ? Math.max(1, (task.memory.toGiga() as int) - 2) : 8)

    """
    set -euo pipefail

    repair.sh \\
        -Xmx${heap_gb}g \\
        in1=${read1} \\
        in2=${read2} \\
        out1=${sample_id}_tmp.r1.fq.gz \\
        out2=${sample_id}_tmp.r2.fq.gz \\
        outs=stdout \\
        overwrite=t \\
        repair

    rename.sh \\
        -Xmx${heap_gb}g \\
        ow=t \\
        in1="${sample_id}_tmp.r1.fq.gz" \\
        in2="${sample_id}_tmp.r2.fq.gz" \\
        out1="${sample_id}_fp1-clmp-fp2-fqscrn-rprd.r1.fq.gz" \\
        out2="${sample_id}_fp1-clmp-fp2-fqscrn-rprd.r2.fq.gz"
    """
}
