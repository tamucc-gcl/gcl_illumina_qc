/*
 * FastQC on raw reads
 */
process fastqc_raw {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_R*_fastqc.*")

    script:
    """
    fastqc -o . --threads ${task.cpus ?: 1} ${reads[0]} ${reads[1]}

    # Rename outputs to predictable sample‑centric names
    mv \$(basename ${reads[0]} .fq.gz)_fastqc.html ${sample_id}_R1_fastqc.html
    mv \$(basename ${reads[0]} .fq.gz)_fastqc.zip  ${sample_id}_R1_fastqc.zip
    mv \$(basename ${reads[1]} .fq.gz)_fastqc.html ${sample_id}_R2_fastqc.html
    mv \$(basename ${reads[1]} .fq.gz)_fastqc.zip  ${sample_id}_R2_fastqc.zip
    """
}
