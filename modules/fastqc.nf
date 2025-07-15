/*
 * FastQC on raw reads
 */
process fastqc_raw {
    label 'fastqc'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        tuple val(sample_id), 
              path("${sample_id}_R1_fastqc.html"), 
              path("${sample_id}_R1_fastqc.zip"), 
              path("${sample_id}_R2_fastqc.html"), 
              path("${sample_id}_R2_fastqc.zip")

    script:
    """
    fastqc -o . --threads ${task.cpus ?: 1} ${read1} ${read2}

    # Rename outputs to predictable sample‑centric names
    mv \$(basename ${read1} .fq.gz)_fastqc.html ${sample_id}_R1_fastqc.html
    mv \$(basename ${read1} .fq.gz)_fastqc.zip  ${sample_id}_R1_fastqc.zip
    mv \$(basename ${read2} .fq.gz)_fastqc.html ${sample_id}_R2_fastqc.html
    mv \$(basename ${read2} .fq.gz)_fastqc.zip  ${sample_id}_R2_fastqc.zip
    """
}

/*
 * Generic FastQC process for use after each QC step
 */
process fastqc_generic {
    label 'fastqc'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        tuple val(sample_id), 
              path("${sample_id}_R1_fastqc.html"), 
              path("${sample_id}_R1_fastqc.zip"), 
              path("${sample_id}_R2_fastqc.html"), 
              path("${sample_id}_R2_fastqc.zip")

    script:
    """
    fastqc -o . --threads ${task.cpus ?: 1} ${read1} ${read2}

    # Rename outputs to predictable sample‑centric names
    mv \$(basename ${read1} .fq.gz)_fastqc.html ${sample_id}_R1_fastqc.html
    mv \$(basename ${read1} .fq.gz)_fastqc.zip  ${sample_id}_R1_fastqc.zip
    mv \$(basename ${read2} .fq.gz)_fastqc.html ${sample_id}_R2_fastqc.html
    mv \$(basename ${read2} .fq.gz)_fastqc.zip  ${sample_id}_R2_fastqc.zip
    """
}