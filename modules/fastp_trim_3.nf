process fastp_trim_3 {
    label 'fastp3'
    tag   "${sample_id}"

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        tuple val(sample_id),
              path("${sample_id}_trim3_1.fq.gz"),
              path("${sample_id}_trim3_2.fq.gz"),
              path("${sample_id}_trim3_fastp.json"),
              path("${sample_id}_trim3_fastp.html")

    script:
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${sample_id}_trim3_1.fq.gz \
        --out2 ${sample_id}_trim3_2.fq.gz \
        --json ${sample_id}_trim3_fastp.json \
        --html ${sample_id}_trim3_fastp.html \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --length_required 33 \
        --low_complexity_filter \
        --complexity_threshold 30 \
        --detect_adapter_for_pe \
        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --cut_tail \
        --cut_tail_window_size 1 \
        --cut_tail_mean_quality 20 \
        --trim_poly_g \
        --poly_g_min_len 10 \
        --trim_poly_x \
        --report_title "First Trim 4 De Novo" \
        -w ${task.cpus ?: 4}
    """
}