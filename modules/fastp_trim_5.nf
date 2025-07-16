process fastp_trim_5 {
    label 'fastp5'
    tag   "${sample_id}"

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        tuple val(sample_id),
              path("${sample_id}_trim5_1.fq.gz"),
              path("${sample_id}_trim5_2.fq.gz"),
              path("${sample_id}_trim5_fastp.json"),
              path("${sample_id}_trim5_fastp.html")

    script:
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${sample_id}_trim5_1.fq.gz \
        --out2 ${sample_id}_trim5_2.fq.gz \
        --json ${sample_id}_trim5_fastp.json \
        --html ${sample_id}_trim5_fastp.html \
        --detect_adapter_for_pe \
        --trim_front1 0 \
        --trim_front2 0 \
        --length_required 140 \
        --cut_front \
        --cut_front_window_size 1 \
        --cut_front_mean_quality 20 \
        --cut_right \
        --cut_right_window_size 10 \
        --cut_right_mean_quality 20 \
        --disable_trim_poly_g \
        --correction \
        --disable_quality_filtering \
        --unqualified_percent_limit 40 \
        --report_title "Second Trim R1R2" \
        -w ${task.cpus ?: 4}
    """
}