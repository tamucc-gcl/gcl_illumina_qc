process multiqc {
    label 'multiqc'
    tag   "$step"

    input:
        path(reports)
        val(step)

    output:
        path("${params.multiqc_dir}/multiqc_${step}_report.html")

    script:
    """
    mkdir -p ${params.multiqc_dir}
    multiqc ${reports.join(' ')} \
           -o ${params.multiqc_dir} \
           -n multiqc_${step}_report.html \
           -f
    """
}