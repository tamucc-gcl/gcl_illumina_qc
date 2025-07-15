process multiqc {
    label 'multiqc'
    tag   "$step"

    input:
        path(reports)
        val(step)

    output:
        path("multiqc_${step}_report.html")

    script:
    """
    multiqc ${reports.join(' ')} \
           -o . \
           -n multiqc_${step}_report.html \
           -f
    """
}