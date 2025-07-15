process multiqc {

    label 'multiqc'
    tag   "$step"

    /*
     * We receive a single channel that already holds:
     *   [ step , outdir , list_of_reports ]
     * `collect: true` turns the *reports* element into a list.
     */
    input:
        tuple val(step), val(outdir), path(reports, collect: true)

    output:
        path "${outdir}/multiqc_${step}_report.html"

    script:
    """
    multiqc ${reports.join(' ')} \
           -o ${outdir} \
           -n multiqc_${step}_report.html
    """
}
