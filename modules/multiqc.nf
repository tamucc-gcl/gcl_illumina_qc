/*
 * Consolidated MultiQC module
 * Assumes upstream channel delivers a tuple:
 *   [ step , outdir , reports ]  where `reports` is a *list*
 */
process multiqc {
    label 'multiqc'
    tag   "$step"

    input:
        tuple val(step), val(outdir), val(reports)

    output:
        path("${outdir}/multiqc_${step}_report.html")

    script:
    """
    mkdir -p ${outdir}
    multiqc ${reports.join(' ')} \
            -o ${outdir} \
            -n multiqc_${step}_report.html \
            --force
    """
}
