/*
 * Consolidated MultiQC module – revised
 * Expects a tuple: [ step , outdir , reports ]  (reports = list or Path)
 */
process multiqc {

    label 'multiqc'
    tag   "$step"                      // shows step name in NF UI

    input:
        tuple val(step), val(outdir), path(reports, collect: true)

    output:
        path "${outdir}/multiqc_${step}_report.html"

    script:
    """
    multiqc \
        ${reports.join(' ')} \
        -o ${outdir} \
        -n multiqc_${step}_report.html
    """
}
