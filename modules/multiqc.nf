/*
 * modules/multiqc.nf
 * Consolidate FastQC/fastp/… reports into one MultiQC page.
 */
process multiqc {

    tag "$step"

    input:
        tuple val(step), val(outdir), path(reports) collect: true

    output:
        path "${outdir}/multiqc_${step}_report.html"

    script:
    """
    multiqc ${reports.join(' ')} \
            -o ${outdir} \
            -n multiqc_${step}_report.html
    """
}
