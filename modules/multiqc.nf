/*
 * Consolidated MultiQC module
 * Groups all result files by processing step and produces a single report per step.
 */
process multiqc {
    label 'multiqc'
    tag   "$step"

    input:
        val  step                        // e.g. 'raw_fastqc'
        val  outdir                      // e.g. 'results/multiqc'
        path reports collect: true       // all HTML / JSON / txt files to include

    output:
        path("${outdir}/multiqc_${step}_report.html")

    script:
    """
    mkdir -p ${outdir}
    multiqc ${reports} \
            -o ${outdir} \
            -n multiqc_${step}_report.html \
            --force
    """
}
