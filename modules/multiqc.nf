process multiqc {

    label 'multiqc'
    tag "${step}"

    input:
        path inputs
        val step
        val outdir

    output:
        path("${outdir}/multiqc_${step}_report.html")

    script:
    """
    mkdir -p ${outdir}/multiqc_${step}
    cp -r ${inputs} ${outdir}/multiqc_${step}/
    multiqc ${outdir}/multiqc_${step} \
           -o ${outdir} \
           -n multiqc_${step}_report.html \
           --force
    """
}
