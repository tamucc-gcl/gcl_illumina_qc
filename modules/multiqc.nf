process multiqc {
  label 'multiqc'

  input:
  path files
  val step
  val outdir

  output:
  path("${outdir}/multiqc_${step}_report.html")

  script:
  """
  mkdir -p ${outdir}/multiqc_${step}
  cp -r $files ${outdir}/multiqc_${step}/
  multiqc ${outdir}/multiqc_${step} -o ${outdir} -n multiqc_${step}_report.html --force
  """
}