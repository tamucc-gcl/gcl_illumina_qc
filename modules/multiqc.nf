/*
 * Consolidated MultiQC module
 * Assumes upstream channel delivers a tuple:
 *   [ step , outdir , reports ]  where `reports` is a *list*
 */
process multiqc {
  label 'multiqc'

  input:
  path files
  val step
  val outdir

  output:
    path "multiqc/${task}_report.html"

  script:
    """
    multiqc \
      ${inputs} \
      -o results/multiqc \
      -n ${task}_report.html
    """
}