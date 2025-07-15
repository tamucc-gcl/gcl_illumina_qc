/*
 * Consolidated MultiQC module
 * Assumes upstream channel delivers a tuple:
 *   [ step , outdir , reports ]  where `reports` is a *list*
 */
process multiqc {
<<<<<<< HEAD
    label 'multiqc'
    tag   "$step"

    input:
        tuple val(step), val(outdir), val(reports)
=======
  label 'multiqc'

  input:
  path files
  val step
  val outdir

  output:
    path "multiqc/${task}_report.html"
>>>>>>> parent of 1fd308f (updates)

  script:
    """
<<<<<<< HEAD
    mkdir -p ${outdir}
    multiqc ${reports.join(' ')} \
            -o ${outdir} \
            -n multiqc_${step}_report.html \
            --force
=======
    multiqc \
      ${inputs} \
      -o results/multiqc \
      -n ${task}_report.html
>>>>>>> parent of 1fd308f (updates)
    """
}