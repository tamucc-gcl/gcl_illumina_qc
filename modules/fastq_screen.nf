process fastq_screen {
  label 'fastq_screen'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)
  val outdir

  output:
  tuple val(sample_id), path("${outdir}/*screen.txt")

  script:
  """
  mkdir -p ${outdir}
  fastq_screen --outdir ${outdir} --threads 4 ${reads[0]}
  fastq_screen --outdir ${outdir} --threads 4 ${reads[1]}
  """
}
