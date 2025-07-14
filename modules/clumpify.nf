process clumpify {
  label 'clumpify'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)
  val outdir

  output:
  tuple val(sample_id), path("${outdir}/*.fq.gz")

  script:
  """
  mkdir -p ${outdir}
  clumpify.sh in=${reads[0]} in2=${reads[1]} out=${outdir}/${sample_id}_clumped_1.fq.gz out2=${outdir}/${sample_id}_clumped_2.fq.gz
  """
}
