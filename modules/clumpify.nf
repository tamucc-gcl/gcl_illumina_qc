process clumpify {
  label 'clumpify'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("data/fq_fp1_clmp/*.fq.gz")

  script:
  """
  mkdir -p data/fq_fp1_clmp
  clumpify.sh in=${reads[0]} in2=${reads[1]} out=data/fq_fp1_clmp/${sample_id}_clumped_1.fq.gz out2=data/fq_fp1_clmp/${sample_id}_clumped_2.fq.gz
  """
}