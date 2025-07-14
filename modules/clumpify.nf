process clumpify {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("*.fq.gz")

  script:
  """
  clumpify.sh in=${reads[0]} in2=${reads[1]} out=${sample_id}_clumped_1.fq.gz out2=${sample_id}_clumped_2.fq.gz
  """
}
