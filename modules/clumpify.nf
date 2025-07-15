process clumpify {
  label 'clumpify'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id),
       path("${sample_id}_clumped_*.fq.gz"),
       path("${sample_id}_clumpify.txt")
  script:
  """
  clumpify.sh in=${reads[0]} in2=${reads[1]} \
           out=${sample_id}_clumped_1.fq.gz out2=${sample_id}_clumped_2.fq.gz \
           dedupe=t subs=0 ecc=t dupesubs=0 \
           stats=${sample_id}_clumpify.txt

  """
}