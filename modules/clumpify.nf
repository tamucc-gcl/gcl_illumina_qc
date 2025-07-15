process clumpify {
  label 'clumpify'
  tag "$sample_id"

  input:
    tuple val(sample_id), path(read1), path(read2)

  output:
    tuple val(sample_id),
          path("${sample_id}_clumped_1.fq.gz"),
          path("${sample_id}_clumped_2.fq.gz"),
          path("${sample_id}_clumpify.txt") emit: stats

  script:
  """
  clumpify.sh in=${read1} in2=${read2} \
             out=${sample_id}_clumped_1.fq.gz out2=${sample_id}_clumped_2.fq.gz \
             dedupe=t subs=0 ecc=t dupesubs=0 \
             stats=${sample_id}_clumpify.txt
  """
}
