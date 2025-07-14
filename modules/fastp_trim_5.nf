process fastp_trim_5 {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("*.fq.gz")

  script:
  """
  fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_trim5_1.fq.gz -O ${sample_id}_trim5_2.fq.gz --trim_front1 5 --trim_front2 5 -w 4
  """
}
