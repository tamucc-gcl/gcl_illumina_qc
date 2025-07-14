process fastp_trim_3 {
  label 'fastp'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("data/fq_fp1/*.fq.gz")

  script:
  """
  mkdir -p data/fq_fp1
  fastp -i ${reads[0]} -I ${reads[1]} -o data/fq_fp1/${sample_id}_trim3_1.fq.gz -O data/fq_fp1/${sample_id}_trim3_2.fq.gz --trim_tail1 1 --trim_tail2 1 -w 4
  """
}