process fastp_trim_5 {
  label 'fastp'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("data/fq_fp1_clmp_fp2/*.fq.gz")

  script:
  """
  mkdir -p data/fq_fp1_clmp_fp2
  fastp -i ${reads[0]} -I ${reads[1]} -o data/fq_fp1_clmp_fp2/${sample_id}_trim5_1.fq.gz -O data/fq_fp1_clmp_fp2/${sample_id}_trim5_2.fq.gz --trim_front1 5 --trim_front2 5 -w 4
  """
}