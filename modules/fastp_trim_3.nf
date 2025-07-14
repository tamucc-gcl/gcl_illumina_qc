process fastp_trim_3 {
  label 'fastp'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)
  val outdir

  output:
  tuple val(sample_id), path("${outdir}/*.fq.gz")

  script:
  """
  mkdir -p ${outdir}
  fastp -i ${reads[0]} -I ${reads[1]} -o ${outdir}/${sample_id}_trim3_1.fq.gz -O ${outdir}/${sample_id}_trim3_2.fq.gz --trim_tail1 1 --trim_tail2 1 -w 4
  """
}
