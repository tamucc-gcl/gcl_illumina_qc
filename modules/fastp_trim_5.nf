process fastp_trim_5 {
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
  fastp -i ${reads[0]} -I ${reads[1]} -o ${outdir}/${sample_id}_trim5_1.fq.gz -O ${outdir}/${sample_id}_trim5_2.fq.gz --trim_front1 5 --trim_front2 5 -w 4
  """
}
