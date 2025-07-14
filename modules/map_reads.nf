process map_reads {
  label 'map_reads'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)
  val genome
  val outdir

  output:
  path("${outdir}/${sample_id}.bam")

  script:
  """
  mkdir -p ${outdir}
  bwa mem ${genome} ${reads[0]} ${reads[1]} | samtools view -Sb - > ${outdir}/${sample_id}.bam
  """
}
