process map_reads {
  label 'map_reads'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)
  val genome

  output:
  path("${sample_id}.bam")

  script:
  """
  bwa mem2 ${genome} ${reads[0]} ${reads[1]} | \
    samtools view -Sb - > ${sample_id}.bam
  """
}