process map_reads {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)
  path genome from params.genome

  output:
  path("${sample_id}.bam")

  script:
  """
  bwa mem $genome ${reads[0]} ${reads[1]} | samtools view -Sb - > ${sample_id}.bam
  """
}
