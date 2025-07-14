process map_reads {
  label 'map_reads'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)
  val genome

  output:
  path("data/fq_fp1_clmp_fp2_scrn_rpr_map/${sample_id}.bam")

  script:
  """
  mkdir -p data/fq_fp1_clmp_fp2_scrn_rpr_map
  bwa mem ${genome} ${reads[0]} ${reads[1]} | samtools view -Sb - > data/fq_fp1_clmp_fp2_scrn_rpr_map/${sample_id}.bam
  """
}