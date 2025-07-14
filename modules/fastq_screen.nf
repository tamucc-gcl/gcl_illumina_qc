process fastq_screen {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("*screen.txt")

  script:
  """
  fastq_screen --outdir . --threads 4 ${reads[0]}
  fastq_screen --outdir . --threads 4 ${reads[1]}
  """
}
