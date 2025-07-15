process fastq_screen {
  label 'fastq_screen'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("${sample_id}_screen_*.fq.gz")

  script:
  """
  fastq_screen --threads 4 ${reads[0]}
  fastq_screen --threads 4 ${reads[1]}
  """
}