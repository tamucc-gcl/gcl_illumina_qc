process fastq_screen {
  label 'fastq_screen'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("data/fq_fp1_clmp_fp2_scrn/*screen.txt")

  script:
  """
  mkdir -p data/fq_fp1_clmp_fp2_scrn
  fastq_screen --outdir data/fq_fp1_clmp_fp2_scrn --threads 4 ${reads[0]}
  fastq_screen --outdir data/fq_fp1_clmp_fp2_scrn --threads 4 ${reads[1]}
  """
}