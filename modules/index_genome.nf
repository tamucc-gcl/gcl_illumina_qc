process index_genome {
  label 'index_genome'
  tag "$genome_path"

  input:
  tuple path(genome), val(genome_path)

  output:
  tuple path(genome), val(genome_path)

  script:
  """
  bwa index ${genome}
  """
}
