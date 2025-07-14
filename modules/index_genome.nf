process index_genome {
  label 'index_genome'
  tag "$genome"

  input:
  path genome

  output:
  path("${genome}.*")
  val(genome)

  script:
  """
  bwa index ${genome}
  """
}
