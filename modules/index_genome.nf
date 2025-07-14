process index_genome {
  label 'map_reads'
  tag "$genome"

  input:
  path genome from genome/*.fa

  output:
  path(genome)
  path("${genome}.*")
  val(genome)

  script:
  """
  bwa index ${genome}
  """
}
