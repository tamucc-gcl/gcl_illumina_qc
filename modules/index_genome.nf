process index_genome {
  label 'index_genome'
  tag "$genome_path"
  
  publishDir "${params.outdir}/genome", mode: 'copy'

  input:
  tuple path(genome), val(genome_path)

  output:
  tuple path(genome), val(genome_path)
  path("${genome}.*"), emit: index_files

  script:
  """
  bwa-mem2 index ${genome}
  """
}