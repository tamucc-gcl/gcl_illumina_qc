process fetch_genome {
  label 'fetch_genome'
  tag "$accession"

  input:
  val accession

  output:
  path("genome/raw")   // All raw files here
  path("genome/genome.fa")

  script:
  """
  mkdir -p genome/raw
  datasets download genome accession ${accession} --filename genome.zip --include genome
  unzip genome.zip -d genome/raw/
  find genome/raw/ -name '*fna' -exec cat {} + > genome/genome.fa
  """
}
