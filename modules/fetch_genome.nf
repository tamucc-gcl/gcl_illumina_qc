process fetch_genome {
  label 'fetch_genome'
  tag "$accession"

  input:
  val accession

  output:
  path "data/genome/genome.fa", emit: genome_fasta

script:
  """
  mkdir -p data/genome
  datasets download genome accession ${accession} --filename temp.zip --include genome
  unzip temp.zip -d temp_dir
  mv temp_dir/ncbi_dataset/data/*/*.fna data/genome/genome.fa
  """
}
