process fetch_genome {
    label 'fetch_genome'
    tag "$accession"

    input:
        val accession

    output:
        tuple path("genome.fa"), val("genome.fa")

    script:
    """
    mkdir -p raw_data
    datasets download genome accession ${accession} --filename genome.zip --include genome
    unzip genome.zip -d raw_data/
    find raw_data/ -name '*.fna' -exec cat {} + > genome.fa
    """
}