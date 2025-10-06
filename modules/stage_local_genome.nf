// modules/stage_local_genome.nf
// Process to stage a local genome file for indexing

process stage_local_genome {
    label 'stage_genome'
    tag "$genome_file"
    
    publishDir "${params.outdir}/genome", mode: 'copy', pattern: "genome.fa"

    input:
        path genome_file

    output:
        tuple path("genome.fa"), val("genome.fa")

    script:
    """
    # Copy or symlink the genome file with a standardized name
    if [[ "${genome_file}" == *.gz ]]; then
        # If the file is gzipped, decompress it
        gunzip -c ${genome_file} > genome.fa
    elif [[ "${genome_file}" == *.fa ]] || [[ "${genome_file}" == *.fasta ]] || [[ "${genome_file}" == *.fna ]]; then
        # If it's already a fasta file, just copy/link it
        cp ${genome_file} genome.fa
    else
        echo "Error: Unrecognized genome file format. Please provide .fa, .fasta, .fna, or .fa.gz file"
        exit 1
    fi
    
    # Quick validation that it looks like a FASTA file
    if ! head -n1 genome.fa | grep -q "^>"; then
        echo "Error: File does not appear to be in FASTA format (first line should start with '>')"
        exit 1
    fi
    
    echo "Staged genome file: \$(wc -l genome.fa | cut -d' ' -f1) lines"
    """
}