// modules/blast_mito_genes.nf
process blast_mito_genes {
    label 'blast_mito'
    tag "$sample_id"
    
    publishDir "${params.outdir}/species_id/blast_results", mode: 'copy', pattern: "*.blast_results.txt"
    
    input:
        tuple val(sample_id), path(mito_genes_fasta)
        val(blast_db)  // Path to BLAST database (or could be a string for remote DB)
    
    output:
        tuple val(sample_id), path("${sample_id}.blast_results.txt"), emit: blast_results
    
    script:
    def blast_dir = blast_db.substring(0, blast_db.lastIndexOf('/'))
    """
    export BLASTDB="${blast_dir}"    #this makes it so the taxonomy database is properly associated
    blastn -query ${mito_genes_fasta} \
           -db ${blast_db} \
           -out ${sample_id}.blast_results.txt \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames" \
           -evalue 1e-10 \
           -num_threads ${task.cpus ?: 4}
    """
}