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
           -out ${sample_id}.blast_raw.txt \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames" \
           -evalue 1e-10 \
           -num_threads ${task.cpus ?: 4}

    # Add full taxonomy to output
    # 1) Build ranked lineage columns from the staxids (take the first taxid if multiple)
    export TAXONKIT_DB="${taxonomy_db}"
    cut -f13 "${sample_id}.blast_raw.txt" \
    | sed 's/;.*//' \
    | taxonkit reformat -I 1 -r NA -f $'{K}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' \
    | cut -f2- \
    > "${sample_id}.taxonkit.tsv"
    #            ^^^^^^
    # cut off the leading taxid column; keep the 7 rank columns

    # 2) Paste lineage columns back onto the BLAST table
    paste "${sample_id}.blast_raw.txt" "${sample_id}.taxonkit.tsv" \
    > "${sample_id}.blast_with_taxa.tsv"

    # 3) (Optional) Add a header
    printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxids\tsscinames\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n" \
    | cat - "${sample_id}.blast_with_taxa.tsv" \
    > "${sample_id}.blast_results.txt"
    """
}