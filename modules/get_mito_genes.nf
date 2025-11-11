process get_mito_genes {
    label 'get_mito_genes'
    tag   "${sample_id}"

    input:
        tuple val(sample_id), path(read1), path(read2)
        path mito_ref from file(params.mito_gene_refs)   // stages the fasta
        val(genetic_code)

    output:
        tuple val(sample_id),
              path("${sample_id}_mito_gene.fasta")

    script:
    """
    #Unzip & remove space from readnames
    gunzip -c ${read1} | sed 's/ /_/' > ${sample_id}.1.fq
    gunzip -c ${read2} | sed 's/ /_/' > ${sample_id}.2.fq

    #Get mito seqs
    MitoGeneExtractor \
        --exonerate_program exonerate \
        --verbosity 3 \
        --genetic_code ${genetic_code} \
        -p ${mito_reference} \
        -q ${sample_id}.1.fq \
        -q ${sample_id}.2.fq \
        -o ${sample_id}_seqMap \
        -c concensus_
    
    #clean-up
    cat concensus_*fas > ${sample_id}_combined.fasta

    #Remove empty sequences
    awk '
    /^>/ {
        if (seq) {
        print header
        print seq
        }
        header = \$0
        seq = ""
        next
    }
    { seq = seq \$0 }
    END {
        if (seq) {
        print header
        print seq
        }
    }
    ' ${sample_id}_combined.fasta > ${sample_id}_mito_gene.fasta
    """
}