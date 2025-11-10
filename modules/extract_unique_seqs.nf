// modules/extract_unique_seqs.nf
process extract_unique_seqs {
    label 'uniq_seqs'
    tag "$sample_id"
    
    //publishDir "${params.outdir}/denovo_assembly/unique_seqs", mode: 'copy', pattern: "*.uniq.seqs"
    
    input:
        tuple val(sample_id), path(read1), path(read2)
    
    output:
        tuple val(sample_id), path("${sample_id}.uniq.seqs"), emit: uniq_seqs
        //path("${sample_id}.forward"), emit: forward
        //path("${sample_id}.reverse"), emit: reverse
    
    script:
    """
    # Define AWK commands as variables
    AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
    AWK2='!/>/'
    AWK3='!/NNN/'
    PERLT='while (<>) {chomp; \$z{\$_}++;} while((\$k,\$v) = each(%z)) {print "\$v\\t\$k\\n";}'
    
    # Convert fastq to fasta for forward reads
    echo "Processing ${sample_id} forward reads..."
    zcat ${read1} | mawk "\$AWK1" | mawk "\$AWK2" > ${sample_id}.forward
    
    # Convert fastq to fasta for reverse reads  
    echo "Processing ${sample_id} reverse reads..."
    zcat ${read2} | mawk "\$AWK1" | mawk "\$AWK2" > ${sample_id}.reverse
    
    # Join forward and reverse reads with N's between them, then get unique sequences with counts
    echo "Joining and finding unique sequences..."
    paste -d '-' ${sample_id}.forward ${sample_id}.reverse | \
        mawk "\$AWK3" | \
        sed 's/-/NNNNNNNNNN/' | \
        perl -e "\$PERLT" > ${sample_id}.uniq.seqs
    
    # Report statistics
    echo "Sample: ${sample_id}" > ${sample_id}.stats
    echo "Forward sequences: \$(wc -l < ${sample_id}.forward)" >> ${sample_id}.stats
    echo "Reverse sequences: \$(wc -l < ${sample_id}.reverse)" >> ${sample_id}.stats
    echo "Unique sequences: \$(wc -l < ${sample_id}.uniq.seqs)" >> ${sample_id}.stats
    
    cat ${sample_id}.stats
    """
}