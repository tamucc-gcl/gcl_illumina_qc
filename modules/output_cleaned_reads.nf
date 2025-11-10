// modules/output_cleaned_reads.nf
process output_cleaned_reads {
    label 'basic'
    tag "$sample_id"
    
    //publishDir "${params.outdir}/cleaned_reads", mode: 'copy'
    
    input:
        tuple val(sample_id), path(read1), path(read2)
    
    output:
        tuple val(sample_id), path("${sample_id}_cleaned.r1.fq.gz"), path("${sample_id}_cleaned.r2.fq.gz")
        path("${sample_id}_cleaned.stats")
    
    script:
    """
    # Copy with standardized naming
    cp ${read1} ${sample_id}_cleaned.r1.fq.gz
    cp ${read2} ${sample_id}_cleaned.r2.fq.gz
    
    # Generate basic stats
    echo "Sample: ${sample_id}" > ${sample_id}_cleaned.stats
    echo "R1 reads: \$(zcat ${sample_id}_cleaned.r1.fq.gz | echo \$((\$(wc -l)/4)))" >> ${sample_id}_cleaned.stats
    echo "R2 reads: \$(zcat ${sample_id}_cleaned.r2.fq.gz | echo \$((\$(wc -l)/4)))" >> ${sample_id}_cleaned.stats
    echo "R1 file size: \$(du -h ${sample_id}_cleaned.r1.fq.gz | cut -f1)" >> ${sample_id}_cleaned.stats
    echo "R2 file size: \$(du -h ${sample_id}_cleaned.r2.fq.gz | cut -f1)" >> ${sample_id}_cleaned.stats
    """
}