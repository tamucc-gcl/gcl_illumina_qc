process map_reads {
    label 'map_reads'
    tag "$sample_id"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
        tuple val(sample_id), path(read1), path(read2)
        tuple path(genome), val(genome_path), path(index_files)

    output:
        tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
        path("${sample_id}_mapping_stats.txt")

    script:
    """
    # Debug: Check input files
    echo "=== INPUT FILE CHECK ===" > ${sample_id}_mapping_stats.txt
    echo "Sample ID: ${sample_id}" >> ${sample_id}_mapping_stats.txt
    echo "Read1 file: ${read1}" >> ${sample_id}_mapping_stats.txt
    echo "Read2 file: ${read2}" >> ${sample_id}_mapping_stats.txt
    echo "Genome file: ${genome}" >> ${sample_id}_mapping_stats.txt
    echo "" >> ${sample_id}_mapping_stats.txt
    
    # Check file sizes
    echo "=== FILE SIZES ===" >> ${sample_id}_mapping_stats.txt
    ls -lh ${read1} ${read2} ${genome} >> ${sample_id}_mapping_stats.txt
    echo "" >> ${sample_id}_mapping_stats.txt
    
    # Count reads
    echo "=== READ COUNTS ===" >> ${sample_id}_mapping_stats.txt
    echo "Read1 count: \$(( \$(wc -l < ${read1}) / 4 ))" >> ${sample_id}_mapping_stats.txt
    echo "Read2 count: \$(( \$(wc -l < ${read2}) / 4 ))" >> ${sample_id}_mapping_stats.txt
    echo "" >> ${sample_id}_mapping_stats.txt
    
    # Check genome index files
    echo "=== GENOME INDEX FILES ===" >> ${sample_id}_mapping_stats.txt
    ls -la ${genome}* >> ${sample_id}_mapping_stats.txt
    echo "" >> ${sample_id}_mapping_stats.txt
    
    # Run BWA-MEM2 with verbose output
    echo "=== STARTING BWA-MEM2 MAPPING ===" >> ${sample_id}_mapping_stats.txt
    date >> ${sample_id}_mapping_stats.txt
    
    bwa-mem2 mem -t ${task.cpus ?: 8} -v 1 \
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
        ${genome} \
        ${read1} \
        ${read2} 2>> ${sample_id}_mapping_stats.txt | \
    samtools view -@ ${task.cpus ?: 8} -Sb - | \
    samtools sort -@ ${task.cpus ?: 8} -o ${sample_id}.bam -
    
    # Check if BAM was created successfully
    if [ -f "${sample_id}.bam" ]; then
        echo "=== BAM FILE CREATED SUCCESSFULLY ===" >> ${sample_id}_mapping_stats.txt
        
        # Index the BAM file
        samtools index ${sample_id}.bam
        
        # Generate mapping statistics
        echo "=== MAPPING STATISTICS ===" >> ${sample_id}_mapping_stats.txt
        samtools flagstat ${sample_id}.bam >> ${sample_id}_mapping_stats.txt
        echo "" >> ${sample_id}_mapping_stats.txt
        
        echo "=== CHROMOSOME MAPPING STATS ===" >> ${sample_id}_mapping_stats.txt
        samtools idxstats ${sample_id}.bam >> ${sample_id}_mapping_stats.txt
        echo "" >> ${sample_id}_mapping_stats.txt
        
        # Check for properly paired reads
        echo "=== PAIRING STATISTICS ===" >> ${sample_id}_mapping_stats.txt
        samtools view -f 2 ${sample_id}.bam | wc -l | awk '{print "Properly paired reads: " \$1}' >> ${sample_id}_mapping_stats.txt
        samtools view -F 4 ${sample_id}.bam | wc -l | awk '{print "Mapped reads: " \$1}' >> ${sample_id}_mapping_stats.txt
        
    else
        echo "=== ERROR: BAM FILE NOT CREATED ===" >> ${sample_id}_mapping_stats.txt
        echo "BWA-MEM2 mapping failed" >> ${sample_id}_mapping_stats.txt
    fi
    
    date >> ${sample_id}_mapping_stats.txt
    echo "=== MAPPING PROCESS COMPLETED ===" >> ${sample_id}_mapping_stats.txt
    """
}