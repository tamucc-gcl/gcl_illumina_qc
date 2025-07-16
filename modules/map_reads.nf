process map_reads {
    label 'map_reads'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(read1), path(read2)
        tuple path(genome), val(genome_path), path(index_files)

    output:
        tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
        path("${sample_id}_mapping_debug.txt")

    script:
    """
    echo "=== MAPPING DEBUG FOR ${sample_id} ===" > ${sample_id}_mapping_debug.txt
    echo "Date: \$(date)" >> ${sample_id}_mapping_debug.txt
    echo "" >> ${sample_id}_mapping_debug.txt
    
    echo "=== INPUT FILES ===" >> ${sample_id}_mapping_debug.txt
    echo "Read1: ${read1}" >> ${sample_id}_mapping_debug.txt
    echo "Read2: ${read2}" >> ${sample_id}_mapping_debug.txt
    echo "Genome: ${genome}" >> ${sample_id}_mapping_debug.txt
    echo "" >> ${sample_id}_mapping_debug.txt
    
    echo "=== FILE SIZES ===" >> ${sample_id}_mapping_debug.txt
    ls -lh ${read1} ${read2} ${genome} >> ${sample_id}_mapping_debug.txt
    echo "" >> ${sample_id}_mapping_debug.txt
    
    echo "=== READ COUNTS IN INPUT FILES ===" >> ${sample_id}_mapping_debug.txt
    read1_count=\$(zcat ${read1} | wc -l)
    read2_count=\$(zcat ${read2} | wc -l)
    read1_sequences=\$((\$read1_count / 4))
    read2_sequences=\$((\$read2_count / 4))
    
    echo "Read1 lines: \$read1_count (sequences: \$read1_sequences)" >> ${sample_id}_mapping_debug.txt
    echo "Read2 lines: \$read2_count (sequences: \$read2_sequences)" >> ${sample_id}_mapping_debug.txt
    echo "" >> ${sample_id}_mapping_debug.txt
    
    echo "=== FIRST FEW READS FROM READ1 ===" >> ${sample_id}_mapping_debug.txt
    zcat ${read1} | head -8 >> ${sample_id}_mapping_debug.txt
    echo "" >> ${sample_id}_mapping_debug.txt
    
    echo "=== FIRST FEW READS FROM READ2 ===" >> ${sample_id}_mapping_debug.txt
    zcat ${read2} | head -8 >> ${sample_id}_mapping_debug.txt
    echo "" >> ${sample_id}_mapping_debug.txt
    
    echo "=== GENOME INFO ===" >> ${sample_id}_mapping_debug.txt
    echo "Genome file size: \$(ls -lh ${genome} | awk '{print \$5}')" >> ${sample_id}_mapping_debug.txt
    echo "Number of sequences in genome: \$(grep -c '^>' ${genome})" >> ${sample_id}_mapping_debug.txt
    echo "First few sequence headers:" >> ${sample_id}_mapping_debug.txt
    head -10 ${genome} | grep '^>' >> ${sample_id}_mapping_debug.txt
    echo "" >> ${sample_id}_mapping_debug.txt
    
    echo "=== INDEX FILES ===" >> ${sample_id}_mapping_debug.txt
    ls -lh genome.fa.* 2>/dev/null || echo "No index files found" >> ${sample_id}_mapping_debug.txt
    echo "" >> ${sample_id}_mapping_debug.txt
    
    # Only proceed with mapping if we have reads
    if [ \$read1_sequences -gt 0 ] && [ \$read2_sequences -gt 0 ]; then
        echo "=== PROCEEDING WITH MAPPING ===" >> ${sample_id}_mapping_debug.txt
        echo "Starting bwa-mem2 mapping..." >> ${sample_id}_mapping_debug.txt
        
        bwa-mem2 mem -t ${task.cpus ?: 8} \
            -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
            ${genome} \
            ${read1} \
            ${read2} | \
        samtools view -@ ${task.cpus ?: 8} -Sb - | \
        samtools sort -@ ${task.cpus ?: 8} -o ${sample_id}.bam -
        
        # Index the BAM file
        samtools index ${sample_id}.bam
        
        echo "Mapping completed." >> ${sample_id}_mapping_debug.txt
        echo "BAM file size: \$(ls -lh ${sample_id}.bam | awk '{print \$5}')" >> ${sample_id}_mapping_debug.txt
        
        # Count mapped reads
        mapped_reads=\$(samtools view -c -F 4 ${sample_id}.bam)
        total_reads=\$(samtools view -c ${sample_id}.bam)
        echo "Total reads in BAM: \$total_reads" >> ${sample_id}_mapping_debug.txt
        echo "Mapped reads: \$mapped_reads" >> ${sample_id}_mapping_debug.txt
        echo "Unmapped reads: \$((\$total_reads - \$mapped_reads))" >> ${sample_id}_mapping_debug.txt
        
    else
        echo "=== NO READS FOUND - CREATING EMPTY BAM ===" >> ${sample_id}_mapping_debug.txt
        echo "Cannot proceed with mapping - no reads in input files!" >> ${sample_id}_mapping_debug.txt
        
        # Create empty BAM file with header only
        echo "@HD\tVN:1.6\tSO:coordinate" | samtools view -bS - > ${sample_id}.bam
        samtools index ${sample_id}.bam
    fi
    
    echo "=== DEBUG COMPLETE ===" >> ${sample_id}_mapping_debug.txt
    """
}