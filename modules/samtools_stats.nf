process samtools_stats {
    label 'samtools_stats'
    tag "$sample_id"

    // publishDir "${params.outdir}/qc/bamstats", mode: params.publish_dir_mode

    input:
        tuple val(sample_id), path(bam), path(bam_index)

    output:
        tuple val(sample_id), path("${sample_id}.stats"), path("${sample_id}.flagstats")
        tuple val(sample_id), path("${sample_id}_soft_clipping_stats.txt")
        tuple val(sample_id), path("${sample_id}_alignment_score_stats.txt")

    script:
    """
    # Generate samtools stats (detailed statistics)
    samtools stats -@ ${task.cpus ?: 4} ${bam} > ${sample_id}.stats
    
    # Generate samtools flagstats (alignment summary)
    samtools flagstats -@ ${task.cpus ?: 4} ${bam} > ${sample_id}.flagstats

    # Get soft clipping stats
    samtools view -@ ${task.cpus ?: 4} ${bam} |
        awk '{
        cigar = $6
        total = 0
        while (match(cigar, /[0-9]+[MIDNSHP=X]/)) {
            op  = substr(cigar, RSTART, RLENGTH)
            len = substr(op, 1, length(op)-1) + 0
            typ = substr(op, length(op), 1)
            if (typ == "S") total += len
            cigar = substr(cigar, RSTART + RLENGTH)
        }
        print total
        }' | sort -n | uniq -c > ${sample_id}_soft_clipping_stats.txt

    # Generate alignment score stats 
    samtools view -@ ${task.cpus ?: 4} ${bam} |
        awk '{
        for(i=12;i<=NF;i++){
            if($i ~ /^AS:i:/){
            split($i,a,":");
            print a[3];
            }
        }
        }' | sort -n | uniq -c > ${sample_id}_alignment_score_stats.txt
    """
}