// modules/sweep_map_score.nf
// Score one candidate reference (meta = c1,c2,sim,id) by mapping a subset of
// cleaned reads back to it and summarizing alignment quality into one row.
// Reads are a flat collected list paired internally by filename prefix.

process sweep_map_score {
    label 'sweep_map_score'
    tag "${meta.id}"

    publishDir "${params.outdir}/denovo_assembly/sweep/${meta.id}", mode: params.publish_dir_mode, pattern: "score_*.tsv"

    input:
        tuple val(meta), path(reference), path(reads)

    output:
        path("score_${meta.id}.tsv"), emit: score

    script:
    """
    set -euo pipefail

    bwa-mem2 index ${reference}
    samtools faidx ${reference}

    : > per_sample_metrics.tsv
    for r1 in *.r1.fq.gz; do
        pre=\$(basename "\$r1" .r1.fq.gz)
        r2="\${pre}.r2.fq.gz"
        [ -f "\$r2" ] || continue
        sid=\$(echo "\$pre" | sed 's/_fp1.*//')

        bwa-mem2 mem -t ${task.cpus ?: 8} ${reference} "\$r1" "\$r2" \
            | samtools sort -@ ${task.cpus ?: 8} -o "\${sid}.bam" -
        samtools index "\${sid}.bam"
        samtools stats -@ ${task.cpus ?: 8} "\${sid}.bam" > "\${sid}.stats"

        raw=\$(grep "^SN.*raw total sequences:" "\${sid}.stats" | cut -f3)
        mapped=\$(grep "^SN.*reads mapped:" "\${sid}.stats" | head -1 | cut -f3)
        pp=\$(grep "^SN.*reads properly paired:" "\${sid}.stats" | cut -f3)
        raw=\${raw:-0}; mapped=\${mapped:-0}; pp=\${pp:-0}

        # soft-clipped bases per primary mapped read (-F 0x904: no secondary/supp/unmapped)
        sc_total=\$(samtools view -@ ${task.cpus ?: 8} -F 0x904 "\${sid}.bam" | awk '{
            cigar=\$6; t=0
            while (match(cigar, /[0-9]+[MIDNSHP=X]/)) {
                op=substr(cigar,RSTART,RLENGTH); len=substr(op,1,length(op)-1)+0; typ=substr(op,length(op),1)
                if (typ=="S") t+=len
                cigar=substr(cigar,RSTART+RLENGTH)
            }
            sum+=t
        } END {print sum+0}')

        as_mean=\$(samtools view -@ ${task.cpus ?: 8} -F 0x904 "\${sid}.bam" | awk '{
            for(i=12;i<=NF;i++) if(\$i ~ /^AS:i:/){split(\$i,a,":"); s+=a[3]; n++}
        } END {if(n>0) printf "%.2f", s/n; else print 0}')

        map_rate=\$(awk -v m=\$mapped -v r=\$raw 'BEGIN{printf "%.4f", (r>0)? m*100/r : 0}')
        pp_rate=\$(awk -v p=\$pp -v r=\$raw 'BEGIN{printf "%.4f", (r>0)? p*100/r : 0}')
        sc_per_read=\$(awk -v s=\$sc_total -v m=\$mapped 'BEGIN{printf "%.4f", (m>0)? s/m : 0}')

        printf "%s\\t%s\\t%s\\t%s\\n" "\$map_rate" "\$pp_rate" "\$sc_per_read" "\$as_mean" >> per_sample_metrics.tsv
    done

    # Aggregate across samples -> one row for this candidate (id c1 c2 sim n map pp sc as)
    awk -v id="${meta.id}" -v c1="${meta.c1}" -v c2="${meta.c2}" -v sim="${meta.sim}" '
        { n++; map+=\$1; pp+=\$2; sc+=\$3; as+=\$4 }
        END {
            if (n>0) printf "%s\\t%s\\t%s\\t%s\\t%d\\t%.4f\\t%.4f\\t%.4f\\t%.4f\\n", id,c1,c2,sim,n,map/n,pp/n,sc/n,as/n
            else     printf "%s\\t%s\\t%s\\t%s\\t0\\t0\\t0\\t0\\t0\\n", id,c1,c2,sim
        }' per_sample_metrics.tsv > score_${meta.id}.tsv

    echo "Candidate ${meta.id}:"; cat score_${meta.id}.tsv
    """
}
