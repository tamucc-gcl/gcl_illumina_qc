// modules/assemble_rainbow_candidate.nf
// Sweep variant of assemble_rainbow. Candidate identity (cutoff1, cutoff2,
// cluster_similarity) is carried in a meta map; outputs/publishDir keyed by meta.id.
//
// meta = [c1:<int>, c2:<int>, sim:<dbl>, id:"c<c1>_k<c2>_s<sim>", ...]

process assemble_rainbow_candidate {
    label 'denovo_assembly'
    tag "${meta.id}"

    publishDir "${params.outdir}/denovo_assembly/optimize/candidates/${meta.id}", mode: params.publish_dir_mode

    input:
        tuple val(meta), path(uniq_fasta), path(totaluniqseq)
        val(div_f)
        val(div_K)
        val(merge_r)
        val(final_similarity)

    output:
        tuple val(meta), path("denovo_reference.fa"), emit: reference
        tuple val(meta), path("assembly_stats.txt"),  emit: stats

    script:
    def cdhit_n = { c -> def t = c as double
                    t >= 0.95 ? 10 : t >= 0.90 ? 8 : t >= 0.88 ? 7 :
                    t >= 0.85 ? 6  : t >= 0.80 ? 5 : 4 }
    """
    echo "Sweep candidate ${meta.id} (cutoff1=${meta.c1}, cutoff2=${meta.c2}, sim=${meta.sim})"

    sed -e 's/NNNNNNNNNN/\\t/g' ${uniq_fasta} | cut -f1 > uniq.F.fasta

    cd-hit-est -i uniq.F.fasta -o xxx -c ${meta.sim} -n ${cdhit_n(meta.sim)} \
               -T ${task.cpus} -M 0 -g 1

    mawk '{if (\$1 ~ /Cl/) clus = clus + 1; else print \$3 "\\t" clus}' xxx.clstr | \
        sed 's/[>Contig_,...]//g' | sort -g -k1 > sort.contig.cluster.ids

    paste sort.contig.cluster.ids ${totaluniqseq} > contig.cluster.totaluniqseq
    sort -k2,2 -g contig.cluster.totaluniqseq | sed -e 's/NNNNNNNNNN/\\t/g' > rcluster

    rainbow div -i rcluster -o rbdiv.out -f ${div_f} -K ${div_K}
    rainbow merge -o rbasm.out -a -i rbdiv.out -r ${merge_r}

    cat rbasm.out <(echo "E") | sed 's/[0-9]*:[0-9]*://g' | mawk '
    {
        if (NR == 1) e=\$2;
        else if (\$1 ~/E/ && lenp > len1) {
            c=c+1; print ">dDocent_Contig_" e "\\n" seq2 "NNNNNNNNNN" seq1;
            seq1=0; seq2=0; lenp=0; e=\$2; fclus=0; len1=0; freqp=0; lenf=0
        }
        else if (\$1 ~/E/ && lenp <= len1) {
            c=c+1; print ">dDocent_Contig_" e "\\n" seq1;
            seq1=0; seq2=0; lenp=0; e=\$2; fclus=0; len1=0; freqp=0; lenf=0
        }
        else if (\$1 ~/C/) clus=\$2;
        else if (\$1 ~/L/) len=\$2;
        else if (\$1 ~/S/) seq=\$2;
        else if (\$1 ~/N/) freq=\$2;
        else if (\$1 ~/R/ && \$0 ~/0/ && \$0 !~/1/ && len > lenf) { seq1 = seq; fclus=clus; lenf=len }
        else if (\$1 ~/R/ && \$0 ~/0/ && \$0 ~/1/) { seq1 = seq; fclus=clus; len1=len }
        else if (\$1 ~/R/ && \$0 !~/0/ && freq > freqp && len >= lenp || \$1 ~/R/ && \$0 !~/0/ && freq == freqp && len > lenp) {
            seq2 = seq; lenp = len; freqp=freq
        }
    }' > rainbow.fasta

    cd-hit-est -i rainbow.fasta -o denovo_reference.fa -M 0 -T ${task.cpus} \
               -c ${final_similarity} -n ${cdhit_n(final_similarity)}

    cat > assembly_stats.txt <<EOF
Rainbow Assembly Statistics (sweep candidate ${meta.id})
========================================================
Assembly Parameters:
  cutoff1 (min reads/individual): ${meta.c1}
  cutoff2 (min individuals): ${meta.c2}
  Initial clustering: ${meta.sim}
  Rainbow div -f: ${div_f}
  Rainbow div -K: ${div_K}
  Rainbow merge -r: ${merge_r}
  Final clustering: ${final_similarity}

Input sequences: \$(grep -c '^>' ${uniq_fasta})
Final reference contigs: \$(grep -c '^>' denovo_reference.fa)

Final Assembly Metrics:
EOF

    # Calculate contig length distribution
    awk '/^>/ {if (seq) print length(seq); seq=""; next} {seq=seq\$0} END {if (seq) print length(seq)}' \
        denovo_reference.fa | sort -rn > contig_lengths.txt

    # N50 and other metrics
    total_length=\$(awk '{sum+=\$1} END {print sum}' contig_lengths.txt)
    num_contigs=\$(wc -l < contig_lengths.txt)
    half_length=\$(echo "\$total_length / 2" | bc)
    n50=\$(awk -v half=\$half_length '{sum+=\$1; if(sum>=half && !found) {print \$1; found=1}}' contig_lengths.txt)

    echo "  Total bases: \$total_length" >> assembly_stats.txt
    echo "  Total contigs: \$num_contigs" >> assembly_stats.txt
    echo "  N50: \$n50" >> assembly_stats.txt
    echo "  Min contig: \$(tail -1 contig_lengths.txt) bp" >> assembly_stats.txt
    echo "  Max contig: \$(head -1 contig_lengths.txt) bp" >> assembly_stats.txt
    echo "  Mean contig: \$(awk -v n=\$num_contigs -v t=\$total_length 'BEGIN{printf "%.0f", t/n}') bp" >> assembly_stats.txt

    # Paired vs forward-only loci (PE spacer NNNNNNNNNN marks paired assembly)
    paired_loci=\$(grep -v '^>' denovo_reference.fa | grep -c 'NNNNNNNNNN' || true)
    echo "  Paired-end loci (with N-spacer): \$paired_loci" >> assembly_stats.txt
    echo "  Forward-only loci: \$(awk -v c=\$num_contigs -v p=\$paired_loci 'BEGIN{print c-p}')" >> assembly_stats.txt

    # Contig size distribution
    echo "" >> assembly_stats.txt
    echo "Contig Size Distribution:" >> assembly_stats.txt
    echo "  >10kb: \$(awk '\$1>10000' contig_lengths.txt | wc -l)" >> assembly_stats.txt
    echo "  5-10kb: \$(awk '\$1>=5000 && \$1<=10000' contig_lengths.txt | wc -l)" >> assembly_stats.txt
    echo "  1-5kb: \$(awk '\$1>=1000 && \$1<5000' contig_lengths.txt | wc -l)" >> assembly_stats.txt
    echo "  500bp-1kb: \$(awk '\$1>=500 && \$1<1000' contig_lengths.txt | wc -l)" >> assembly_stats.txt
    echo "  <500bp: \$(awk '\$1<500' contig_lengths.txt | wc -l)" >> assembly_stats.txt

    cat assembly_stats.txt
    rm -f contig_lengths.txt xxx xxx.clstr
    """
}
