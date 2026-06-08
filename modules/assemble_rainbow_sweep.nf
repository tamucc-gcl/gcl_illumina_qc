// modules/assemble_rainbow_sweep.nf
// Sweep variant of assemble_rainbow. Candidate identity (cutoff1, cutoff2,
// cluster_similarity) is carried in a meta map; outputs/publishDir keyed by meta.id.
//
// meta = [c1:<int>, c2:<int>, sim:<dbl>, id:"c<c1>_k<c2>_s<sim>", ...]

process assemble_rainbow_sweep {
    label 'denovo_assembly'
    tag "${meta.id}"

    publishDir "${params.outdir}/denovo_assembly/sweep/${meta.id}", mode: params.publish_dir_mode

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
EOF
    cat assembly_stats.txt
    rm -f xxx xxx.clstr
    """
}
