// modules/assemble_rainbow_sweep.nf
// Sweep variant of assemble_rainbow: identical Rainbow logic, but the initial
// CD-HIT clustering similarity is carried as a per-candidate tag (`sim`) so the
// process can fan out over a grid of cluster_similarity values. Outputs and the
// publish dir are keyed by `sim` to avoid collisions between candidates.
//
// NOTE: the assembly body is duplicated from assemble_rainbow.nf. If this proves
// stable we can unify both into one tagged module + shared bin/ script.

process assemble_rainbow_sweep {
    label 'denovo_assembly'
    tag "sweep_sim_${sim}"

    publishDir "${params.outdir}/denovo_assembly/sweep/sim_${sim}", mode: params.publish_dir_mode

    input:
        tuple val(sim), path(uniq_fasta), path(totaluniqseq)
        val(div_f)
        val(div_K)
        val(merge_r)
        val(final_similarity)

    output:
        tuple val(sim), path("denovo_reference.fa"), emit: reference
        tuple val(sim), path("assembly_stats.txt"),  emit: stats

    script:
    // CD-HIT-EST word size must track the similarity threshold
    def cdhit_n = { c -> def t = c as double
                    t >= 0.95 ? 10 : t >= 0.90 ? 8 : t >= 0.88 ? 7 :
                    t >= 0.85 ? 6  : t >= 0.80 ? 5 : 4 }
    """
    echo "Sweep candidate cluster_similarity = ${sim}"

    # Extract forward sequences only for initial clustering
    sed -e 's/NNNNNNNNNN/\\t/g' ${uniq_fasta} | cut -f1 > uniq.F.fasta

    # Initial clustering with CD-HIT (word size matched to similarity)
    cd-hit-est -i uniq.F.fasta \
               -o xxx \
               -c ${sim} \
               -n ${cdhit_n(sim)} \
               -T ${task.cpus} \
               -M 0 \
               -g 1

    # Convert from CD-HIT to Rainbow input
    mawk '{if (\$1 ~ /Cl/) clus = clus + 1; else print \$3 "\\t" clus}' xxx.clstr | \
        sed 's/[>Contig_,...]//g' | \
        sort -g -k1 > sort.contig.cluster.ids

    paste sort.contig.cluster.ids ${totaluniqseq} > contig.cluster.totaluniqseq

    sort -k2,2 -g contig.cluster.totaluniqseq | \
        sed -e 's/NNNNNNNNNN/\\t/g' > rcluster

    # Rainbow div / merge
    rainbow div -i rcluster -o rbdiv.out -f ${div_f} -K ${div_K}
    rainbow merge -o rbasm.out -a -i rbdiv.out -r ${merge_r}

    # Extract optimal sequences for de novo genome
    cat rbasm.out <(echo "E") | sed 's/[0-9]*:[0-9]*://g' | mawk '
    {
        if (NR == 1) e=\$2;
        else if (\$1 ~/E/ && lenp > len1) {
            c=c+1;
            print ">dDocent_Contig_" e "\\n" seq2 "NNNNNNNNNN" seq1;
            seq1=0; seq2=0; lenp=0; e=\$2; fclus=0; len1=0; freqp=0; lenf=0
        }
        else if (\$1 ~/E/ && lenp <= len1) {
            c=c+1;
            print ">dDocent_Contig_" e "\\n" seq1;
            seq1=0; seq2=0; lenp=0; e=\$2; fclus=0; len1=0; freqp=0; lenf=0
        }
        else if (\$1 ~/C/) clus=\$2;
        else if (\$1 ~/L/) len=\$2;
        else if (\$1 ~/S/) seq=\$2;
        else if (\$1 ~/N/) freq=\$2;
        else if (\$1 ~/R/ && \$0 ~/0/ && \$0 !~/1/ && len > lenf) {
            seq1 = seq; fclus=clus; lenf=len
        }
        else if (\$1 ~/R/ && \$0 ~/0/ && \$0 ~/1/) {
            seq1 = seq; fclus=clus; len1=len
        }
        else if (\$1 ~/R/ && \$0 !~/0/ && freq > freqp && len >= lenp || \$1 ~/R/ && \$0 !~/0/ && freq == freqp && len > lenp) {
            seq2 = seq; lenp = len; freqp=freq
        }
    }' > rainbow.fasta

    # Final clustering to remove redundancy (word size matched to similarity)
    cd-hit-est -i rainbow.fasta \
               -o denovo_reference.fa \
               -M 0 \
               -T ${task.cpus} \
               -c ${final_similarity} \
               -n ${cdhit_n(final_similarity)}

    # Assembly statistics
    cat > assembly_stats.txt <<-EOF
	Rainbow Assembly Statistics (sweep candidate)
	=============================================
	Assembly Parameters:
	  Initial clustering: ${sim}
	  Rainbow div -f: ${div_f}
	  Rainbow div -K: ${div_K}
	  Rainbow merge -r: ${merge_r}
	  Final clustering: ${final_similarity}
	
	Input sequences: \$(grep -c '^>' ${uniq_fasta})
	Final reference contigs: \$(grep -c '^>' denovo_reference.fa)
	
	Final Assembly Metrics:
	EOF

    awk '/^>/ {if (seq) print length(seq); seq=""; next} {seq=seq\$0} END {if (seq) print length(seq)}' \
        denovo_reference.fa | sort -rn > contig_lengths.txt

    total_length=\$(awk '{sum+=\$1} END {print sum}' contig_lengths.txt)
    num_contigs=\$(wc -l < contig_lengths.txt)
    half_length=\$(echo "\$total_length / 2" | bc)
    n50=\$(awk -v half=\$half_length '{sum+=\$1; if(sum>=half && !found) {print \$1; found=1}}' contig_lengths.txt)

    echo "  Total bases: \$total_length"  >> assembly_stats.txt
    echo "  Total contigs: \$num_contigs" >> assembly_stats.txt
    echo "  N50: \$n50"                   >> assembly_stats.txt
    echo "  Min contig: \$(tail -1 contig_lengths.txt) bp" >> assembly_stats.txt
    echo "  Max contig: \$(head -1 contig_lengths.txt) bp" >> assembly_stats.txt
    echo "  Mean contig: \$(awk -v n=\$num_contigs -v t=\$total_length 'BEGIN{printf "%.0f", t/n}') bp" >> assembly_stats.txt

    cat assembly_stats.txt
    rm -f contig_lengths.txt xxx xxx.clstr
    """
}
