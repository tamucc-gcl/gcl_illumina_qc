// modules/assemble_rainbow.nf
process assemble_rainbow {
    label 'denovo_assembly'
    tag "rainbow_assembly"
    
    publishDir "${params.outdir}/denovo_assembly", mode: 'copy'
    
    input:
        path(uniq_fasta)
        path(totaluniqseq)
        val(cluster_similarity)  // CD-HIT clustering threshold (default: 0.8)
        val(div_f)               // Rainbow div -f parameter (default: 0.5)
        val(div_K)               // Rainbow div -K parameter (default: 10)
        val(merge_r)             // Rainbow merge -r parameter (default: 2)
        val(final_similarity)    // Final CD-HIT clustering (default: 0.9)
    
    output:
        path("denovo_reference.fa"), emit: reference
        path("assembly_stats.txt"), emit: stats
        //path("rainbow.fasta"), emit: rainbow_assembly
    
    script:
    """
    echo "Starting Rainbow assembly pipeline"
    echo "Parameters:"
    echo "  CD-HIT initial similarity: ${cluster_similarity}"
    echo "  Rainbow div -f: ${div_f}"
    echo "  Rainbow div -K: ${div_K}"
    echo "  Rainbow merge -r: ${merge_r}"
    echo "  Final clustering similarity: ${final_similarity}"
    
    # Extract forward sequences only for initial clustering
    sed -e 's/NNNNNNNNNN/\\t/g' ${uniq_fasta} | cut -f1 > uniq.F.fasta
    
    # Initial clustering with CD-HIT
    echo "Running initial CD-HIT clustering..."
    cd-hit-est -i uniq.F.fasta \
               -o xxx \
               -c ${cluster_similarity} \
               -T ${task.cpus} \
               -M 0 \
               -g 1
    
    # Convert from CD-HIT to Rainbow input
    echo "Converting CD-HIT output to Rainbow format..."
    mawk '{if (\$1 ~ /Cl/) clus = clus + 1; else print \$3 "\\t" clus}' xxx.clstr | \
        sed 's/[>Contig_,...]//g' | \
        sort -g -k1 > sort.contig.cluster.ids
    
    paste sort.contig.cluster.ids ${totaluniqseq} > contig.cluster.totaluniqseq
    
    sort -k2,2 -g contig.cluster.totaluniqseq | \
        sed -e 's/NNNNNNNNNN/\\t/g' > rcluster
    
    # Rainbow div - redivide clusters
    echo "Running Rainbow div..."
    rainbow div -i rcluster \
                -o rbdiv.out \
                -f ${div_f} \
                -K ${div_K}
    
    # Rainbow merge - merge redivided clusters
    echo "Running Rainbow merge..."
    rainbow merge -o rbasm.out \
                  -a \
                  -i rbdiv.out \
                  -r ${merge_r}
    
    # Extract optimal sequences for de novo genome
    echo "Extracting optimal sequences..."
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
    
    # Final clustering to remove redundancy
    echo "Final CD-HIT clustering to remove redundancy..."
    cd-hit-est -i rainbow.fasta \
               -o denovo_reference.fa \
               -M 0 \
               -T ${task.cpus} \
               -c ${final_similarity}
    
    # Generate assembly statistics
    echo "Generating assembly statistics..."
    
    # Calculate statistics for final assembly
    cat > assembly_stats.txt <<-EOF
	Rainbow Assembly Statistics
	===========================
	Assembly Parameters:
	  Initial clustering: ${cluster_similarity}
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
    
    # Calculate N50 and other metrics
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
    
    # Contig size distribution
    echo "" >> assembly_stats.txt
    echo "Contig Size Distribution:" >> assembly_stats.txt
    echo "  >10kb: \$(awk '\$1>10000' contig_lengths.txt | wc -l)" >> assembly_stats.txt
    echo "  5-10kb: \$(awk '\$1>=5000 && \$1<=10000' contig_lengths.txt | wc -l)" >> assembly_stats.txt
    echo "  1-5kb: \$(awk '\$1>=1000 && \$1<5000' contig_lengths.txt | wc -l)" >> assembly_stats.txt
    echo "  500bp-1kb: \$(awk '\$1>=500 && \$1<1000' contig_lengths.txt | wc -l)" >> assembly_stats.txt
    echo "  <500bp: \$(awk '\$1<500' contig_lengths.txt | wc -l)" >> assembly_stats.txt
    
    cat assembly_stats.txt
    
    # Cleanup
    rm -f contig_lengths.txt xxx xxx.clstr
    """
}