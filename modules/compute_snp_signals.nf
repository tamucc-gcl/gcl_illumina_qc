// modules/compute_snp_signals.nf  (chunk 5b — the non-circular genotype signals)
// Runs on EVERY candidate (1-pass design; no cheap-signal gate). Two bcftools
// passes per candidate, both mapping reads to that candidate then calling:
//
//  PASS A — pseudo-rep CONCORDANCE (primary): for each pseudo-rep individual, map
//    its two 50/50 half-replicates SEPARATELY, call genotypes on each, and compute
//    DEPTH-WEIGHTED genotype agreement at sites covered in both halves:
//        concordance = sum_i w_i * [gt_a_i == gt_b_i] / sum_i w_i,  w_i = dp_a+dp_b
//    Averaged across the pseudo-rep individuals. A good reference genotypes the
//    same individual consistently regardless of which read-half; a paralog-
//    collapsed reference yields spurious hets that disagree between halves.
//    This is NON-CIRCULAR: not monotone in reference size.
//
//  PASS B — r80 + SNP density (secondary): map the snp_sample_pct subset, joint-
//    call, then:
//      r80_loci       = # contigs with >=1 SNP genotyped in >= r80_threshold of samples
//                       (STACKS r80 rule; size-robust — junk loci fail the share test)
//      snps_per_locus = total SNPs / n_contigs (density)
//
// Emits one row: id c1 c2 sim concordance r80_loci snps_per_locus n_snps n_called_contigs

process compute_snp_signals {
    label 'stage2_call'      // bwa-mem2 + samtools + bcftools
    tag "${meta.id}"

    publishDir "${params.outdir}/denovo_assembly/optimize/snp_signals", mode: params.publish_dir_mode, pattern: "snp_${meta.id}.tsv"

    input:
        // candidate reference + pseudo-rep half fastqs (flat bag) + snp-subset fastqs (flat bag)
        tuple val(meta), path(reference), path(pseudo_reads), path(snp_reads)
        val(r80_threshold)

    output:
        path("snp_${meta.id}.tsv"), emit: snp_row

    script:
    """
    set -euo pipefail
    bwa-mem2 index ${reference}
    samtools faidx ${reference}

    # ---------------------------------------------------------------
    # PASS A — pseudo-rep concordance
    # Pseudo-rep halves are named <parent>_a.r1.fq.gz / <parent>_a.r2.fq.gz and
    # <parent>_b.r1.fq.gz / <parent>_b.r2.fq.gz. Map each half, call, compare.
    # ---------------------------------------------------------------
    CONC_NUM=0; CONC_DEN=0; N_INDIV=0

    # discover pseudo-rep parents by stripping _a.r1 suffix from a-half R1 files
    for ar1 in \$(ls *_a.r1.fq.gz 2>/dev/null | sort); do
        parent=\${ar1%_a.r1.fq.gz}
        ar2=\${parent}_a.r2.fq.gz
        br1=\${parent}_b.r1.fq.gz
        br2=\${parent}_b.r2.fq.gz
        [ -f "\$br1" ] || continue

        for half in a b; do
            r1=\${parent}_\${half}.r1.fq.gz; r2=\${parent}_\${half}.r2.fq.gz
            bwa-mem2 mem -t ${task.cpus} ${reference} "\$r1" "\$r2" 2>/dev/null \\
                | samtools sort -@ ${task.cpus} -o \${half}.bam -
            samtools index \${half}.bam
            bcftools mpileup -f ${reference} -a FORMAT/DP -Ou \${half}.bam 2>/dev/null \\
                | bcftools call -mv -Oz -o \${half}.vcf.gz 2>/dev/null
            bcftools index \${half}.vcf.gz
            # site -> GT,DP table
            bcftools query -f '%CHROM\\t%POS[\\t%GT\\t%DP]\\n' \${half}.vcf.gz > \${half}.gt
        done

        # join a.gt and b.gt on CHROM,POS; depth-weighted agreement
        read NUM DEN <<< \$(awk '
            FNR==NR { key=\$1"_"\$2; gtA[key]=\$3; dpA[key]=(\$4=="."?0:\$4); next }
            { key=\$1"_"\$2; if (key in gtA) {
                  dp = dpA[key] + (\$4=="."?0:\$4)
                  agree = (gtA[key]==\$3) ? 1 : 0
                  num += dp*agree; den += dp
              } }
            END { printf "%.6f %.6f", num, den }' a.gt b.gt)
        if awk "BEGIN{exit !(\$DEN>0)}"; then
            CONC_NUM=\$(awk -v a=\$CONC_NUM -v b=\$NUM 'BEGIN{print a+b}')
            CONC_DEN=\$(awk -v a=\$CONC_DEN -v b=\$DEN 'BEGIN{print a+b}')
            N_INDIV=\$((N_INDIV+1))
        fi
        rm -f a.bam b.bam a.bam.bai b.bam.bai a.vcf.gz b.vcf.gz a.vcf.gz.csi b.vcf.gz.csi a.gt b.gt
    done

    if awk "BEGIN{exit !(\$CONC_DEN>0)}"; then
        CONCORDANCE=\$(awk -v n=\$CONC_NUM -v d=\$CONC_DEN 'BEGIN{printf "%.6f", n/d}')
    else
        CONCORDANCE=NA
    fi

    # ---------------------------------------------------------------
    # PASS B — r80 loci + SNP density on the snp-subset
    # ---------------------------------------------------------------
    SNP_BAMS=()
    for r1 in \$(ls *.r1.fq.gz 2>/dev/null | grep -v '_a.r1.fq.gz' | grep -v '_b.r1.fq.gz' | sort); do
        r2=\${r1/.r1./.r2.}
        [ -f "\$r2" ] || continue
        bn=\$(basename "\$r1" .r1.fq.gz)
        bwa-mem2 mem -t ${task.cpus} ${reference} "\$r1" "\$r2" 2>/dev/null \\
            | samtools sort -@ ${task.cpus} -o \${bn}.snp.bam -
        samtools index \${bn}.snp.bam
        SNP_BAMS+=(\${bn}.snp.bam)
    done

    N_CONTIGS=\$(grep -c '^>' ${reference} || echo 0)

    if [ \${#SNP_BAMS[@]} -ge 2 ]; then
        N_SAMPLES=\${#SNP_BAMS[@]}
        bcftools mpileup -f ${reference} -a FORMAT/DP -Ou \${SNP_BAMS[@]} 2>/dev/null \\
            | bcftools call -mv -Oz -o snps.vcf.gz 2>/dev/null
        bcftools index snps.vcf.gz

        N_SNPS=\$(bcftools view -H snps.vcf.gz 2>/dev/null | wc -l)

        # r80: per SNP, fraction of samples with a non-missing GT; a LOCUS (contig)
        # counts if it has >=1 SNP meeting the r80 share threshold. Count distinct
        # contigs passing.
        bcftools query -f '%CHROM\\t[%GT,]\\n' snps.vcf.gz 2>/dev/null \\
          | awk -v thr=${r80_threshold} -v ns=\$N_SAMPLES '
              {
                n=split(\$2, g, ",")
                called=0
                for (i=1;i<=n;i++) if (g[i]!="" && g[i]!="./." && g[i]!=".|.") called++
                if (ns>0 && called/ns >= thr) pass_contig[\$1]=1
              }
              END { c=0; for (k in pass_contig) c++; print c }' > r80.count
        R80_LOCI=\$(cat r80.count)
    else
        N_SAMPLES=\${#SNP_BAMS[@]}
        N_SNPS=0; R80_LOCI=0
    fi

    if [ "\$N_CONTIGS" -gt 0 ]; then
        SNPS_PER_LOCUS=\$(awk -v s=\$N_SNPS -v c=\$N_CONTIGS 'BEGIN{printf "%.6f", s/c}')
    else
        SNPS_PER_LOCUS=0
    fi

    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${meta.id}" "${meta.c1}" "${meta.c2}" "${meta.fsim}" \\
        "\$CONCORDANCE" "\$R80_LOCI" "\$SNPS_PER_LOCUS" "\$N_SNPS" "\$N_CONTIGS" \\
        > snp_${meta.id}.tsv

    echo "[snp ${meta.id}] conc=\$CONCORDANCE r80=\$R80_LOCI snps/locus=\$SNPS_PER_LOCUS (snps=\$N_SNPS contigs=\$N_CONTIGS)"
    rm -f *.snp.bam *.snp.bam.bai snps.vcf.gz snps.vcf.gz.csi ${reference}.* ${reference}.fai 2>/dev/null || true
    """
}
