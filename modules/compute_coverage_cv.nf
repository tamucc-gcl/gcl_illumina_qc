// modules/compute_coverage_cv.nf  (chunk 5a — per-base depth, size-robust)
// Per-candidate coverage-uniformity signal. Maps a small SEEDED subset of samples
// to ONE candidate and measures how UNIFORM per-contig depth is. Paralog-collapsed
// references pile multiple genomic loci onto a contig -> a fat high-depth tail.
//
// CHANGED in 5a: uses TRUE per-base depth (samtools depth -a), not idxstats read
// counts. The idxstats proxy was size-confounded (dominated by how many junk
// contigs existed). Per-base mean depth per contig, MEDIAN-NORMALIZED, makes the
// statistics scale-free. Validation (synthetic): with median-normalization the
// statistics match across small/large assemblies at equal paralog level, and
// frac>Nx separates clean (2% paralog) from collapsed (15%) ~7.5x.
//
// Emits THREE statistics (prune later after inspecting on real data):
//   frac_hi  = fraction of contigs with mean depth > mult * median  (PRIMARY:
//              direct paralog-pileup count; most discriminating + size-robust)
//   gini     = Gini of per-contig mean depth (size-robust on per-base depth)
//   cv       = coefficient of variation of per-contig mean depth (weakest)
//
// Row: id c1 c2 sim frac_hi gini cv   (NA's if too few contigs)

process compute_coverage_cv {
    label 'map_reads'        // bwa-mem2 + samtools
    tag "${meta.id}"

    input:
        tuple val(meta), path(reference), path(reads)   // reads = flat list of subset fastqs
        val(pileup_mult)     // paralog-pileup threshold multiplier (e.g. 2 => >2x median)

    output:
        path("cv_${meta.id}.tsv"), emit: cv_row

    script:
    """
    set -euo pipefail

    bwa-mem2 index ${reference}

    R1S=( \$(ls *.r1.fq.gz 2>/dev/null | sort) )
    if [ \${#R1S[@]} -eq 0 ]; then
        echo "[cv ${meta.id}] no subset reads staged" >&2
        printf "%s\\t%s\\t%s\\t%s\\tNA\\tNA\\tNA\\n" "${meta.id}" "${meta.c1}" "${meta.c2}" "${meta.sim}" > cv_${meta.id}.tsv
        exit 0
    fi

    : > all.sam
    FIRST=1
    for r1 in "\${R1S[@]}"; do
        r2=\${r1/.r1./.r2.}
        if [ "\$FIRST" -eq 1 ]; then
            bwa-mem2 mem -t ${task.cpus} ${reference} "\$r1" "\$r2" > all.sam
            FIRST=0
        else
            bwa-mem2 mem -t ${task.cpus} ${reference} "\$r1" "\$r2" | samtools view -h | grep -v '^@' >> all.sam
        fi
    done

    samtools sort -@ ${task.cpus} -o all.bam all.sam
    samtools index all.bam

    # Per-base depth (-a = include zero-depth positions) piped straight into a
    # collapse to per-contig MEAN depth (never materialize the full per-base file).
    samtools depth -a all.bam \\
        | awk '{ sum[\$1]+=\$3; cnt[\$1]++ } END{ for (c in sum) if (cnt[c]>0) print sum[c]/cnt[c] }' \\
        > contig_means.txt

    # Stats on median-normalized per-contig means (awk = computation only).
    awk -v id="${meta.id}" -v c1="${meta.c1}" -v c2="${meta.c2}" -v sim="${meta.sim}" -v mult="${pileup_mult}" '
    { v[NR]=\$1; vals[NR]=\$1; n=NR; sum+=\$1; sumsq+=\$1*\$1 }
    END{
        if (n < 2) { printf "%s\\t%s\\t%s\\t%s\\tNA\\tNA\\tNA\\n", id, c1, c2, sim; exit }
        # median (sort)
        for(i=1;i<=n;i++) for(j=i+1;j<=n;j++) if(vals[j]<vals[i]){t=vals[i];vals[i]=vals[j];vals[j]=t}
        med = (n%2) ? vals[int(n/2)+1] : (vals[n/2]+vals[n/2+1])/2
        mean = sum/n; var = sumsq/n - mean*mean; if(var<0)var=0
        cv = (mean>0) ? sqrt(var)/mean : 0
        # frac of contigs with mean depth > mult*median (paralog pileups)
        hi=0; for(i=1;i<=n;i++) if(v[i] > mult*med) hi++
        frac_hi = hi/n
        # Gini on raw per-contig means
        cum=0; g=0; for(i=1;i<=n;i++){ cum+=vals[i]; g+=cum }
        gini = (mean>0) ? (n+1-2*g/(mean*n))/n : 0
        printf "%s\\t%s\\t%s\\t%s\\t%.6f\\t%.6f\\t%.6f\\n", id, c1, c2, sim, frac_hi, gini, cv
    }' contig_means.txt > cv_${meta.id}.tsv

    echo "[cv ${meta.id}] \$(cat cv_${meta.id}.tsv)"
    rm -f all.sam all.bam all.bam.bai contig_means.txt ${reference}.*
    """
}
