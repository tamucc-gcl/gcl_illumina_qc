// modules/compute_coverage_cv.nf
// CHUNK 3c — per-candidate coverage-uniformity signal (Option 2).
// Maps a small SEEDED subset of samples to ONE candidate reference and measures
// how UNIFORM the per-contig depth is. Paralog-collapsed references pile many
// genomic loci onto a few contigs -> a fat high-depth tail -> high CV / high Gini.
// A clean reference spreads reads evenly across loci -> low CV / low Gini.
//
// This is the quality signal that discriminates along the cutoff2 / similarity
// axes, which the NB-cutoff1 signal is blind to. It is the (light) mapping that
// lives in STAGE 1 (cheap signals), distinct from the heavy bcftools STAGE 2.
//
// Per-contig depth proxy: `samtools idxstats` mapped-read count per contig,
// normalized to reads-per-bp (so contig-length variation doesn't masquerade as
// depth non-uniformity). NOTE: idxstats is a cheap proxy; if CV/Gini look noisy,
// upgrade to true per-base `samtools depth` (see provisional_rank discussion).
//
// Reads are a fixed cv_sample_n-sample subset, pooled and mapped together (we only
// need aggregate per-contig depth, not per-sample). Emits one row per candidate.

process compute_coverage_cv {
    label 'map_reads'        // bwa-mem2 + samtools env
    tag "${meta.id}"

    input:
        tuple val(meta), path(reference), path(reads)   // reads = flat list of subset fastqs (r1/r2 interleaved by pairs)

    output:
        path("cv_${meta.id}.tsv"), emit: cv_row

    script:
    """
    set -euo pipefail

    bwa-mem2 index ${reference}

    # Map each R1/R2 pair in the subset, accumulate into one sorted BAM.
    # reads are staged as <sid>.r1/<sid>.r2 pairs; glob them deterministically.
    R1S=( \$(ls *.r1.fq.gz 2>/dev/null | sort) )
    if [ \${#R1S[@]} -eq 0 ]; then
        echo "[cv ${meta.id}] no subset reads staged" >&2
        printf "%s\\t%s\\t%s\\t%s\\tNA\\tNA\\n" "${meta.id}" "${meta.c1}" "${meta.c2}" "${meta.sim}" > cv_${meta.id}.tsv
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

    # Per-contig mapped-read counts (idxstats: name, length, mapped, unmapped)
    # -> reads-per-bp -> CV and Gini computed in awk over contigs with length>0.
    samtools idxstats all.bam | awk -v id="${meta.id}" -v c1="${meta.c1}" -v c2="${meta.c2}" -v sim="${meta.sim}" '
        \$1 != "*" && \$2 > 0 {
            rpb = \$3 / \$2          # reads per bp for this contig
            n++; vals[n] = rpb; sum += rpb; sumsq += rpb*rpb
        }
        END {
            if (n < 2) { printf "%s\\t%s\\t%s\\t%s\\tNA\\tNA\\n", id, c1, c2, sim; exit }
            mean = sum / n
            # CV = sd/mean (population sd)
            var = sumsq/n - mean*mean; if (var < 0) var = 0
            cv = (mean > 0) ? sqrt(var)/mean : 0
            # Gini via sorted cumulative (mean absolute difference / (2*mean))
            # sort vals ascending
            for (i = 1; i <= n; i++) for (j = i+1; j <= n; j++) if (vals[j] < vals[i]) { t = vals[i]; vals[i] = vals[j]; vals[j] = t }
            cum = 0; g = 0
            for (i = 1; i <= n; i++) { cum += vals[i]; g += cum }
            # Gini = (n + 1 - 2*sum_i(cum_i)/sum_vals) / n
            gini = (mean > 0) ? (n + 1 - 2 * g / (mean * n)) / n : 0
            printf "%s\\t%s\\t%s\\t%s\\t%.6f\\t%.6f\\n", id, c1, c2, sim, cv, gini
        }' > cv_${meta.id}.tsv

    echo "[cv ${meta.id}] \$(cat cv_${meta.id}.tsv)"
    rm -f all.sam all.bam all.bam.bai ${reference}.*
    """
}
