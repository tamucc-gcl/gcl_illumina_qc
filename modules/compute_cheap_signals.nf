// modules/compute_cheap_signals.nf  (chunk 3c — redundancy removed)
// Per-candidate CONTIG STATS used by the provisional rank step:
//   - n_contigs, total_len, mean_len  ->  inflection (1a), REPORT-ONLY (excluded
//        from the aggregate rank because it duplicates NB on the c1 axis), and
//        n_contigs feeds the optional biological anchor (signal 2).
// No mapping here. Redundancy (1b) was REMOVED in chunk 3c: self-clustering a
// reference that is already a CD-HIT product returns ~0 for every candidate, so
// it was a dead signal. Assembly QUALITY is now measured by coverage-uniformity
// (compute_coverage_cv.nf: CV + Gini), which discriminates along c2/similarity.
//
// Emits ONE tsv row per candidate:
//   id  c1  c2  sim  n_contigs  total_len  mean_len

process compute_cheap_signals {
    label 'basic'
    tag "${meta.id}"

    input:
        tuple val(meta), path(reference)

    output:
        path("cheap_${meta.id}.tsv"), emit: cheap_row

    script:
    """
    set -euo pipefail
    N_CONTIGS=\$(grep -c '^>' ${reference} || echo 0)
    read TOTAL_LEN MEAN_LEN <<< \$(awk '/^>/{next} {l+=length(\$0); n++} END{ if(n>0) printf "%d %d", l, l/n; else printf "0 0" }' ${reference})

    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${meta.id}" "${meta.c1}" "${meta.c2}" "${meta.sim}" \\
        "\$N_CONTIGS" "\$TOTAL_LEN" "\$MEAN_LEN" > cheap_${meta.id}.tsv

    echo "[cheap_signals ${meta.id}] contigs=\$N_CONTIGS total=\$TOTAL_LEN mean=\$MEAN_LEN"
    """
}
