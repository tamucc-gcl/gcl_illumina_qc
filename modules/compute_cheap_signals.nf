// modules/compute_cheap_signals.nf
// Per-candidate CONTIG COUNT for the selection curve. Emits n_contigs only:
//   - n_contigs is the x-axis of the r80-vs-n_contigs curve (rank_and_select.R),
//     i.e. the "assembly size" against which the diminishing-returns knee is found.
//   - n_contigs also feeds the optional biological anchor (expected-loci proximity).
// No mapping here; this is the cheap (grep-only) per-candidate signal.
//
// Full per-candidate contig-length stats (total bases, N50, mean, size distribution)
// already live in assemble_rainbow_candidate's assembly_stats.txt, so they are NOT
// duplicated here.
//
// Emits ONE tsv row per candidate:
//   id  c1  c2  isim  divf  mr  fsim  n_contigs

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

    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${meta.id}" "${meta.c1}" "${meta.c2}" "${meta.isim}" "${meta.divf}" "${meta.mr}" "${meta.fsim}" \\
        "\$N_CONTIGS" > cheap_${meta.id}.tsv

    echo "[cheap_signals ${meta.id}] contigs=\$N_CONTIGS"
    """
}
