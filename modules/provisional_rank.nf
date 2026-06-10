// modules/provisional_rank.nf  (chunk 3c)
// Weight-free rank aggregation over cheap signals: NB cutoff1 proximity (c1 axis)
// + coverage uniformity CV & Gini (c2/similarity quality) + optional biological
// anchor. Inflection is reported but excluded from the aggregate. Emits the
// ranked table + survivors list for stage-2.

process provisional_rank {
    label 'optimize_rank'
    tag "provisional_rank"

    publishDir "${params.outdir}/denovo_assembly/optimize", mode: params.publish_dir_mode, pattern: "provisional_rank.tsv"

    input:
        path(cheap_rows)            // many cheap_<id>.tsv (id c1 c2 sim n_contigs total_len mean_len)
        path(cv_rows)               // many cv_<id>.tsv    (id c1 c2 sim cv gini)
        path(nb_cutoff1)            // nb_cutoff1.value
        val(expected_loci)          // integer or "NA"
        val(min_surv)               // survivors to emit

    output:
        path("provisional_rank.tsv"), emit: provisional
        path("survivors.txt"),        emit: survivors

    script:
    """
    set -euo pipefail
    cat cheap_*.tsv > cheap_all.tsv
    cat cv_*.tsv    > cv_all.tsv
    echo "Aggregating \$(wc -l < cheap_all.tsv) candidates (\$(wc -l < cv_all.tsv) with CV rows)"

    Rscript ${projectDir}/r_scripts/provisional_rank.R \\
        cheap_all.tsv cv_all.tsv ${nb_cutoff1} ${expected_loci} ${min_surv}
    """
}
