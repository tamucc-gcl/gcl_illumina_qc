// modules/rank_and_select.nf  (chunk 5b)
// Single weight-free rank aggregation over ALL candidate signals (1-pass; no
// gate). Concatenates the per-candidate cheap rows + SNP rows, runs the rank R
// script, and emits the ranked table, the winning id, and a signals plot.

process rank_and_select {
    label 'optimize_rank'
    tag "rank_and_select"

    publishDir "${params.outdir}/denovo_assembly/optimize", mode: params.publish_dir_mode, pattern: "{final_rank.tsv,optimize_plot.png}"

    input:
        path(cheap_rows)        // many cheap_<id>.tsv (id c1 c2 sim n_contigs total_len mean_len)
        path(snp_rows)          // many snp_<id>.tsv   (id c1 c2 sim concordance r80_loci snps_per_locus n_snps n_called_contigs)
        path(nb_cutoff1)        // nb_cutoff1.value
        val(expected_loci)      // integer or "NA"

    output:
        path("final_rank.tsv"),   emit: final_rank
        path("best_id.value"),    emit: best_id
        path("optimize_plot.png"), emit: plot

    script:
    """
    set -euo pipefail
    cat cheap_*.tsv > cheap_all.tsv
    cat snp_*.tsv   > snp_all.tsv
    echo "Aggregating \$(wc -l < cheap_all.tsv) candidates (\$(wc -l < snp_all.tsv) with SNP rows)"

    Rscript ${projectDir}/r_scripts/rank_and_select.R \\
        cheap_all.tsv snp_all.tsv ${nb_cutoff1} ${expected_loci}
    """
}
