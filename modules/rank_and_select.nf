// modules/rank_and_select.nf  (PIVOT: r80-vs-n_contigs elbow selector)
// Concatenates per-candidate cheap rows + SNP rows and runs the elbow selector.
// Winner = fewest contigs at the r80 plateau (Kneedle on the r80-vs-size envelope).

process rank_and_select {
    label 'optimize_rank'
    tag "rank_and_select"

    publishDir "${params.outdir}/denovo_assembly/optimize", mode: params.publish_dir_mode, pattern: "{final_rank.tsv,optimize_plot.png,optimize_params_plot.png}"

    input:
        path(cheap_rows)        // many cheap_<id>.tsv (id c1 c2 sim n_contigs total_len mean_len)
        path(snp_rows)          // many snp_<id>.tsv   (id c1 c2 sim concordance r80_loci snps_per_locus n_snps n_called_contigs)
        path(nb_cutoff1)        // nb_cutoff1.value (reported)
        val(expected_loci)      // integer or "NA" (reported as anchor)
        val(knee_frac)          // leveling-off threshold: fraction of initial fitted slope (default 0.10)

    output:
        path("final_rank.tsv"),          emit: final_rank
        path("best_id.value"),           emit: best_id
        path("optimize_plot.png"),       emit: plot
        path("optimize_params_plot.png"), emit: params_plot

    script:
    """
    set -euo pipefail
    cat cheap_*.tsv > cheap_all.tsv
    cat snp_*.tsv   > snp_all.tsv
    echo "Selecting from \$(wc -l < cheap_all.tsv) candidates (\$(wc -l < snp_all.tsv) with SNP rows)"

    Rscript ${projectDir}/r_scripts/rank_and_select.R \\
        cheap_all.tsv snp_all.tsv ${nb_cutoff1} ${expected_loci} ${knee_frac}
    """
}
