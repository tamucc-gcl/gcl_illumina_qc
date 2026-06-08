// modules/score_assemblies.nf
// Collect per-candidate score rows, rank by a composite metric, pick the best
// cluster_similarity, and emit the winning value + a comparison table and plot.

process score_assemblies {
    label 'assembly_diagnostics'   // reuse the R env (tidyverse)
    tag "sweep_scoring"

    publishDir "${params.outdir}/denovo_assembly/sweep", mode: params.publish_dir_mode

    input:
        path(score_rows)

    output:
        path("best_sim.value"),        emit: best_sim
        path("sweep_summary.tsv"),     emit: summary
        path("sweep_comparison.png"),  emit: plot

    script:
    """
    cp ${projectDir}/r_scripts/score_assemblies.R .

    # Combine candidate rows (no header in score_*.tsv) into one table with a header
    echo -e "sim\\tn_samples\\tmap_rate\\tpp_rate\\tsoftclip_per_read\\tmean_AS" > sweep_scores.tsv
    cat score_*.tsv >> sweep_scores.tsv

    Rscript score_assemblies.R sweep_scores.tsv
    echo "Best cluster_similarity = \$(cat best_sim.value)"
    """
}
