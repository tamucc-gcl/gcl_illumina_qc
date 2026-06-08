// modules/score_assemblies.nf
// Collect per-candidate score rows (c1,c2,sim), rank by composite, pick the best
// candidate id, and emit a comparison heatmap + ranked table.

process score_assemblies {
    label 'assembly_diagnostics'   // reuse the R env (tidyverse)
    tag "sweep_scoring"

    publishDir "${params.outdir}/denovo_assembly/sweep", mode: params.publish_dir_mode

    input:
        path(score_rows)

    output:
        path("best_id.value"),         emit: best_id
        path("sweep_summary.tsv"),     emit: summary
        path("sweep_comparison.png"),  emit: plot

    script:
    """
    cp ${projectDir}/r_scripts/score_assemblies.R .

    echo -e "id\\tc1\\tc2\\tsim\\tn_samples\\tmap_rate\\tpp_rate\\tsoftclip_per_read\\tmean_AS" > sweep_scores.tsv
    cat score_*.tsv >> sweep_scores.tsv

    Rscript score_assemblies.R sweep_scores.tsv
    echo "Best candidate = \$(cat best_id.value)"
    """
}
