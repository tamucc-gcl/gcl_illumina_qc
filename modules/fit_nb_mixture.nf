// modules/fit_nb_mixture.nf
// CHUNK 3b — GLOBAL negative-binomial mixture fit for a principled cutoff1.
// Fits a 2-component NB mixture to the pooled within-individual coverage
// distribution (coverage_freq.txt, already built by assembly_diagnostics):
//   component 1 = low-coverage "error/noise" sequences
//   component 2 = real-locus coverage
// The NB cutoff1 is the coverage at which the true-locus component becomes more
// probable than the error component (posterior crossover). This is a data-driven
// alternative to the geometric knee. It is used as a RANKING signal in
// provisional_rank (candidates whose c1 is near this value rank better) — it does
// NOT change grid construction.
//
// Emits a single value file. On non-convergence/insufficient data, falls back to
// the geometric knee value (cutoff1.value) so the signal degrades gracefully.

process fit_nb_mixture {
    label 'optimize_rank'    // light R env (r-base + tidyverse); add if not present
    tag "nb_mixture_cutoff1"

    publishDir "${params.outdir}/denovo_assembly/optimize", mode: params.publish_dir_mode, pattern: "nb_mixture_fit.txt"

    input:
        path(coverage_freq)      // coverage_freq.txt: "<count> <coverage>" per line
        path(knee_fallback)      // cutoff1.value from assembly_diagnostics (fallback)

    output:
        path("nb_cutoff1.value"), emit: nb_cutoff1
        path("nb_mixture_fit.txt"), emit: fit_summary

    script:
    """
    Rscript ${projectDir}/r_scripts/fit_nb_mixture.R ${coverage_freq} ${knee_fallback}
    """
}
