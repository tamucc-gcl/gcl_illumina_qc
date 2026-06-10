// modules/provisional_rank.nf
// CHUNK 3b — provisional (cheap-signal) ranking via weight-free rank aggregation.
// Concatenates per-candidate cheap rows, runs provisional_rank.R, emits the
// ranked table + a survivors list for the expensive stage-2 step.
//
// expected_loci: integer expected RAD locus count (signal 2 anchor), or "NA" to
// disable. min_surv: number of survivors to hand to stage-2 (chunk 4 will swap
// this fixed N for the gap-based rule).

process provisional_rank {
    label 'optimize_rank'
    tag "provisional_rank"

    input:
        path(cheap_rows)            // many cheap_<id>.tsv
        path(nb_cutoff1)            // nb_cutoff1.value
        val(expected_loci)          // integer or "NA"
        val(min_surv)               // survivors to emit

    output:
        path("provisional_rank.tsv"), emit: provisional
        path("survivors.txt"),        emit: survivors

    script:
    """
    set -euo pipefail
    # Concatenate all per-candidate rows into one headerless table
    cat cheap_*.tsv > cheap_all.tsv
    echo "Aggregating \$(wc -l < cheap_all.tsv) candidate rows"

    Rscript ${projectDir}/r_scripts/provisional_rank.R \\
        cheap_all.tsv ${nb_cutoff1} ${expected_loci} ${min_surv}
    """
}
