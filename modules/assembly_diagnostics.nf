// modules/assembly_diagnostics.nf
// Auto-detect de novo assembly cutoffs (cutoff1, cutoff2) from the per-sample
// unique-sequence files, and produce diagnostic curves for the report.
//
// cutoff1 = minimum within-individual coverage to keep a unique sequence
// cutoff2 = minimum number of individuals a sequence must appear in (after cutoff1)
//
// Both are chosen by knee detection on the "sequences retained vs. threshold"
// curve (see r_scripts/assembly_diagnostics.R). The emitted .value files are
// read by the denovo_assembly subworkflow and fed to filter_unique_seqs unless
// the user has explicitly set params.cutoff1 / params.cutoff2.

process assembly_diagnostics {
    label 'assembly_diagnostics'
    tag "cutoff_selection"

    publishDir "${params.outdir}/denovo_assembly/diagnostics", mode: params.publish_dir_mode

    input:
        path(uniq_seq_files)

    output:
        path("cutoff1.value"),           emit: cutoff1_value
        path("cutoff2.value"),           emit: cutoff2_value
        path("cutoff1_curve.png"),       emit: cutoff1_plot
        path("cutoff2_curve.png"),       emit: cutoff2_plot
        path("diagnostics_summary.txt"), emit: summary
        path("coverage_freq.txt"),       emit: coverage_freq

    script:
    """
    cp ${projectDir}/r_scripts/assembly_diagnostics.R .

    # ---- Cutoff 1: pooled within-individual coverage distribution ----
    echo "Building within-individual coverage distribution (cutoff1)..."
    cat *.uniq.seqs | cut -f1 | sort -n | uniq -c | sed 's/^ *//' > coverage_freq.txt
    Rscript assembly_diagnostics.R cutoff1 coverage_freq.txt
    CUTOFF1=\$(cat cutoff1.value)
    echo "Auto-selected cutoff1 = \$CUTOFF1"

    # ---- Cutoff 2: per-sequence individual counts AT the chosen cutoff1 ----
    # (mirrors the first-stage logic in filter_unique_seqs.nf)
    echo "Building per-sequence individual-count distribution (cutoff2)..."
    for f in *.uniq.seqs; do
        mawk -v x="\$CUTOFF1" '\$1 >= x' "\$f" | cut -f2
    done | perl -e 'while (<>) {chomp; \$z{\$_}++;} while((\$k,\$v) = each(%z)) {print "\$v\\t\$k\\n";}' > uniqCperindv

    cut -f1 uniqCperindv | sort -n | uniq -c | sed 's/^ *//' > indiv_freq.txt
    Rscript assembly_diagnostics.R cutoff2 indiv_freq.txt
    CUTOFF2=\$(cat cutoff2.value)
    echo "Auto-selected cutoff2 = \$CUTOFF2"

    # ---- Summary ----
    N_SAMPLES=\$(ls *.uniq.seqs | wc -l)
    cat > diagnostics_summary.txt <<EOF
De Novo Assembly Cutoff Diagnostics
===================================
Samples analyzed: \$N_SAMPLES
Auto-selected cutoff1 (min reads per individual): \$CUTOFF1
Auto-selected cutoff2 (min individuals): \$CUTOFF2

Selection method: knee detection (Kneedle-style) on the
"unique sequences retained vs. threshold" curves.
See cutoff1_curve.png and cutoff2_curve.png.

These values seed the reference assembly only; final genotype
filtering should be done downstream. Override the auto values
by setting params.cutoff1 / params.cutoff2.
EOF
    cat diagnostics_summary.txt
    """
}
