// modules/write_anchor_calc.nf
// Computes AND writes the biological-anchor (expected ddRAD locus count) using the
// rare-cutter-anchored exponential-spacing model, and emits both the expected-loci
// value (for the rank step) and a full human-readable calculation file (published).
//
// Self-contained: all arithmetic is done here in awk (exp/power/format, no deps),
// so the model is the single source of truth and is testable in isolation. The
// workflow only decides ENABLED vs DISABLED and passes raw params through.
//
// Model: a usable ddRAD locus is a RARE-cutter site (larger recognition sequence,
// cuts less often) whose NEAREST COMMON-cutter site falls within the INSERT size
// window. Common-cutter inter-site spacing ~ Exponential(rate = 1/4^common_len).
//   P(nearest common site in [min,max], one side) = exp(-rate*min) - exp(-rate*max)
//   P(window satisfied on >=1 of 2 sides)         = 1 - (1-p_one)^2
//   expected_loci = (rare-cutter sites genome-wide) * P_window
// IMPORTANT: insert window = genomic DNA between cut sites, NOT the fragment-
// analyzer trace window (which includes adapters).

process write_anchor_calc {
    label 'basic'
    tag "anchor_calc"

    publishDir "${params.outdir}/denovo_assembly/optimize", mode: params.publish_dir_mode, pattern: "expected_loci_calculation.txt"

    input:
        val(enabled)            // true => compute; false => write DISABLED note
        val(genome_size)        // bp (or null/"" when disabled)
        val(site_len_1)         // recognition-site length, enzyme 1
        val(site_len_2)         // recognition-site length, enzyme 2
        val(insert_min)         // insert window lower bound (bp)
        val(insert_max)         // insert window upper bound (bp)
        val(enzyme_pair)        // informational label (or null)

    output:
        path("expected_loci.value"),            emit: expected_loci   // integer, or "NA"
        path("expected_loci_calculation.txt"),  emit: calc

    script:
    if (enabled)
        """
        set -euo pipefail
        awk -v G='${genome_size}' -v s1='${site_len_1}' -v s2='${site_len_2}' \\
            -v imin='${insert_min}' -v imax='${insert_max}' -v epair='${enzyme_pair ?: "unspecified"}' '
        BEGIN {
            G = G + 0; s1 = s1 + 0; s2 = s2 + 0; imin = imin + 0; imax = imax + 0
            rare_len   = (s1 > s2) ? s1 : s2
            common_len = (s1 < s2) ? s1 : s2
            rare_period   = 4 ^ rare_len
            common_period = 4 ^ common_len
            n_rare   = G / rare_period
            n_common = G / common_period
            common_rate = 1.0 / common_period
            p_one = exp(-common_rate * imin) - exp(-common_rate * imax)
            p_win = 1 - (1 - p_one) ^ 2
            exp_loci = int(n_rare * p_win + 0.5)

            print exp_loci > "expected_loci.value"

            cf = "expected_loci_calculation.txt"
            printf "Expected ddRAD Locus Count — Biological Anchor (signal 2)\n" > cf
            printf "=========================================================\n" >> cf
            printf "Model: rare-cutter-anchored. A usable ddRAD locus is a RARE-cutter site\n" >> cf
            printf "whose nearest COMMON-cutter site falls within the INSERT size-selection\n" >> cf
            printf "window. Common-cutter inter-site spacing ~ Exponential(rate = 1/4^common_len).\n\n" >> cf
            printf "Inputs:\n" >> cf
            printf "  Genome size estimate (G)          : %.3e bp\n", G >> cf
            printf "  Enzyme pair                       : %s\n", epair >> cf
            printf "  Recognition-site lengths          : %d and %d bp\n", s1, s2 >> cf
            printf "  Rare cutter (larger site)         : %d-cutter -> 1 site per 4^%d = %.0f bp\n", rare_len, rare_len, rare_period >> cf
            printf "  Common cutter (smaller site)      : %d-cutter -> 1 site per 4^%d = %.0f bp\n", common_len, common_len, common_period >> cf
            printf "  Insert window (genomic, NOT trace): %.0f - %.0f bp\n\n", imin, imax >> cf
            printf "Derivation:\n" >> cf
            printf "  Rare-cutter sites genome-wide     : G / 4^%d   = %.0f\n", rare_len, n_rare >> cf
            printf "  Common-cutter sites genome-wide   : G / 4^%d   = %.0f\n", common_len, n_common >> cf
            printf "  Common-cutter rate (per bp)       : 1 / 4^%d   = %.3e\n", common_len, common_rate >> cf
            printf "  P(nearest common site in window, one side):\n" >> cf
            printf "      exp(-rate*%.0f) - exp(-rate*%.0f)        = %.5f\n", imin, imax, p_one >> cf
            printf "  P(window satisfied on >=1 of 2 sides):\n" >> cf
            printf "      1 - (1 - %.5f)^2                       = %.5f\n\n", p_one, p_win >> cf
            printf "  EXPECTED LOCI = rare_sites * P_window = %.0f * %.5f = %d\n\n", n_rare, p_win, exp_loci >> cf
            printf "Note: candidate assemblies typically contain many MORE contigs than this\n" >> cf
            printf "(repeats, paralogs, allelic/error contigs the cutoffs prune toward the true\n" >> cf
            printf "single-copy count). The anchor ranks candidates by |n_contigs - expected|;\n" >> cf
            printf "it is a soft signal, not a hard target.\n" >> cf
        }'
        cat expected_loci_calculation.txt
        echo "expected_loci = \$(cat expected_loci.value)"
        """
    else
        """
        set -euo pipefail
        echo "NA" > expected_loci.value
        cat > expected_loci_calculation.txt <<'DISABLED_EOF'
Expected ddRAD Locus Count — Biological Anchor (signal 2)
=========================================================
DISABLED. The expected-locus anchor was not computed because one or more of the
required parameters was not supplied:
  genome_size_est : ${genome_size ?: 'MISSING'}
  size_select_min : ${insert_min ?: 'MISSING'}  (INSERT window, not fragment-trace)
  size_select_max : ${insert_max ?: 'MISSING'}
Supply all three (plus optional enzyme*_site_len) to enable signal 2.
DISABLED_EOF
        cat expected_loci_calculation.txt
        """
}
