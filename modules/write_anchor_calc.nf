// modules/write_anchor_calc.nf
// Computes AND writes the biological-anchor (expected ddRAD locus count) using the
// rare-cutter-anchored exponential-spacing model. Emits both expected_loci.value
// (consumed by the rank step) and a published human-readable calculation file.
//
// Design: awk does ONLY the arithmetic and prints the numeric intermediates to a
// small file (no string-literal newlines in awk -> no multi-layer escaping bugs).
// The human-readable report is then assembled with a plain bash heredoc, where
// newlines are literal and need no escaping.
//
// Model: a usable ddRAD locus is a RARE-cutter site (larger recognition sequence,
// cuts less often) whose NEAREST COMMON-cutter site falls within the INSERT size
// window. Common-cutter inter-site spacing ~ Exponential(rate = 1/4^common_len).
//   P(nearest common site in [min,max], one side) = exp(-rate*min) - exp(-rate*max)
//   P(window satisfied on >=1 of 2 sides)         = 1 - (1-p_one)^2
//   expected_loci = (rare-cutter sites genome-wide) * P_window
// INSERT window = genomic DNA between cut sites, NOT the fragment-analyzer trace.

process write_anchor_calc {
    label 'basic'
    tag "anchor_calc"

    publishDir "${params.outdir}/denovo_assembly/optimize", mode: params.publish_dir_mode, pattern: "expected_loci_calculation.txt"

    input:
        val(enabled)
        val(genome_size)
        val(site_len_1)
        val(site_len_2)
        val(insert_min)
        val(insert_max)
        val(enzyme_pair)

    output:
        path("expected_loci.value"),            emit: expected_loci
        path("expected_loci_calculation.txt"),  emit: calc

    script:
    if (enabled)
        """
        set -euo pipefail

        # --- awk: arithmetic ONLY. Print key=value lines; no newline escapes. ---
        awk -v G='${genome_size}' -v s1='${site_len_1}' -v s2='${site_len_2}' \\
            -v imin='${insert_min}' -v imax='${insert_max}' '
        BEGIN {
            G += 0; s1 += 0; s2 += 0; imin += 0; imax += 0
            rare_len      = (s1 > s2) ? s1 : s2
            common_len    = (s1 < s2) ? s1 : s2
            rare_period   = 4 ^ rare_len
            common_period = 4 ^ common_len
            n_rare        = G / rare_period
            n_common      = G / common_period
            common_rate   = 1.0 / common_period
            p_one         = exp(-common_rate * imin) - exp(-common_rate * imax)
            p_win         = 1 - (1 - p_one) ^ 2
            exp_loci      = int(n_rare * p_win + 0.5)
            print "rare_len="      rare_len
            print "common_len="    common_len
            print "rare_period="   rare_period
            print "common_period=" common_period
            print "n_rare="        sprintf("%.0f", n_rare)
            print "n_common="      sprintf("%.0f", n_common)
            print "common_rate="   sprintf("%.3e", common_rate)
            print "p_one="         sprintf("%.5f", p_one)
            print "p_win="         sprintf("%.5f", p_win)
            print "exp_loci="      exp_loci
        }' > calc_vars.txt

        # load the computed values into shell variables
        while IFS='=' read -r k v; do eval "CV_\$k=\$v"; done < calc_vars.txt

        # expected_loci.value (consumed by the rank step)
        printf '%s\\n' "\$CV_exp_loci" > expected_loci.value

        # human-readable report via heredoc (literal newlines, no escaping)
        cat > expected_loci_calculation.txt <<EOF
Expected ddRAD Locus Count — Biological Anchor (signal 2)
=========================================================
Model: rare-cutter-anchored. A usable ddRAD locus is a RARE-cutter site whose
nearest COMMON-cutter site falls within the INSERT size-selection window.
Common-cutter inter-site spacing ~ Exponential(rate = 1/4^common_len).

Inputs:
  Genome size estimate (G)          : ${genome_size} bp
  Enzyme pair                       : ${enzyme_pair ?: 'unspecified'}
  Recognition-site lengths          : ${site_len_1} and ${site_len_2} bp
  Rare cutter (larger site)         : \${CV_rare_len}-cutter -> 1 site per 4^\${CV_rare_len} = \${CV_rare_period} bp
  Common cutter (smaller site)      : \${CV_common_len}-cutter -> 1 site per 4^\${CV_common_len} = \${CV_common_period} bp
  Insert window (genomic, NOT trace): ${insert_min} - ${insert_max} bp

Derivation:
  Rare-cutter sites genome-wide     : G / 4^\${CV_rare_len}   = \${CV_n_rare}
  Common-cutter sites genome-wide   : G / 4^\${CV_common_len}   = \${CV_n_common}
  Common-cutter rate (per bp)       : 1 / 4^\${CV_common_len}   = \${CV_common_rate}
  P(nearest common site in window, one side):
      exp(-rate*${insert_min}) - exp(-rate*${insert_max}) = \${CV_p_one}
  P(window satisfied on >=1 of 2 sides):
      1 - (1 - p_one)^2 = \${CV_p_win}

  EXPECTED LOCI = rare_sites * P_window = \${CV_n_rare} * \${CV_p_win} = \${CV_exp_loci}

Note: candidate assemblies typically contain many MORE contigs than this
(repeats, paralogs, allelic/error contigs the cutoffs prune toward the true
single-copy count). The anchor ranks candidates by |n_contigs - expected|;
it is a soft signal, not a hard target.
EOF

        cat expected_loci_calculation.txt
        echo "expected_loci = \$CV_exp_loci"
        """
    else
        """
        set -euo pipefail
        echo "NA" > expected_loci.value
        cat > expected_loci_calculation.txt <<EOF
Expected ddRAD Locus Count — Biological Anchor (signal 2)
=========================================================
DISABLED. The expected-locus anchor was not computed because one or more of the
required parameters was not supplied:
  genome_size_est : ${genome_size ?: 'MISSING'}
  size_select_min : ${insert_min ?: 'MISSING'}  (INSERT window, not fragment-trace)
  size_select_max : ${insert_max ?: 'MISSING'}
Supply all three (plus optional enzyme*_site_len) to enable signal 2.
EOF
        cat expected_loci_calculation.txt
        """
}
