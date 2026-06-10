// modules/compute_cheap_signals.nf
// CHUNK 3b — per-candidate CHEAP signals for de novo optimization.
// Computes signals that are NOT monotone in reference size (or whose raw inputs
// the rank R script turns into non-monotone signals across the grid):
//   1a inflection INPUTS: n_contigs, total_len, mean_len  (curve computed across
//        the grid by the rank R script — a single candidate can't see the curve)
//   1b internal redundancy: fraction of contigs that are near-duplicates of
//        another contig, detected by self-clustering the reference at a HIGHER
//        identity than the assembly used (the reference is already a CD-HIT
//        product at its own sim, so re-clustering above that surfaces residual
//        redundancy / collapsed-paralog leftovers).
//   2 anchor INPUT: n_contigs is emitted; the rank R script computes
//        |n_contigs - expected_loci| if the biological anchor is supplied.
//
// Emits ONE tsv row per candidate. No scoring here — ranking happens globally in
// the provisional_rank R step so signals like 1a (which need the whole grid) and
// 5b (NB cutoff1, a global value) can be joined in.

process compute_cheap_signals {
    label 'denovo_assembly'      // reuses cd-hit-est (already in this env)
    tag "${meta.id}"

    input:
        tuple val(meta), path(reference)
        val(redundancy_identity)   // self-cluster identity for 1b (e.g. 0.98)

    output:
        path("cheap_${meta.id}.tsv"), emit: cheap_row

    script:
    def cdhit_n = { c -> def t = c as double
                    t >= 0.95 ? 10 : t >= 0.90 ? 8 : t >= 0.88 ? 7 :
                    t >= 0.85 ? 6  : t >= 0.80 ? 5 : 4 }
    """
    set -euo pipefail

    # ---- contig stats (1a inputs + 2 input) ----
    N_CONTIGS=\$(grep -c '^>' ${reference} || echo 0)
    # total length and mean length
    read TOTAL_LEN MEAN_LEN <<< \$(awk '/^>/{next} {l+=length(\$0); n++} END{ if(n>0) printf "%d %d", l, l/n; else printf "0 0" }' ${reference})
    # NOTE: n above counts sequence LINES; references here are 2-line records
    # (header + single seq line) so lines == contigs. Guard anyway:
    SEQ_LINES=\$(grep -vc '^>' ${reference} || echo 0)

    # ---- 1b internal redundancy ----
    # Self-cluster the reference at a higher identity than assembly. The fraction
    # of contigs ABSORBED into clusters (i.e. 1 - clusters/contigs) is the residual
    # redundancy. A clean reference self-clusters ~1:1 (low redundancy); a bloated
    # one collapses many near-duplicates (high redundancy).
    REDUN=0
    if [ "\$N_CONTIGS" -gt 1 ]; then
        cd-hit-est -i ${reference} -o selfclust -c ${redundancy_identity} \\
                   -n ${cdhit_n(redundancy_identity)} -T ${task.cpus} -M 0 -g 1 >/dev/null 2>&1 || true
        if [ -f selfclust ]; then
            N_CLUST=\$(grep -c '^>' selfclust || echo \$N_CONTIGS)
            # redundancy = 1 - (clusters / contigs); 0 = no dupes, ->1 = very redundant
            REDUN=\$(awk -v c=\$N_CLUST -v t=\$N_CONTIGS 'BEGIN{ if(t>0) printf "%.6f", 1 - (c/t); else print 0 }')
        fi
    fi

    # ---- emit one row ----
    # columns: id c1 c2 sim n_contigs total_len mean_len redundancy
    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${meta.id}" "${meta.c1}" "${meta.c2}" "${meta.sim}" \\
        "\$N_CONTIGS" "\$TOTAL_LEN" "\$MEAN_LEN" "\$REDUN" > cheap_${meta.id}.tsv

    echo "[cheap_signals ${meta.id}] contigs=\$N_CONTIGS total=\$TOTAL_LEN mean=\$MEAN_LEN redundancy=\$REDUN"
    rm -f selfclust selfclust.clstr
    """
}
