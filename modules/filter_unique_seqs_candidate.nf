// modules/filter_unique_seqs_candidate.nf
// Candidate variant of filter_unique_seqs: filters pooled unique sequences at a
// per-candidate (cutoff1, cutoff2) carried in a meta map, so it fans out over the
// cutoff grid. Filtering depends only on (c1,c2) — NOT similarity — so this runs
// once per (c1,c2) pair and the result is later crossed with cluster_similarity.
// Mirrors modules/filter_unique_seqs.nf logic.
//
// meta = [c1: <int>, c2: <int>, id_cutoff: "c<c1>_k<c2>"]

process filter_unique_seqs_candidate {
    label 'denovo_assembly'
    tag "filter_${meta.id_cutoff}"

    input:
        tuple val(meta), path(uniq_seq_files)

    output:
        tuple val(meta), path("uniq.fasta"), path("totaluniqseq"), path("filter_stats.txt"), emit: filtered

    script:
    """
    echo "Filtering with cutoff1=${meta.c1}, cutoff2=${meta.c2}"

    # Total unique sequences across all samples (for stats)
    cat *.uniq.seqs > all.uniq.seqs
    TOTAL_UNIQ=\$(wc -l < all.uniq.seqs)

    # Keep sequences present at least cutoff1 times per individual, then count individuals
    parallel --no-notice -j ${task.cpus} 'mawk -v x=${meta.c1} "\\\$1 >= x" {} | cut -f2' ::: *.uniq.seqs | \\
        perl -e 'while (<>) {chomp; \$z{\$_}++;} while((\$k,\$v) = each(%z)) {print "\$v\\t\$k\\n";}' > uniqCperindv
    AFTER_INDIV=\$(wc -l < uniqCperindv)

    # Keep sequences present in at least cutoff2 individuals
    mawk -v x=${meta.c2} '\$1 >= x' uniqCperindv > uniq.k.${meta.c2}.c.${meta.c1}.seqs
    AFTER_COUNT=\$(wc -l < uniq.k.${meta.c2}.c.${meta.c1}.seqs)

    # Back to fasta
    cut -f2 uniq.k.${meta.c2}.c.${meta.c1}.seqs > totaluniqseq
    mawk '{c = c + 1; print ">Contig_" c "\\n" \$1}' totaluniqseq > uniq.fasta
    FINAL=\$(grep -c '^>' uniq.fasta || echo 0)

    cat > filter_stats.txt <<EOF
Filtering Statistics (candidate ${meta.id_cutoff})
==================================================
Cutoff 1 (min reads per individual): ${meta.c1}
Cutoff 2 (min individuals): ${meta.c2}

Total unique sequences: \$TOTAL_UNIQ
Sequences after per-individual filter: \$AFTER_INDIV
Sequences after individual count filter: \$AFTER_COUNT
Final sequences for assembly: \$FINAL
EOF
    echo "Final sequences for assembly: \$FINAL"
    rm -f all.uniq.seqs
    """
}
