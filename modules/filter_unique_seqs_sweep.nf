// modules/filter_unique_seqs_sweep.nf
// Sweep variant of filter_unique_seqs: filters the pooled unique sequences at a
// per-candidate (cutoff1, cutoff2) carried in a meta map, so it can fan out over
// a grid of cutoff combinations. Mirrors filter_unique_seqs.nf logic.
//
// meta = [c1: <int>, c2: <int>, ...]

process filter_unique_seqs_sweep {
    label 'denovo_assembly'
    tag "filter_${meta.id_cutoff}"

    input:
        tuple val(meta), path(uniq_seq_files)

    output:
        tuple val(meta), path("uniq.fasta"), path("totaluniqseq"), emit: filtered

    script:
    """
    echo "Filtering with cutoff1=${meta.c1}, cutoff2=${meta.c2}"

    # Keep sequences present at least cutoff1 times per individual, count individuals
    parallel --no-notice -j ${task.cpus} 'mawk -v x=${meta.c1} "\\\$1 >= x" {} | cut -f2' ::: *.uniq.seqs | \\
        perl -e 'while (<>) {chomp; \$z{\$_}++;} while((\$k,\$v) = each(%z)) {print "\$v\\t\$k\\n";}' > uniqCperindv

    # Keep sequences present in at least cutoff2 individuals
    mawk -v x=${meta.c2} '\$1 >= x' uniqCperindv > uniq.k.${meta.c2}.c.${meta.c1}.seqs

    # Back to fasta
    cut -f2 uniq.k.${meta.c2}.c.${meta.c1}.seqs > totaluniqseq
    mawk '{c = c + 1; print ">Contig_" c "\\n" \$1}' totaluniqseq > uniq.fasta

    echo "Final sequences for assembly: \$(grep -c '^>' uniq.fasta)"
    """
}
