// modules/split_pseudo_rep.nf
// Split one individual's repaired paired-end reads into two ~equal halves
// ("pseudo-replicates") for genotype-concordance scoring during de novo
// optimization. These halves are used ONLY inside the optimization subworkflow
// (stage-2 concordance); they are never published and never enter the production
// sample set — the intact individual flows separately to final output/mapping.
//
// Splitting is deterministic given params.sweep_seed so -resume stays cache-stable.
// Read pairs are kept together (R1/R2 of a pair go to the same half) by assigning
// halves on the read-pair index, so mate pairing is preserved in each pseudo-rep.

process split_pseudo_rep {
    label 'split_pseudo_rep'
    tag "${sample_id}"

    // No publishDir: pseudo-replicates are internal to optimization only.

    input:
        tuple val(sample_id), path(read1), path(read2)
        val(seed)

    output:
        // Two pseudo-replicate pairs, tagged _a / _b, carrying the parent id in meta
        tuple val(sample_id), val("${sample_id}_a"),
              path("${sample_id}_a.r1.fq.gz"), path("${sample_id}_a.r2.fq.gz"),
              val("${sample_id}_b"),
              path("${sample_id}_b.r1.fq.gz"), path("${sample_id}_b.r2.fq.gz"),
              emit: pseudo_reps

    script:
    """
    set -euo pipefail

    # Deterministic 50/50 split on read-pair index, mates kept together.
    # paste joins the 4 lines of an R1 record with the 4 lines of its R2 mate
    # into one TSV row; awk assigns whole pairs to half a or b using a seeded
    # multiplicative hash of the record number, so the partition is well-mixed
    # AND reproducible across -resume (no reliance on rand()).
    paste <(zcat ${read1} | paste - - - -) <(zcat ${read2} | paste - - - -) \\
      | awk -v seed=${seed} '
          {
            # 64-bit-ish LCG-style hash of (record index XOR seed); take low bit.
            x = (NR * 2654435761 + seed * 40503) % 2147483647
            if (x % 2 == 0) { print > "a.tsv" } else { print > "b.tsv" }
          }'

    for half in a b; do
        # Split the joined TSV back into R1 (cols 1-4) and R2 (cols 5-8) fastq
        if [ -s \${half}.tsv ]; then
            cut -f1-4 \${half}.tsv | tr '\\t' '\\n' | gzip > ${sample_id}_\${half}.r1.fq.gz
            cut -f5-8 \${half}.tsv | tr '\\t' '\\n' | gzip > ${sample_id}_\${half}.r2.fq.gz
        else
            # Empty half (shouldn't happen with real data) -> emit empty gzip so
            # the channel cardinality stays correct and downstream can skip it.
            : | gzip > ${sample_id}_\${half}.r1.fq.gz
            : | gzip > ${sample_id}_\${half}.r2.fq.gz
        fi
    done

    rm -f a.tsv b.tsv
    """
}
