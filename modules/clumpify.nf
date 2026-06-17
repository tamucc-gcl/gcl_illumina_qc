process clumpify {
    label 'clumpify'
    tag   "$sample_id"

    // publishDir "${params.outdir}/fqgz/clumpify", mode: params.publish_dir_mode

    input:
        tuple val(sample_id), 
              path(read1), 
              path(read2), 
              path(json),   // fastp report (ignored)
              path(html)    // fastp report (ignored)
        val(sequencing_type)

    output:
        tuple val(sample_id),
              path("${sample_id}_fp1-clmp.r1.fq.gz"),
              path("${sample_id}_fp1-clmp.r2.fq.gz"),
              path("${sample_id}_clumpify_stats.txt")

    script:
    // Set parameters based on sequencing type
    def subs_param = sequencing_type == 'ddrad' ? 'subs=0' : 'subs=2'
    def containment_param = sequencing_type == 'ddrad' ? 'containment=f' : 'containment=t'
    def want_optical = sequencing_type == 'ddrad'   // ddRAD wants optical dedup *if possible*

    """
    set -euo pipefail

    # ---- Decide the deduplication mode -------------------------------------
    # ddRAD CAVEAT: in double-digest RAD both read ends sit at fixed restriction
    # cut sites, so every molecule of a locus is byte-identical. PCR/biological
    # "duplicates" are then INDISTINGUISHABLE from independent observations, and
    # that multiplicity IS the within-individual coverage signal the de novo
    # assembly (dDocent-style) and genotyping rely on. Therefore for ddRAD we
    # must NOT remove sequence duplicates. We may still remove OPTICAL duplicates
    # (a true sequencer artifact: same cluster miscalled at adjacent x/y), but
    # only when the header carries flowcell coordinates — Clumpify needs the
    # …:lane:tile:x:y fields (>=6 colons). Old/demultiplexed reads lack them.
    #
    #   ddRAD + coordinate headers  -> dedupe=t optical=t  (strip ONLY optical artifacts)
    #   ddRAD + no coordinate headers-> dedupe=f           (clump for compression, REMOVE NOTHING)
    #   non-ddRAD (e.g. whole_genome)-> dedupe=t optical=f  (full duplicate removal — correct
    #                                   for randomly-sheared libraries)
    #
    # NOTE: the previous logic fell back to "optical=f" WITH dedupe=t when headers
    # were missing, which silently performed FULL duplicate removal on ddRAD and
    # collapsed within-individual coverage to ~1 (dead assembly). That is the bug
    # this block fixes.
    HDR=\$( { zcat ${read1} 2>/dev/null || true; } | head -n1 )
    NCOLON=\$(echo "\$HDR" | tr -cd ':' | wc -c)

    if [ "${want_optical}" = "true" ]; then
        # ddRAD: never remove sequence duplicates.
        if [ "\$NCOLON" -ge 6 ]; then
            DEDUPE_ARGS="dedupe=t optical=t dupedist=12000 addcount=t"
            echo "ddRAD + coordinate-bearing header (\$NCOLON colons) -> OPTICAL-only dedup (sequence duplicates kept as coverage)"
        else
            DEDUPE_ARGS="dedupe=f"
            echo "ddRAD + header lacks flowcell coordinates (\$NCOLON colons) -> NO deduplication (clump only); sequence duplicates kept as coverage"
        fi
    else
        # non-ddRAD (randomly-sheared libraries): full duplicate removal is correct.
        DEDUPE_ARGS="dedupe=t optical=f addcount=t"
        echo "Non-ddRAD library -> full duplicate removal (dedupe=t optical=f)"
    fi

    clumpify.sh \
        in=${read1} \
        in2=${read2} \
        out=${sample_id}_fp1-clmp.r1.fq.gz \
        out2=${sample_id}_fp1-clmp.r2.fq.gz \
        overwrite=t \
        \${DEDUPE_ARGS} \
        ${subs_param} \
        ${containment_param} \
        consensus=f \
        2>&1 | tee ${sample_id}_clumpify_stats.txt
    """
}