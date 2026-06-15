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

    # Decide whether optical dedup is even possible: Clumpify needs 7 colon-fields
    # (…:lane:tile:x:y) in the read header. Old/demultiplexed reads often don't have them.
    HDR=\$(zcat ${read1} | head -n1)
    NCOLON=\$(echo "\$HDR" | tr -cd ':' | wc -c)
    if [ "${want_optical}" = "true" ] && [ "\$NCOLON" -ge 6 ]; then
        OPTICAL="optical=t dupedist=12000"
        echo "Header looks coordinate-bearing (\$NCOLON colons) -> optical dedup ON"
    else
        OPTICAL="optical=f"
        echo "Header lacks flowcell coordinates (\$NCOLON colons) or optical not requested -> optical dedup OFF"
    fi

    clumpify.sh \
        in=${read1} \
        in2=${read2} \
        out=${sample_id}_fp1-clmp.r1.fq.gz \
        out2=${sample_id}_fp1-clmp.r2.fq.gz \
        overwrite=t \
        dedupe=t \
        \${OPTICAL} \
        addcount=t \
        ${subs_param} \
        ${containment_param} \
        consensus=f \
        2>&1 | tee ${sample_id}_clumpify_stats.txt
    """
}