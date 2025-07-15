process clumpify {
    label 'clumpify'
    tag   "$sample_id"

    input:
        tuple val(sample_id), 
              path(read1), 
              path(read2), 
              path(json),   // fastp report (ignored)
              path(html)    // fastp report (ignored)

    output:
        tuple val(sample_id),
              path("${sample_id}_clumped_1.fq.gz"),
              path("${sample_id}_clumped_2.fq.gz")

    script:
    """
    temp_dir=$(mktemp -d)
    clumpify.sh \
        in=${read1} in2=${read2} \
        out=${sample_id}_clumped_1.fq.gz \
        out2=${sample_id}_clumped_2.fq.gz \
        overwrite=t \
	    usetmpdir=t \
	    tmpdir=${temp_dir} \
        deletetemp=t \
        dedupe=t \
	    addcount=t \
	    subs=2 \
	    containment=t \
	    consensus=f
    """
}