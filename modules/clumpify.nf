process clumpify {
    label 'clumpify'
    tag   "$sample_id"

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
    def optical_param = sequencing_type == 'ddrad' ? 'optical=t' : 'optical=f'

    """
    temp_dir=\$(mktemp -d)
    trap "rm -rf \$temp_dir" EXIT

    clumpify.sh \
        in=${read1} \
        in2=${read2} \
        out=${sample_id}_fp1-clmp.r1.fq.gz \
        out2=${sample_id}_fp1-clmp.r2.fq.gz \
        overwrite=t \
        usetmpdir=t \
        deletetemp=t \
        dedupe=t \
        ${optical_param} \
        dupedist=12000 \
        addcount=t \
        ${subs_param} \
        ${containment_param} \
        consensus=f \
        2>&1 | tee ${sample_id}_clumpify_stats.txt
    """
}