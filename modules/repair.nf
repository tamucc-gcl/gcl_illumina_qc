process repair {

    label 'repair'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        tuple val(sample_id),
              path("${sample_id}_repaired_1.fq.gz"),
              path("${sample_id}_repaired_2.fq.gz")

    script:
    """
    repair.sh \
      in1=${read1} \
      in2=${read2} \
      out1=${sample_id}_repaired_1.fq.gz \
      out2=${sample_id}_repaired_2.fq.gz \
      outs=stdout \
      overwrite=t \
      repair
    """
}