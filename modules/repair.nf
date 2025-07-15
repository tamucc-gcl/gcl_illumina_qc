process repair {
  label 'repair'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("${sample_id}_repaired_*.fq.gz")

  script:
  """
  repair.sh \
    in1=${reads[0]} \
    in2=${reads[1]} \
    out1=${sample_id}_repaired_1.fq.gz \
    out2=${sample_id}_repaired_2.fq.gz \
    outs=stdout \
    overwrite=t \
    repair
  """
}