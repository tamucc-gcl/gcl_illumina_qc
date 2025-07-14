process repair {
  label 'repair'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)
  val outdir

  output:
  tuple val(sample_id), path("${outdir}/*.fq.gz")

  script:
  """
  mkdir -p ${outdir}
  repair.sh in=${reads[0]} in2=${reads[1]} out=${outdir}/${sample_id}_repaired_1.fq.gz out2=${outdir}/${sample_id}_repaired_2.fq.gz outs=stdout
  """
}
