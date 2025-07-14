process repair {
  label 'repair'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("data/fq_fp1_clmp_fp2_scrn_rpr/*.fq.gz")

  script:
  """
  mkdir -p data/fq_fp1_clmp_fp2_scrn_rpr
  repair.sh in=${reads[0]} in2=${reads[1]} out=data/fq_fp1_clmp_fp2_scrn_rpr/${sample_id}_repaired_1.fq.gz out2=data/fq_fp1_clmp_fp2_scrn_rpr/${sample_id}_repaired_2.fq.gz outs=stdout
  """
}