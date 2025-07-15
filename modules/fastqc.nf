
process fastqc_raw {
  label 'fastqc'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id),
        path("${sample_id}_R1_fastqc.html"),
        path("${sample_id}_R1_fastqc.zip"),
        path("${sample_id}_R2_fastqc.html"),
        path("${sample_id}_R2_fastqc.zip")

  script:
  """
    fastqc -o . ${reads[0]} ${reads[1]}
    # Rename outputs for consistent downstream naming
    mv $(basename ${reads[0]} .fq.gz)_fastqc.html ${sample_id}_R1_fastqc.html
    mv $(basename ${reads[0]} .fq.gz)_fastqc.zip  ${sample_id}_R1_fastqc.zip
    mv $(basename ${reads[1]} .fq.gz)_fastqc.html ${sample_id}_R2_fastqc.html
    mv $(basename ${reads[1]} .fq.gz)_fastqc.zip  ${sample_id}_R2_fastqc.zip
  """
}
