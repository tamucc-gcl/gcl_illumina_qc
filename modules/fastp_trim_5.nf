process fastp_trim_5 {
  label 'fastp'
  tag "$sample_id"

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}_trim5_1.fq.gz"), path("${sample_id}_trim5_2.fq.gz") emit: reads
    path("${sample_id}_trim5_fastp.json") emit: json
    path("${sample_id}_trim5_fastp.html") emit: html

  script:
  """
  fastp -i ${reads[0]} -I ${reads[1]} \
        -o ${sample_id}_trim5_1.fq.gz -O ${sample_id}_trim5_2.fq.gz \
        --trim_front1 5 --trim_front2 5 -w 4 \
        --json ${sample_id}_trim5_fastp.json \
        --html ${sample_id}_trim5_fastp.html
  """
}
