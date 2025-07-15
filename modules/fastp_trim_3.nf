process fastp_trim_3 {
  label 'fastp'
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id),
       path("${sample_id}_trim3_1.fq.gz"),
       path("${sample_id}_trim3_2.fq.gz"),
       path("${sample_id}_trim3_fastp.json"),
       path("${sample_id}_trim3_fastp.html")

  script:
  """
    fastp -i ${reads[0]} -I ${reads[1]} \
      -o ${sample_id}_trim3_1.fq.gz -O ${sample_id}_trim3_2.fq.gz \
      --trim_tail1 1 --trim_tail2 1 -w 4 \
      --json ${sample_id}_trim3_fastp.json \
      --html ${sample_id}_trim3_fastp.html

  """
}