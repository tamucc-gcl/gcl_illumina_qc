// main.nf

nextflow.enable.dsl=2

params.reads = "data/*_{1,2}.fq.gz"
params.accession = "GCF_002234675.1"
params.outdir = "results"
params.multiqc_dir = "results/multiqc"

workflow {
  Channel.fromFilePairs(params.reads, flat: true) 
         .set { read_pairs }

  // Run MultiQC on raw data
  read_pairs.collect().set { raw_reads }
  multiqc(raw_reads, "raw", params.multiqc_dir)

  // Download and prepare genome
  prepare_genome(params.accession)
    .set { genome_info }

  fastp_trim_3(read_pairs, outdir: "data/fq_fp1")
    .into { trimmed_reads_for_clumpify; fastp3_out }

  trimmed_reads_for_clumpify
    .combine(Channel.value("data/fq_fp1_clmp"))
    | clumpify
    .into { clumped_reads_for_fp5; clumpify_out }

  clumped_reads_for_fp5
    .combine(Channel.value("data/fq_fp1_clmp_fp5"))
    | fastp_trim_5
    .into { fp5_reads_for_screen; fastp5_out }

  fp5_reads_for_screen
    .combine(Channel.value("data/fq_fp1_clmp_fp5_scrn"))
    | fastq_screen
    .into { screened_reads_for_repair; fastqscreen_out }

  screened_reads_for_repair
    .combine(Channel.value("data/fq_fp1_clmp_fp5_scrn_rpr"))
    | repair
    .into { repaired_reads_for_mapping; repair_out }

  repaired_reads_for_mapping
    .combine(genome_info)
    .combine(Channel.value("data/fq_fp1_clmp_fp5_scrn_rpr_map"))
    | map_reads

  multiqc(fastp3_out, "fastp_trim_3", params.multiqc_dir)
  multiqc(clumpify_out, "clumpify", params.multiqc_dir)
  multiqc(fastp5_out, "fastp_trim_5", params.multiqc_dir)
  multiqc(fastqscreen_out, "fastq_screen", params.multiqc_dir)
  multiqc(repair_out, "repair", params.multiqc_dir)
}

workflow prepare_genome {
  take: accession
  main:
    fetch_genome(accession) | index_genome
  emit:
    index_genome.out
}

include { fastp_trim_3 } from './modules/fastp_trim_3.nf'
include { clumpify } from './modules/clumpify.nf'
include { fastp_trim_5 } from './modules/fastp_trim_5.nf'
include { fastq_screen } from './modules/fastq_screen.nf'
include { repair } from './modules/repair.nf'
include { map_reads } from './modules/map_reads.nf'
include { multiqc } from './modules/multiqc.nf'
include { fetch_genome } from './modules/fetch_genome.nf'
include { index_genome } from './modules/index_genome.nf'