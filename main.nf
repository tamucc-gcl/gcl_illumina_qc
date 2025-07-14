// main.nf

nextflow.enable.dsl=2

params.reads = "data/*_{1,2}.fq.gz"
params.accession = "GCA_042920385.1"
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

  fastp_trim_3(read_pairs)
    | clumpify(outdir: "data/fq_fp1_clmp")
    | fastp_trim_5(outdir: "data/fq_fp1_clmp_fp5")
    | fastq_screen(outdir: "data/fq_fp1_clmp_fp5_scrn")
    | repair(outdir: "data/fq_fp1_clmp_fp5_scrn_rpr")
    | map_reads(genome_info)

  fastp_trim_3.out.collect().set { fastp3_out }
  multiqc(fastp3_out, "fastp_trim_3", params.multiqc_dir)

  clumpify.out.collect().set { clumpify_out }
  multiqc(clumpify_out, "clumpify", params.multiqc_dir)

  fastp_trim_5.out.collect().set { fastp5_out }
  multiqc(fastp5_out, "fastp_trim_5", params.multiqc_dir)

  fastq_screen.out.collect().set { fastqscreen_out }
  multiqc(fastqscreen_out, "fastq_screen", params.multiqc_dir)

  repair.out.collect().set { repair_out }
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
