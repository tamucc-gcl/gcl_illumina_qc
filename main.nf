/*
 * main.nf – Illumina QC & mapping pipeline (consolidated MultiQC)
 */

nextflow.enable.dsl = 2

//--------------------------------------------------------------------
// PARAMETERS
//--------------------------------------------------------------------
params.reads       = "data/fq_raw/*_{1,2}.fq.gz"
params.accession   = "GCA_042920385.1"
params.outdir      = "results"
params.multiqc_dir = "${params.outdir}/multiqc"

//--------------------------------------------------------------------
// WORKFLOW
//--------------------------------------------------------------------
workflow {

    //---------------------------------------------------------------
    // 1. INPUT  (raw read pairs)
    //---------------------------------------------------------------
    Channel.fromFilePairs( params.reads )
           .set { raw_reads_pairs }

    //---------------------------------------------------------------
    // 2. RAW READ QC  (FastQC → MultiQC)
    //---------------------------------------------------------------
    raw_reads_pairs | fastqc_raw

    // build channel:  step label, outdir, and each FastQC report file
    fastqc_raw.out
        .flatMap { sid, html1, html2, zip1, zip2 ->
            [
              ['raw_fastqc', params.multiqc_dir, html1],
              ['raw_fastqc', params.multiqc_dir, html2],
              ['raw_fastqc', params.multiqc_dir, zip1 ],
              ['raw_fastqc', params.multiqc_dir, zip2 ]
            ]
        }
        .set { multiqc_raw_ch }

    // consolidated MultiQC report
    multiqc_raw_ch | multiqc

    //---------------------------------------------------------------
    // 3. GENOME PREPARATION  (download → index)
    //---------------------------------------------------------------
    Channel.value( params.accession ) | fetch_genome | index_genome | set { genome_ch }

    //---------------------------------------------------------------
    // 4. READ‑PREP CHAIN
    //---------------------------------------------------------------
    raw_reads_pairs         | fastp_trim_3         | clumpify         | fastp_trim_5         | fastq_screen         | repair         | set { cleaned_reads_ch }

    //---------------------------------------------------------------
    // 5. MAPPING
    //---------------------------------------------------------------
    map_reads( cleaned_reads_ch, genome_ch )
}

//--------------------------------------------------------------------
// MODULE IMPORTS
//--------------------------------------------------------------------
include { fastqc_raw }     from './modules/fastqc.nf'
include { fastp_trim_3 }   from './modules/fastp_trim_3.nf'
include { clumpify }       from './modules/clumpify.nf'
include { fastp_trim_5 }   from './modules/fastp_trim_5.nf'
include { fastq_screen }   from './modules/fastq_screen.nf'
include { repair }         from './modules/repair.nf'
include { fetch_genome }   from './modules/fetch_genome.nf'
include { index_genome }   from './modules/index_genome.nf'
include { map_reads }      from './modules/map_reads.nf'
include { multiqc }        from './modules/multiqc.nf'