// main.nf – GCL Illumina QC pipeline
// July 2025 — consolidated MultiQC per step

nextflow.enable.dsl = 2

//--------------------------------------------------------------------
// USER PARAMETERS
//--------------------------------------------------------------------
params.reads       = "data/fq_raw/*.{1,2}.fq.gz"    // paired‑end,  sampleID.1.fq.gz / .2.fq.gz
params.accession   = "GCA_042920385.1"              // NCBI assembly accession
params.outdir      = "results"
params.multiqc_dir = "${params.outdir}/multiqc"

//--------------------------------------------------------------------
// WORKFLOW DEFINITION
//--------------------------------------------------------------------
workflow {

    //----------------------------------------------------------------
    // 1. RAW READ INPUT
    //----------------------------------------------------------------
    Channel
        .fromFilePairs( params.reads, flat: true )
        .map { sid, reads -> tuple( sid, reads[0], reads[1] ) }
        .set { raw_reads_pairs }

    //----------------------------------------------------------------
    // 2. FASTQC ON RAW READS  ➜  MULTIQC (raw_fastqc)
    //----------------------------------------------------------------
    raw_reads_pairs | fastqc_raw

    // collect all FastQC HTML/ZIP files for MultiQC
    fastqc_raw.out
        .map{ sid, html1, html2, zip1, zip2 -> [html1, html2, zip1, zip2] }
        .flatten()
        .map { f -> tuple('raw_fastqc', params.multiqc_dir, f) }
        .set { multiqc_raw_in }

    multiqc_raw_in | multiqc            // emits multiqc_${step}_report.html


    //----------------------------------------------------------------
    // 3. GENOME PREP  (download ➜ index)
    //----------------------------------------------------------------
    Channel.value( params.accession ) | prepare_genome | set { genome_ch }

    //----------------------------------------------------------------
    // 4. READ PREP  (3′ trim → clumpify → 5′ trim → FastQ‑Screen → repair)
    //----------------------------------------------------------------
    raw_reads_pairs                
    | fastp_trim_3             
    | clumpify                 
    | fastp_trim_5             
    | fastq_screen            
    | repair                   
    | set { repaired_reads_ch }

/*
    //----------------------------------------------------------------
    // 5. MULTIQC FOR EACH QC STEP
    //----------------------------------------------------------------
    // fastp_trim_3
    fastp_trim_3.out
        .map{ sid, r1, r2, json, html -> [json, html] }
        .flatten()
        .map { f -> tuple('fastp_trim_3', params.multiqc_dir, f) }
        .set { multiqc_fastp3_in }

    // clumpify (stats file)
    clumpify.out
        .map{ sid, r1, r2, stats -> stats }
        .map { f -> tuple('clumpify', params.multiqc_dir, f) }
        .set { multiqc_clumpify_in }

    // fastp_trim_5
    fastp_trim_5.out
        .map{ sid, r1, r2, json, html -> [json, html] }
        .flatten()
        .map { f -> tuple('fastp_trim_5', params.multiqc_dir, f) }
        .set { multiqc_fastp5_in }

    // fastq_screen (txt reports)
    fastq_screen.out
        .map{ sid, r1, r2, txt1, txt2 -> [txt1, txt2] }
        .flatten()
        .map { f -> tuple('fastq_screen', params.multiqc_dir, f) }
        .set { multiqc_fastq_in }

    // repair – nothing QC‑related, but include to track reads length/GC if desired
    repair.out
        .map{ sid, r1, r2 -> [r1, r2] }
        .flatten()
        .map { f -> tuple('repair', params.multiqc_dir, f) }
        .set { multiqc_repair_in }
/*
    // Merge all QC inputs into one channel and run MultiQC once per step
    multiqc_fastp3_in.mix( multiqc_clumpify_in )
                    .mix( multiqc_fastp5_in )
                    .mix( multiqc_fastq_in )
                    .mix( multiqc_repair_in ) | multiqc
*/
/*
    fastp_trim_5.out
              .map{ sid, f1, f2, json, html -> [json, html] }
              .flatten()
              .set { fastp5_qc_ch }

    multiqc_fastp5(
        fastp5_qc_ch,
        "fastp_trim_5",
        params.multiqc_dir
    )


    fastq_screen.out
               .map{ sid, r1, r2 -> [r1, r2] }
               .flatten()
               .set { fastqscreen_files_ch }

    multiqc_fastqscreen(
        fastqscreen_files_ch,
        "fastq_screen",
        params.multiqc_dir
    )


    repair.out
               .map{ sid, r1, r2 -> [r1, r2] }
               .flatten()
               .set { repair_files_ch }

    multiqc_repair(
        repair_files_ch,
        "repair",
        params.multiqc_dir
    )
*/
}

//--------------------------------------------------------------------
// SUB‑WORKFLOW: fetch + index genome
//--------------------------------------------------------------------
workflow prepare_genome {
    take:
        accession
    main:
        fetch_genome( accession ) | index_genome
    emit:
        index_genome.out
}

//--------------------------------------------------------------------
// MODULE IMPORTS
//--------------------------------------------------------------------
include { fastp_trim_3 }      from './modules/fastp_trim_3.nf'
include { clumpify }          from './modules/clumpify.nf'
include { fastp_trim_5 }      from './modules/fastp_trim_5.nf'
include { fastq_screen }      from './modules/fastq_screen.nf'
include { repair }            from './modules/repair.nf'
include { map_reads }         from './modules/map_reads.nf'
include { fetch_genome }      from './modules/fetch_genome.nf'
include { index_genome }      from './modules/index_genome.nf'

include { fastqc_raw }        from './modules/fastqc.nf'

/*
// one MultiQC process, aliased as simply `multiqc`
include { multiqc }           from './modules/multiqc.nf'
*/