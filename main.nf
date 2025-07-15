// main.nf – sequential read‑prep + per‑step MultiQC

nextflow.enable.dsl = 2

//--------------------------------------------------------------------
// USER PARAMETERS
//--------------------------------------------------------------------
params.reads       = "data/fq_raw/*.{1,2}.fq.gz"    // paired‑end files like sample.1.fq.gz & sample.2.fq.gz
params.accession   = "GCA_042920385.1"              // NCBI assembly to fetch
params.outdir      = "results"
params.multiqc_dir = "${params.outdir}/multiqc"

//--------------------------------------------------------------------
workflow {
    //----------------------------------------------------------------
    // 1. INPUT  (raw read pairs)
    //----------------------------------------------------------------
    Channel.fromFilePairs( params.reads )               // emits: [ id , [r1,r2] ]
           .set { raw_reads_pairs }

    //── quick sanity‑check (remove once happy) ─────────────────────
    raw_reads_pairs.view { "FOUND ▶️  $it" }

    //----------------------------------------------------------------
    // 2. QC on raw reads (FastQC + MultiQC)
    //----------------------------------------------------------------
    // run FastQC on each raw read pair
    raw_reads_pairs | fastqc_raw

    // consolidate raw FastQC reports into one MultiQC HTML
    // gather FastQC HTML+ZIP files into one channel for MultiQC
    fastqc_raw.out
              .map{ sid, f1, f2, f3, f4 -> [f1,f2,f3,f4] }
              .flatten()
              .set { raw_fastqc_files_ch }

    multiqc_raw(
        raw_fastqc_files_ch,
        "raw_fastqc",
        params.multiqc_dir
    )





    //----------------------------------------------------------------
    // 3. GENOME PREPARATION (download ➜ index)
    //----------------------------------------------------------------
    //----------------------------------------------------------------
    // 3. GENOME PREPARATION (download ➜ index)
    //----------------------------------------------------------------
    Channel.value(params.accession) \
        | fetch_genome \
        | index_genome \
        | set { genome_ch }


    //----------------------------------------------------------------
    // 4. READ‑PREP CHAIN: 3′trim ➜ clumpify ➜ 5′trim ➜ fastq_screen ➜ repair
    //----------------------------------------------------------------
    raw_reads_pairs \
        | fastp_trim_3 \
        | clumpify     \
        | fastp_trim_5 \
        | fastq_screen \
        | repair       \
        | set { repaired_reads_ch }

    //----------------------------------------------------------------
    // 5. MAPPING
    //----------------------------------------------------------------
    map_reads( repaired_reads_ch, genome_ch )

    //----------------------------------------------------------------
    // 6. MultiQC reports for each step
    //----------------------------------------------------------------
    fastp_trim_3.out
              .map{ sid, f1, f2, json, html -> [json, html] }
              .flatten()
              .set { fastp3_qc_ch }

    multiqc_fastp3(
        fastp3_qc_ch,
        "fastp_trim_3",
        params.multiqc_dir
    )

/*
    clumpify.out
              .map{ sid, r1, r2, stats -> stats }
              .set { clumpify_stats_ch }

    multiqc_clumpify(
        clumpify_stats_ch,
        "clumpify",
        params.multiqc_dir
    )
*/



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

}

//--------------------------------------------------------------------
// SUB‑WORKFLOW: fetch + index genome
//--------------------------------------------------------------------
workflow prepare_genome {
    take:
        accession
    main:
        fetch_genome(accession) | index_genome
    emit:
        index_genome.out
}

//--------------------------------------------------------------------
// MODULE IMPORTS
//--------------------------------------------------------------------
include { fastp_trim_3 }   from './modules/fastp_trim_3.nf'
include { clumpify }       from './modules/clumpify.nf'
include { fastp_trim_5 }   from './modules/fastp_trim_5.nf'
include { fastq_screen }   from './modules/fastq_screen.nf'
include { repair }         from './modules/repair.nf'
include { map_reads }      from './modules/map_reads.nf'
include { fetch_genome }   from './modules/fetch_genome.nf'
include { index_genome }   from './modules/index_genome.nf'

include { fastqc_raw }        from './modules/fastqc.nf'
include { multiqc as multiqc_raw }        from './modules/multiqc.nf'

// MultiQC (aliased per context to avoid duplicate process names)
include { multiqc as multiqc_fastp3 }      from './modules/multiqc.nf'
include { multiqc as multiqc_clumpify }    from './modules/multiqc.nf'
include { multiqc as multiqc_fastp5 }      from './modules/multiqc.nf'
include { multiqc as multiqc_fastqscreen } from './modules/multiqc.nf'
include { multiqc as multiqc_repair }      from './modules/multiqc.nf'