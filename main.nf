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
    
//----------------------------------------------------------------
// 2. QC on raw reads (FastQC + MultiQC)
//----------------------------------------------------------------
fastqc_raw_out = fastqc_raw( raw_reads_pairs )

multiqc_raw(
    fastqc_raw_out.out
                  .map{ sid, files -> files }   // drop the ID
                  .flatten()
                  .collect(),
    "raw_fastqc",
    params.multiqc_dir
)



    //----------------------------------------------------------------
    // 3. GENOME PREPARATION (download ➜ index)
    //----------------------------------------------------------------
    prepare_genome( params.accession )
        .set { genome_ch }

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
    multiqc_fastp3(
    fastp_trim_3.out
               .map{ sid, reads -> reads }   // drop the ID
               .flatten()                    // turn [[R1,R2],[R1,R2]…] → [R1,R2,R1,R2…]
               .collect(),
    "fastp_trim_3",
    params.multiqc_dir
    )

    multiqc_clumpify(
    clumpify.out
               .map{ sid, reads -> reads }   // drop the ID
               .flatten()                    // turn [[R1,R2],[R1,R2]…] → [R1,R2,R1,R2…]
               .collect(),
    "clumpify",
    params.multiqc_dir
    )

    multiqc_fastp5(
    fastp_trim_5.out
               .map{ sid, reads -> reads }   // drop the ID
               .flatten()                    // turn [[R1,R2],[R1,R2]…] → [R1,R2,R1,R2…]
               .collect(),
    "fastp_trim_5",
    params.multiqc_dir
    )

    multiqc_fastqscreen(
    fastq_screen.out
               .map{ sid, reads -> reads }   // drop the ID
               .flatten()                    // turn [[R1,R2],[R1,R2]…] → [R1,R2,R1,R2…]
               .collect(),
    "fastq_screen",
    params.multiqc_dir
    )

    multiqc_repair(
    repair.out
               .map{ sid, reads -> reads }   // drop the ID
               .flatten()                    // turn [[R1,R2],[R1,R2]…] → [R1,R2,R1,R2…]
               .collect(),
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