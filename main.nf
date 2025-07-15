// main.nf – sequential read‑prep + per‑step MultiQC

nextflow.enable.dsl = 2

params.reads       = "data/raw_fq/*.{1,2}.fq.gz"   // paired‑end R1/R2 files
params.accession   = "GCA_042920385.1"              // NCBI assembly to fetch
params.outdir      = "results"
params.multiqc_dir = "${params.outdir}/multiqc"

/* ========================================================================== */
workflow {
    /* ----------------------------------------------------------------------
     * INPUT (raw pairs + early MultiQC)
     * ------------------------------------------------------------------- */
    Channel.fromFilePairs( params.reads, flat: true )
            .set { raw_reads_pairs }

    // ── MultiQC on raw data ----------------------------------------------
    raw_reads_pairs.collect().set { raw_reads_for_multiqc }
    multiqc_raw( raw_reads_for_multiqc, "raw", params.multiqc_dir )

    /* ----------------------------------------------------------------------
     * GENOME PREPARATION (download ➜ index)
     * ------------------------------------------------------------------- */
    prepare_genome( params.accession ).set { genome_ch }

    /* ----------------------------------------------------------------------
     * READ‑PREP CHAIN
     *   3′‑trim ➜ clumpify ➜ 5′‑trim ➜ fastq_screen ➜ repair
     * ------------------------------------------------------------------- */
    raw_reads_pairs \
        | fastp_trim_3 \
        | clumpify     \
        | fastp_trim_5 \
        | fastq_screen \
        | repair       \
        | set { repaired_reads_ch }

    /* ----------------------------------------------------------------------
     * MAPPING
     * ------------------------------------------------------------------- */
    map_reads( repaired_reads_ch, genome_ch )

    /* ----------------------------------------------------------------------
     * MULTIQC REPORTS FOR EACH PREP STEP
     * (pass the collected output channel as the first argument)
     * ------------------------------------------------------------------- */
    multiqc_fastp3(      fastp_trim_3.out.collect(),  "fastp_trim_3",   params.multiqc_dir )
    multiqc_clumpify(    clumpify.out.collect(),      "clumpify",       params.multiqc_dir )
    multiqc_fastp5(      fastp_trim_5.out.collect(),  "fastp_trim_5",   params.multiqc_dir )
    multiqc_fastqscreen( fastq_screen.out.collect(),  "fastq_screen",   params.multiqc_dir )
    multiqc_repair(      repair.out.collect(),        "repair",         params.multiqc_dir )
}

/* ========================================================================== */
workflow prepare_genome {
    // sub‑workflow: fetch genome FASTA → build BWA‑MEM2 index
    take:
        accession
    main:
        fetch_genome(accession) | index_genome
    emit:
        index_genome.out
}

/* ========================================================================== */
// MODULE IMPORTS ------------------------------------------------------------
include { fastp_trim_3 }   from './modules/fastp_trim_3.nf'
include { clumpify }       from './modules/clumpify.nf'
include { fastp_trim_5 }   from './modules/fastp_trim_5.nf'
include { fastq_screen }   from './modules/fastq_screen.nf'
include { repair }         from './modules/repair.nf'
include { map_reads }      from './modules/map_reads.nf'
include { fetch_genome }   from './modules/fetch_genome.nf'
include { index_genome }   from './modules/index_genome.nf'

// MultiQC (aliased per context to avoid duplicate process names)
include { multiqc as multiqc_raw }         from './modules/multiqc.nf'
include { multiqc as multiqc_fastp3 }      from './modules/multiqc.nf'
include { multiqc as multiqc_clumpify }    from './modules/multiqc.nf'
include { multiqc as multiqc_fastp5 }      from './modules/multiqc.nf'
include { multiqc as multiqc_fastqscreen } from './modules/multiqc.nf'
include { multiqc as multiqc_repair }      from './modules/multiqc.nf'
