// main.nf

nextflow.enable.dsl=2

params.reads        = "data/raw_fq/*.{1,2}.fq.gz"
params.accession    = "GCA_042920385.1"
params.outdir       = "results"
params.multiqc_dir  = "results/multiqc"

workflow {
    /* -----------------------------------------------------------------------------
     * INPUT
     * -------------------------------------------------------------------------- */
    Channel.fromFilePairs( params.reads, flat: true )
            .set { raw_reads_pairs }

    // ── MultiQC on raw data ──────────────────────────────────────────────────────
    raw_reads_pairs.collect()
                   .set { raw_reads_for_multiqc }
    multiqc_raw( raw_reads_for_multiqc, "raw", params.multiqc_dir )

    /* -----------------------------------------------------------------------------
     * GENOME PREPARATION (download  ➜  index)
     * -------------------------------------------------------------------------- */
    prepare_genome( params.accession )
        .set { genome_ch }

    /* -----------------------------------------------------------------------------
     * READ‑PREP CHAIN
     *  3′ trim  ➜  clumpify  ➜  5′ trim  ➜  fastq_screen  ➜  repair
     * -------------------------------------------------------------------------- */
    raw_reads_pairs                        |
        fastp_trim_3                      |
        clumpify                          |
        fastp_trim_5                      |
        fastq_screen                      |
        repair                            |
        set { repaired_reads_ch }         // cleaned, synchronised read pairs

    /* -----------------------------------------------------------------------------
     * MAPPING
     * -------------------------------------------------------------------------- */
    map_reads( repaired_reads_ch, genome_ch )

    /* -----------------------------------------------------------------------------
     * MULTIQC REPORTS FOR EACH STEP
     * -------------------------------------------------------------------------- */
    fastp_trim_3.out.collect()   | multiqc_fastp3(   it, "fastp_trim_3",  params.multiqc_dir )
    clumpify.out.collect()       | multiqc_clumpify( it, "clumpify",      params.multiqc_dir )
    fastp_trim_5.out.collect()   | multiqc_fastp5(   it, "fastp_trim_5",  params.multiqc_dir )
    fastq_screen.out.collect()   | multiqc_fastqscreen( it, "fastq_screen", params.multiqc_dir )
    repair.out.collect()         | multiqc_repair(   it, "repair",        params.multiqc_dir )
}

/* -------------------------------------------------------------------------------
 * SUB‑WORKFLOW: genome fetch + index
 * ----------------------------------------------------------------------------- */
workflow prepare_genome {
    take:
        accession
    main:
        fetch_genome( accession ) | index_genome
    emit:
        index_genome.out
}

/* -------------------------------------------------------------------------------
 * MODULE IMPORTS
 * ----------------------------------------------------------------------------- */
include { fastp_trim_3 }   from "./modules/fastp_trim_3.nf"
include { clumpify }       from "./modules/clumpify.nf"
include { fastp_trim_5 }   from "./modules/fastp_trim_5.nf"
include { fastq_screen }   from "./modules/fastq_screen.nf"
include { repair }         from "./modules/repair.nf"
include { map_reads }      from "./modules/map_reads.nf"
include { fetch_genome }   from "./modules/fetch_genome.nf"
include { index_genome }   from "./modules/index_genome.nf"

include { multiqc as multiqc_raw }         from "./modules/multiqc.nf"
include { multiqc as multiqc_fastp3 }      from "./modules/multiqc.nf"
include { multiqc as multiqc_clumpify }    from "./modules/multiqc.nf"
include { multiqc as multiqc_fastp5 }      from "./modules/multiqc.nf"
include { multiqc as multiqc_fastqscreen } from "./modules/multiqc.nf"
include { multiqc as multiqc_repair }      from "./modules/multiqc.nf"}
