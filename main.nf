/*
 * gcl_illumina_qc / main.nf
 * Full QC chain with consolidated MultiQC reporting and genome download,
 * indexing, and read mapping.
 *
 * QC stages consolidated:
 *   • raw FastQC
 *   • fastp_trim_3
 *   • clumpify
 *   • fastp_trim_5
 *   • fastq_screen
 *   • repair
 *
 * Mapping stage:
 *   • fetch_genome -> index_genome -> map_reads
 *
 * Edit the glob in `fromFilePairs` and the `params.accession` as needed.
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------
//  Module imports
// ---------------------------------------------------------------------
include { fastqc_raw      } from './modules/fastqc.nf'
include { fastp_trim_3    } from './modules/fastp_trim_3.nf'
include { clumpify        } from './modules/clumpify.nf'
include { fastp_trim_5    } from './modules/fastp_trim_5.nf'
include { fastq_screen    } from './modules/fastq_screen.nf'
include { repair          } from './modules/repair.nf'
include { multiqc         } from './modules/multiqc.nf'

include { fetch_genome    } from './modules/fetch_genome.nf'
include { index_genome    } from './modules/index_genome.nf'
include { map_reads       } from './modules/map_reads.nf'

// ---------------------------------------------------------------------
//  Parameters
// ---------------------------------------------------------------------
//params.accession = params.accession ?: 'GCF_PLACEHOLDER'   // replace or set via CLI
params.reads       = "data/fq_raw/*.{1,2}.fq.gz"    // paired‑end,  sampleID.1.fq.gz / .2.fq.gz
params.accession   = "GCA_042920385.1"              // NCBI assembly accession

// ---------------------------------------------------------------------
//  Inputs
// ---------------------------------------------------------------------
// Paired FASTQ reads
Channel
    .fromFilePairs(params.reads, flat: true)
    .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }
    .set { reads_ch }

// Single accession (or multiple, comma‑separated) for genome retrieval
Channel
    .from( params.accession instanceof String
           ? params.accession.split(/[,\s]+/)
           : params.accession )
    .set { accession_ch }

// ---------------------------------------------------------------------
//  Workflow definition
// ---------------------------------------------------------------------
workflow {

    // ---------- RAW FastQC ----------
    fastqc_raw_out = fastqc_raw(reads_ch)

    // ---------- fastp 3′ trim ----------
    trim3_out      = fastp_trim_3(reads_ch)

    // ---------- clumpify ----------
    clumpify_out   = clumpify(trim3_out)

    // ---------- fastp 5′ trim ----------
    trim5_out      = fastp_trim_5(clumpify_out)

    // ---------- FastQ‑Screen ----------
    fqscreen_out   = fastq_screen(trim5_out)

    // ---------- BBduk repair ----------
    repair_out     = repair(fqscreen_out)

    // ----------------------------------------------------------------
    //  Genome download and indexing
    // ----------------------------------------------------------------
    genome_fa_ch   = accession_ch
                     | fetch_genome
                     | index_genome

    // Broadcast genome so every sample gets it
    genome_fa_ch.broadcast().set { genome_bcast_ch }

    // ---------- Read mapping ----------
    reads_plus_genome = reads_ch.combine(genome_bcast_ch)
                          .map { read_tup, genome ->
                              tuple(read_tup[0], read_tup[1], read_tup[2], genome)
                          }
    map_reads_out   = map_reads(reads_plus_genome)

    // ----------------------------------------------------------------
    //  Consolidated MultiQC channels
    // ----------------------------------------------------------------
    def consolidate = { ch, step ->
        ch.flatMap { it[1..-1] }
          .collect()
          .map { files -> tuple(step, 'results/multiqc', files) }
    }

    multiqc_raw_ch       = consolidate(fastqc_raw_out      , 'raw_fastqc')
    multiqc_trim3_ch     = consolidate(trim3_out           , 'fastp_trim_3')
    multiqc_clumpify_ch  = consolidate(clumpify_out        , 'clumpify')
    multiqc_trim5_ch     = consolidate(trim5_out           , 'fastp_trim_5')
    multiqc_screen_ch    = consolidate(fqscreen_out        , 'fastq_screen')
    multiqc_repair_ch    = consolidate(repair_out          , 'repair')

    // Merge all QC channels
    Channel.merge(
        multiqc_raw_ch,
        multiqc_trim3_ch,
        multiqc_clumpify_ch,
        multiqc_trim5_ch,
        multiqc_screen_ch,
        multiqc_repair_ch
    )
    | multiqc
}
