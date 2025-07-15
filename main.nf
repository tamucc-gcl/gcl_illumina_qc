/*
 * gcl_illumina_qc / main.nf
 * Full pipeline with genome download/indexing, mapping, QC stages, and consolidated MultiQC reports.
 */

nextflow.enable.dsl = 2

// ------------------  Module imports  ------------------
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

// ------------------  Parameters  ------------------
//params.accession = params.accession ?: 'GCF_PLACEHOLDER'   // replace or set via CLI
params.reads       = "data/fq_raw/*.{1,2}.fq.gz"    // paired‑end,  sampleID.1.fq.gz / .2.fq.gz
params.accession   = "GCA_042920385.1"              // NCBI assembly accession

// ------------------  Input channels  ------------------
// Paired FASTQ reads
Channel
    .fromFilePairs(params.reads, flat: true)
    .map { id, reads -> tuple(id, reads[0], reads[1]) }
    .set { reads_ch }

// Genome accessions
Channel
    .from( (params.accession instanceof String)
           ? params.accession.split(/[,\s]+/)
           : params.accession )
    .set { accession_ch }

// ------------------  Workflow  ------------------
workflow {

    // Genome download → index
    genome_fa_ch = accession_ch | fetch_genome | index_genome

    // QC processing chain
    fastqc_raw_out = fastqc_raw(reads_ch)
    trim3_out      = fastp_trim_3(reads_ch)
    clumpify_out   = clumpify(trim3_out)
    trim5_out      = fastp_trim_5(clumpify_out)
    fqscreen_out   = fastq_screen(trim5_out)
    repair_out     = repair(fqscreen_out)

    // Read mapping (each sample paired with the single genome)
    map_reads_out  = map_reads(reads_ch, genome_fa_ch)

    // ---------- Consolidated MultiQC ----------
    def consolidate = { channel, step ->
        channel
            .flatMap { it[1..-1] }     // keep only file paths
            .collect()
            .map { files -> tuple(step, 'results/multiqc', files) }
    }

    Channel.merge(
        consolidate(fastqc_raw_out , 'raw_fastqc'),
        consolidate(trim3_out      , 'fastp_trim_3'),
        consolidate(clumpify_out   , 'clumpify'),
        consolidate(trim5_out      , 'fastp_trim_5'),
        consolidate(fqscreen_out   , 'fastq_screen'),
        consolidate(repair_out     , 'repair')
    )
    | multiqc
}
