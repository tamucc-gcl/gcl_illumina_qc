/*
 * gcl_illumina_qc / main.nf
 * Pipeline compatible with Nextflow 25.04.2
 * - Genome download & index
 * - QC chain
 * - Mapping
 * - Consolidated MultiQC reports (one task per stage using unique aliases)
 */

nextflow.enable.dsl = 2

/* ----------------  Module imports  ---------------- */
include { fastqc_raw      } from './modules/fastqc.nf'
include { fastp_trim_3    } from './modules/fastp_trim_3.nf'
include { clumpify        } from './modules/clumpify.nf'
include { fastp_trim_5    } from './modules/fastp_trim_5.nf'
include { fastq_screen    } from './modules/fastq_screen.nf'
include { repair          } from './modules/repair.nf'

// MultiQC module imported with unique aliases
include { multiqc as multiqc_raw      } from './modules/multiqc.nf'
include { multiqc as multiqc_trim3    } from './modules/multiqc.nf'
include { multiqc as multiqc_clumpify } from './modules/multiqc.nf'
include { multiqc as multiqc_trim5    } from './modules/multiqc.nf'
include { multiqc as multiqc_screen   } from './modules/multiqc.nf'
include { multiqc as multiqc_repair   } from './modules/multiqc.nf'

include { fetch_genome    } from './modules/fetch_genome.nf'
include { index_genome    } from './modules/index_genome.nf'
include { map_reads       } from './modules/map_reads.nf'

/* ----------------  Parameters  ---------------- */
//params.accession = params.accession ?: 'GCF_PLACEHOLDER'   // replace or set via CLI
params.reads       = "data/fq_raw/*.{1,2}.fq.gz"    // paired‑end,  sampleID.1.fq.gz / .2.fq.gz
params.accession   = "GCA_042920385.1"              // NCBI assembly accession

/* ----------------  Input channels  ---------------- */
Channel
    .fromFilePairs(params.reads, flat: true)
    .map { sid, reads -> tuple(sid, reads[0], reads[1]) }
    .set { reads_ch }

Channel
    .from( (params.accession instanceof String)
           ? params.accession.split(/[,\s]+/)
           : params.accession )
    .set { accession_ch }

/* ----------------  Helper to pack MultiQC inputs  ---------------- */
def consolidate = { ch, step ->
    ch.flatMap { it[1..-1] }              // drop sample_id
      .collect()
      .map { files -> tuple(step, 'results/multiqc', files) }
}

/* ----------------  Workflow  ---------------- */
workflow {

    // Genome retrieval & indexing
    genome_fa_ch = accession_ch | fetch_genome | index_genome

    // QC processes
    fastqc_raw_out = fastqc_raw(reads_ch)
    trim3_out      = fastp_trim_3(reads_ch)
    clumpify_out   = clumpify(trim3_out)
    trim5_out      = fastp_trim_5(clumpify_out)
    fqscreen_out   = fastq_screen(trim5_out)
    repair_out     = repair(fqscreen_out)

    // Mapping
    map_reads_out  = map_reads(reads_ch, genome_fa_ch)

    // Launch MultiQC once per stage
    multiqc_raw   ( consolidate(fastqc_raw_out , 'raw_fastqc') )
    multiqc_trim3 ( consolidate(trim3_out      , 'fastp_trim_3') )
    multiqc_clumpify( consolidate(clumpify_out , 'clumpify') )
    multiqc_trim5 ( consolidate(trim5_out      , 'fastp_trim_5') )
    multiqc_screen( consolidate(fqscreen_out   , 'fastq_screen') )
    multiqc_repair( consolidate(repair_out     , 'repair') )
}
