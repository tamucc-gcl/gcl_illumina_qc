nextflow.enable.dsl = 2

/*
 * Full pipeline compatible with Nextflow 25.04.2
 */

/* QC modules */
include { fastqc_raw      } from './modules/fastqc.nf'
include { fastp_trim_3    } from './modules/fastp_trim_3.nf'
include { clumpify        } from './modules/clumpify.nf'
include { fastp_trim_5    } from './modules/fastp_trim_5.nf'
include { fastq_screen    } from './modules/fastq_screen.nf'
include { repair          } from './modules/repair.nf'

/* MultiQC aliases */
include { multiqc as multiqc_raw      } from './modules/multiqc.nf'
include { multiqc as multiqc_trim3    } from './modules/multiqc.nf'
include { multiqc as multiqc_clumpify } from './modules/multiqc.nf'
include { multiqc as multiqc_trim5    } from './modules/multiqc.nf'
include { multiqc as multiqc_screen   } from './modules/multiqc.nf'
include { multiqc as multiqc_repair   } from './modules/multiqc.nf'

/* Genome modules */
include { fetch_genome } from './modules/fetch_genome.nf'
include { index_genome } from './modules/index_genome.nf'
include { map_reads    } from './modules/map_reads.nf'

/* Parameters */
params.fastq_glob = params.fastq_glob ?: 'data/fq_raw/*.{1,2}.fq.gz'
//params.accession  = params.accession  ?: 'GCF_PLACEHOLDER'
params.accession   = "GCA_042920385.1"              // NCBI assembly accession

/* Input channels */
Channel
    .fromFilePairs(params.fastq_glob, flat: true)
    .map { sid, reads -> tuple(sid, reads[0], reads[1]) }
    .set { reads_ch }

Channel
    .from( (params.accession instanceof String) ? params.accession.split(/[,\s]+/) : params.accession )
    .set { accession_ch }

/* Helper for MultiQC */
def consolidate = { ch, step ->
    ch.flatMap { it[1..-1] }
      .collect()
      .map { files -> tuple(step, 'results/multiqc', files) }
}

workflow {

    /* Genome download + index */
    genome_fa_ch = accession_ch | fetch_genome | index_genome

    /* QC chain */
    fastqc_raw_out = fastqc_raw(reads_ch)
    trim3_out      = fastp_trim_3(reads_ch)
    clumpify_out   = clumpify(trim3_out)
    trim5_out      = fastp_trim_5(clumpify_out)
    fqscreen_out   = fastq_screen(trim5_out)
    repair_out     = repair(fqscreen_out)

    /* Prepare reads tuple list for map_reads */
    reads_list_ch = reads_ch.map { sid, r1, r2 -> tuple(sid, [r1, r2]) }

    /* Map reads */
    map_reads_out = map_reads(reads_list_ch, genome_fa_ch)

    /* MultiQC reports */
    multiqc_raw      ( consolidate(fastqc_raw_out , 'raw_fastqc') )
    multiqc_trim3    ( consolidate(trim3_out      , 'fastp_trim_3') )
    multiqc_clumpify ( consolidate(clumpify_out   , 'clumpify') )
    multiqc_trim5    ( consolidate(trim5_out      , 'fastp_trim_5') )
    multiqc_screen   ( consolidate(fqscreen_out   , 'fastq_screen') )
    multiqc_repair   ( consolidate(repair_out     , 'repair') )
}
