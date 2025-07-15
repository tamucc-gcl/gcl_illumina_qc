nextflow.enable.dsl = 2

/* ---------- modules ---------- (unchanged) */
include { fastqc_raw      } from './modules/fastqc.nf'
include { fastp_trim_3    } from './modules/fastp_trim_3.nf'
include { clumpify        } from './modules/clumpify.nf'
include { fastp_trim_5    } from './modules/fastp_trim_5.nf'
include { fastq_screen    } from './modules/fastq_screen.nf'
include { repair          } from './modules/repair.nf'

include { multiqc as mq_raw      } from './modules/multiqc.nf'
include { multiqc as mq_trim3    } from './modules/multiqc.nf'
include { multiqc as mq_clumpify } from './modules/multiqc.nf'
include { multiqc as mq_trim5    } from './modules/multiqc.nf'
include { multiqc as mq_screen   } from './modules/multiqc.nf'
include { multiqc as mq_repair   } from './modules/multiqc.nf'

include { fetch_genome } from './modules/fetch_genome.nf'
include { index_genome } from './modules/index_genome.nf'
include { map_reads    } from './modules/map_reads.nf'

/* ---------- parameters ---------- */
//params.accession = params.accession ?: 'GCF_PLACEHOLDER'   // replace or set via CLI
params.fastq_glob       = "data/fq_raw/*.{1,2}.fq.gz"    // paired‑end,  sampleID.1.fq.gz / .2.fq.gz
params.accession   = "GCA_042920385.1"              // NCBI assembly accession

/* ---------- channels ---------- */
Channel
    .fromFilePairs(params.fastq_glob, flat: true)
    .map { id, reads -> tuple(id, reads[0], reads[1]) }
    .set { reads_ch }

Channel
    .from( (params.accession instanceof String) ? params.accession.split(/[,\\s]+/) : params.accession )
    .set { accession_ch }

/* helper for MultiQC */
def collectForMQ = { ch, step ->
    ch.flatMap{ it[1..-1] }.collect()
      .map { files -> tuple(step, 'results/multiqc', files) }
}

/* ---------- workflow ---------- */
workflow {

    genome_fa_ch = accession_ch | fetch_genome | index_genome

    raw_out   = fastqc_raw(reads_ch)
    t3_out    = fastp_trim_3(reads_ch)
    clp_out   = clumpify(t3_out)
    t5_out    = fastp_trim_5(clp_out)
    scr_out   = fastq_screen(t5_out)
    rep_out   = repair(scr_out)

    /* reads + genome → map_reads (two‑channel call) */
    reads_list_ch = reads_ch.map { sid, r1, r2 -> tuple(sid, [r1, r2]) }
    map_reads_out = map_reads(reads_list_ch, genome_fa_ch)

    /* MultiQC tasks (unique aliases) */
    mq_raw   ( collectForMQ(raw_out , 'raw_fastqc') )
    mq_trim3 ( collectForMQ(t3_out  , 'fastp_trim_3') )
    mq_clumpify( collectForMQ(clp_out , 'clumpify') )
    mq_trim5 ( collectForMQ(t5_out  , 'fastp_trim_5') )
    mq_screen( collectForMQ(scr_out , 'fastq_screen') )
    mq_repair( collectForMQ(rep_out , 'repair') )
}
