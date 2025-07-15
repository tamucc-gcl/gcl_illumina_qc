/*
 * gcl_illumina_qc / main.nf
 * ------------------------------------------------------------------
 * Compatible with Nextflow 25.04.2  (no .mix(), .merge(), broadcast())
 * ------------------------------------------------------------------
 * Pipeline outline
 *   0.  Input paired FASTQ files
 *   1.  raw           : FastQC on raw reads
 *   2.  fastp_trim_3  : 3′ trimming with fastp
 *   3.  clumpify      : optical duplicate removal
 *   4.  fastp_trim_5  : 5′ trimming with fastp
 *   5.  fastq_screen  : screen against contaminant genomes
 *   6.  repair        : BBduk repair to synchronise read pairs
 *   7.  Mapping       : bwa‑mem2 map to reference genome
 *   8.  MultiQC       : one consolidated report for each QC stage
 *   9.  Outputs       : BAM + MultiQC HTMLs under results/
 *
 * Each MultiQC task is imported under its own alias to avoid
 * "process already used" errors on older Nextflow versions.
 */

nextflow.enable.dsl = 2

/* ----------------  QC module imports  ---------------- */
include { fastqc_raw      } from './modules/fastqc.nf'
include { fastp_trim_3    } from './modules/fastp_trim_3.nf'
include { clumpify        } from './modules/clumpify.nf'
include { fastp_trim_5    } from './modules/fastp_trim_5.nf'
include { fastq_screen    } from './modules/fastq_screen.nf'
include { repair          } from './modules/repair.nf'

/*  MultiQC module imported six times with unique aliases  */
include { multiqc as multiqc_raw      } from './modules/multiqc.nf'
include { multiqc as multiqc_trim3    } from './modules/multiqc.nf'
include { multiqc as multiqc_clumpify } from './modules/multiqc.nf'
include { multiqc as multiqc_trim5    } from './modules/multiqc.nf'
include { multiqc as multiqc_screen   } from './modules/multiqc.nf'
include { multiqc as multiqc_repair   } from './modules/multiqc.nf'

/* ----------------  Genome module imports  ---------------- */
include { fetch_genome } from './modules/fetch_genome.nf'
include { index_genome } from './modules/index_genome.nf'
include { map_reads    } from './modules/map_reads.nf'

/* ----------------  Parameters  ---------------- */
// Genome assembly accession(s) (comma or space‑separated) e.g. --accession GCF_000001405.40
//params.accession = params.accession ?: 'GCF_PLACEHOLDER'   // replace or set via CLI
params.accession   = "GCA_042920385.1"              // NCBI assembly accession
// Directory pattern for raw FASTQ pairs (change as needed)
params.fastq_glob = params.fastq_glob ?: 'data/fq_raw/*.{1,2}.fq.gz'    // paired‑end,  sampleID.1.fq.gz / .2.fq.gz



/* ----------------  Input channels  ---------------- */
// Paired FASTQ channel: (sample_id, R1, R2)
Channel
    .fromFilePairs(params.fastq_glob, flat: true)
    .map { sid, reads -> tuple(sid, reads[0], reads[1]) }
    .set { reads_ch }

// Genome accessions channel
Channel
    .from( (params.accession instanceof String)
           ? params.accession.split(/[,\s]+/)
           : params.accession )
    .set { accession_ch }

/* ----------------  Helper: consolidate for MultiQC  ---------------- */
/*
 * Takes a QC output channel whose first element is sample_id and remaining
 * elements are files. Flattens to files only, collects across samples,
 * and packages into a tuple (step, outdir, files) expected by multiqc.nf
 */
def consolidate = { ch, step ->
    ch.flatMap { it[1..-1] }
      .collect()
      .map { files -> tuple(step, 'results/multiqc', files) }
}

/* ----------------  Workflow  ---------------- */
workflow {

    /* --- Genome download + indexing --- */
    genome_fa_ch = accession_ch
                   | fetch_genome
                   | index_genome

    /* --- QC chain --- */
    fastqc_raw_out = fastqc_raw(reads_ch)

    trim3_out      = fastp_trim_3(reads_ch)

    clumpify_out   = clumpify(trim3_out)

    trim5_out      = fastp_trim_5(clumpify_out)

    fqscreen_out   = fastq_screen(trim5_out)

    repair_out     = repair(fqscreen_out)

    /* --- Mapping --- */
    // Combine each read tuple with the (single) genome fasta
    map_in_ch = reads_ch.combine(genome_fa_ch)
                        .map { r, genome ->
                            tuple(r[0], r[1], r[2], genome)
                        }

    map_reads_out = map_reads(map_in_ch)

    /* --- MultiQC reports (one per stage) --- */
    multiqc_raw      ( consolidate(fastqc_raw_out , 'raw_fastqc') )
    multiqc_trim3    ( consolidate(trim3_out      , 'fastp_trim_3') )
    multiqc_clumpify ( consolidate(clumpify_out   , 'clumpify') )
    multiqc_trim5    ( consolidate(trim5_out      , 'fastp_trim_5') )
    multiqc_screen   ( consolidate(fqscreen_out   , 'fastq_screen') )
    multiqc_repair   ( consolidate(repair_out     , 'repair') )
}
