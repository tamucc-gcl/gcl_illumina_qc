/*
 * gcl_illumina_qc / main.nf
 * Complete QC chain with consolidated MultiQC reporting for:
 *   • raw FastQC
 *   • fastp_trim_3
 *   • clumpify
 *   • fastp_trim_5
 *   • fastq_screen
 *   • repair
 *
 * Assumes each module emits a tuple whose first element is `val(sample_id)`
 * followed by one or more `path` items (reads, logs, html, json, etc.).
 * MultiQC ignores file types it doesn't recognise, so we simply feed all
 * non‑sample‑id paths to the report step.
 *
 * Update the glob in `fromFilePairs` to match your raw fastq layout.
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

// ---------------------------------------------------------------------
//  Inputs
// ---------------------------------------------------------------------
Channel
    .fromFilePairs('fq_raw/*_{1,2}.fq.gz', flat: true)
    .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }
    .set { reads_ch }

// ---------------------------------------------------------------------
//  Workflow definition
// ---------------------------------------------------------------------
workflow {

    // ----------------  Step 0: raw FastQC  ---------------------------
    fastqc_raw_out = fastqc_raw(reads_ch)

    // ----------------  Step 1: fastp 3′ trim  ------------------------
    trim3_out      = fastp_trim_3(reads_ch)

    // ----------------  Step 2: clumpify  -----------------------------
    clumpify_out   = clumpify(trim3_out)

    // ----------------  Step 3: fastp 5′ trim  ------------------------
    trim5_out      = fastp_trim_5(clumpify_out)

    // ----------------  Step 4: FastQ‑Screen  -------------------------
    fqscreen_out   = fastq_screen(trim5_out)

    // ----------------  Step 5: BBduk repair  -------------------------
    repair_out     = repair(fqscreen_out)

    // ----------------------------------------------------------------
    //  Consolidated MultiQC channels
    // ----------------------------------------------------------------
    multiqc_raw_ch       = fastqc_raw_out                                        \
                            .flatMap{ it[1..-1] }                                \
                            .collect()                                           \
                            .map{ files -> tuple('raw_fastqc', 'results/multiqc', files) }

    multiqc_trim3_ch     = trim3_out                                             \
                            .flatMap{ it[1..-1] }                                \
                            .collect()                                           \
                            .map{ files -> tuple('fastp_trim_3', 'results/multiqc', files) }

    multiqc_clumpify_ch  = clumpify_out                                          \
                            .flatMap{ it[1..-1] }                                \
                            .collect()                                           \
                            .map{ files -> tuple('clumpify', 'results/multiqc', files) }

    multiqc_trim5_ch     = trim5_out                                             \
                            .flatMap{ it[1..-1] }                                \
                            .collect()                                           \
                            .map{ files -> tuple('fastp_trim_5', 'results/multiqc', files) }

    multiqc_screen_ch    = fqscreen_out                                          \
                            .flatMap{ it[1..-1] }                                \
                            .collect()                                           \
                            .map{ files -> tuple('fastq_screen', 'results/multiqc', files) }

    multiqc_repair_ch    = repair_out                                            \
                            .flatMap{ it[1..-1] }                                \
                            .collect()                                           \
                            .map{ files -> tuple('repair', 'results/multiqc', files) }

    // Merge all QC channels and spawn MultiQC once per `step`
    [ multiqc_raw_ch,
      multiqc_trim3_ch,
      multiqc_clumpify_ch,
      multiqc_trim5_ch,
      multiqc_screen_ch,
      multiqc_repair_ch ]                                                         \
      .mix()                                                                      \
      | multiqc
}
