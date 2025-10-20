// main.nf – GCL Illumina QC pipeline
// July 2025 – consolidated MultiQC per step
// Oct 2025 - implement local genome option

nextflow.enable.dsl = 2

//--------------------------------------------------------------------
// USER PARAMETERS
//--------------------------------------------------------------------
params.reads       = "data/fq_raw/*.{1,2}.fq.gz"    // paired‑end,  sampleID.1.fq.gz / .2.fq.gz
params.accession   = null                            // NCBI assembly accession (optional)
params.genome      = null                            // Path to local genome file (optional)
params.decontam_conffile    = "/work/birdlab/fastq_screen_databases/example_fastq-screen.conf"  // FastQ Screen config file
params.sequencing_type = "whole_genome"  // Options: "ddrad" or "whole_genome"
params.outdir      = "results"

//--------------------------------------------------------------------
// WORKFLOW DEFINITION
//--------------------------------------------------------------------
workflow {
    // Validate genome input parameters
    if (!params.genome && !params.accession) {
        error "Error: Either --genome (local file) or --accession (NCBI) must be specified"
    }
    if (params.genome && params.accession) {
        log.warn "Warning: Both --genome and --accession specified. Using local genome: ${params.genome}"
    }

    //----------------------------------------------------------------
    // 1. RAW READ INPUT
    //----------------------------------------------------------------
    Channel
        .fromFilePairs( params.reads, flat: true )
        .set { raw_reads_pairs }

    //----------------------------------------------------------------
    // 2. FASTQC ON RAW READS  ➜  MULTIQC (raw_fastqc)
    //----------------------------------------------------------------
    fastqc_raw( raw_reads_pairs )

    // MultiQC for raw FastQC
    multiqc_raw_out = multiqc_raw( 
        fastqc_raw.out
            .map{ sid, html1, zip1, html2, zip2 -> [html1, zip1, html2, zip2] }
            .flatten()
            .collect(),
        Channel.value('raw_fastqc')
    )

    //----------------------------------------------------------------
    // 3. GENOME PREP  (download ➜ index OR local ➜ index)
    //----------------------------------------------------------------
    if (params.genome) {
        // Use local genome file
        genome_file = file(params.genome, checkIfExists: true)
        prepare_genome_local( Channel.value(genome_file) )
        genome_indexed = prepare_genome_local.out.genome
    } else {
        // Download genome from NCBI
        prepare_genome( Channel.value( params.accession ) )
        genome_indexed = prepare_genome.out.genome
    }

    //----------------------------------------------------------------
    // 4. QC PIPELINE STEPS
    //----------------------------------------------------------------
    
    // Step 1: 3' trimming with fastp
    fastp_trim_3( raw_reads_pairs )
    
    // FastQC after 3' trimming
    fastqc_trim3( 
        fastp_trim_3.out.map{ sid, r1, r2, json, html -> tuple(sid, r1, r2) }
    )
    
    // MultiQC for 3' trimming (fastp + FastQC)
    multiqc_trim3_out = multiqc_trim3(
        fastp_trim_3.out
            .map{ sid, r1, r2, json, html -> [json, html] }
            .flatten()
            .mix( 
                fastqc_trim3.out
                    .map{ sid, html1, zip1, html2, zip2 -> [html1, zip1, html2, zip2] }
                    .flatten()
            )
            .collect(),
        Channel.value('fastp_trim_3')
    )
    
    // Step 2: Clumpify
    clumpify( fastp_trim_3.out,
              Channel.value(params.sequencing_type))
    
    // FastQC after clumpify
    fastqc_clumpify( 
        clumpify.out.map{ sid, r1, r2, stats -> tuple(sid, r1, r2) }
    )
    
    // MultiQC for clumpify (stats + FastQC)
    multiqc_clumpify_out = multiqc_clumpify(
        clumpify.out
            .map{ sid, r1, r2, stats -> stats }
            .mix( 
                fastqc_clumpify.out
                    .map{ sid, html1, zip1, html2, zip2 -> [html1, zip1, html2, zip2] }
                    .flatten()
            )
            .collect(),
        Channel.value('clumpify')
    )
    
    // Step 3: 5' trimming with fastp
    fastp_trim_5( clumpify.out.map{ sid, r1, r2, stats -> tuple(sid, r1, r2) } )
    
    // FastQC after 5' trimming
    fastqc_trim5( 
        fastp_trim_5.out.map{ sid, r1, r2, json, html -> tuple(sid, r1, r2) }
    )
    
    // MultiQC for 5' trimming (fastp + FastQC)
    multiqc_trim5_out = multiqc_trim5(
        fastp_trim_5.out
            .map{ sid, r1, r2, json, html -> [json, html] }
            .flatten()
            .mix( 
                fastqc_trim5.out
                    .map{ sid, html1, zip1, html2, zip2 -> [html1, zip1, html2, zip2] }
                    .flatten()
            )
            .collect(),
        Channel.value('fastp_trim_5')
    )
    
    // Step 4: FastQ Screen
    fastp_trim_5.out
        .flatMap{ sid, r1, r2, json, html -> [
            [sid, r1, "1"],
            [sid, r2, "2"]
        ]}
        .set { individual_reads }
    
    // Create a value channel for the config file (broadcast to all processes)
    config_ch = Channel.value(file(params.decontam_conffile))
    
    fastq_screen( individual_reads, config_ch )

    // Group fastq_screen results back together for repair
    fastq_screen.out
        .groupTuple(by: 0)
        .map{ sid, reads, reports, read_nums -> 
            // Sort by read number to ensure R1, R2 order
            def sorted = [reads, reports, read_nums].transpose().sort{ it[2] }
            [sid, sorted[0][0], sorted[1][0], sorted[0][1], sorted[1][1]]
        }
        .set { screen_paired }
    
    // FastQC after fastq_screen
    fastqc_screen( 
        screen_paired.map{ sid, r1, r2, txt1, txt2 -> tuple(sid, r1, r2) }
    )
    
    // MultiQC for fastq_screen (screen reports + FastQC)
    multiqc_screen_out = multiqc_screen(
        screen_paired
            .map{ sid, r1, r2, txt1, txt2 -> [txt1, txt2] }
            .flatten()
            .mix( 
                fastqc_screen.out
                    .map{ sid, html1, zip1, html2, zip2 -> [html1, zip1, html2, zip2] }
                    .flatten()
            )
            .collect(),
        Channel.value('fastq_screen')
    )
    
    // Step 5: Repair
    repair( screen_paired.map{ sid, r1, r2, txt1, txt2 -> tuple(sid, r1, r2) } )
    
    
    // FastQC after repair
    fastqc_repair( repair.out )
    
    // MultiQC for repair (FastQC only)
    multiqc_repair_out = multiqc_repair(
        fastqc_repair.out
            .map{ sid, html1, zip1, html2, zip2 -> [html1, zip1, html2, zip2] }
            .flatten()
            .collect(),
        Channel.value('repair')
    )
    
    // Step 6: Map reads to genome (use the indexed genome from either path)
    map_reads( repair.out, genome_indexed )
    
    // Step 7: Generate BAM statistics
    samtools_stats( map_reads.out )
    
    // MultiQC for mapping (samtools stats + flagstats)
    multiqc_mapping_out = multiqc_mapping(
        samtools_stats.out
            .map{ sid, stats, flagstats -> [stats, flagstats] }
            .flatten()
            .collect(),
        Channel.value('mapping')
    )

    // Create mapping summary instead of problematic MultiQC final step
    samtools_summary(
        samtools_stats.out
            .map{ sid, stats, flagstats -> stats }
            .collect()
    )

    // Collect all the general stats files using indexed outputs
    // multiqc_*_out[0] = HTML file, multiqc_*_out[1] = stats file
    all_stats = Channel.empty()
        .mix(
            multiqc_raw_out[1],        // Get the stats file (index 1)
            multiqc_trim3_out[1],
            multiqc_clumpify_out[1],
            multiqc_trim5_out[1],
            multiqc_screen_out[1],
            multiqc_repair_out[1],
            multiqc_mapping_out[1],    // Add mapping MultiQC stats
            samtools_summary.out       // Add the mapping summary file
        )
        .collect()
    
    // Run R analysis on collected stats
    analyze_read_stats(all_stats)
    
}

//--------------------------------------------------------------------
// SUB‑WORKFLOW: fetch + index genome (from NCBI)
//--------------------------------------------------------------------
workflow prepare_genome {
    take:
        accession
    
    main:
        fetch_genome(accession)
        index_genome(fetch_genome.out)  
        
    emit:
        genome = index_genome.out[0]  // Output: tuple path(genome), val(genome_path), path(index_files)
        index_files = index_genome.out.index_files  // Named output: index files
}

//--------------------------------------------------------------------
// SUB‑WORKFLOW: prepare local genome (copy + index)
//--------------------------------------------------------------------
workflow prepare_genome_local {
    take:
        genome_file
    
    main:
        // Create a process to stage the local genome file
        stage_local_genome(genome_file)
        index_genome(stage_local_genome.out)
        
    emit:
        genome = index_genome.out[0]  // Output: tuple path(genome), val(genome_path), path(index_files)
        index_files = index_genome.out.index_files  // Named output: index files
}

//--------------------------------------------------------------------
// MODULE IMPORTS
//--------------------------------------------------------------------
include { fastp_trim_3 }      from './modules/fastp_trim_3.nf'
include { clumpify }          from './modules/clumpify.nf'
include { fastp_trim_5 }      from './modules/fastp_trim_5.nf'
include { fastq_screen }      from './modules/fastq_screen.nf'
include { repair }            from './modules/repair.nf'
include { map_reads }         from './modules/map_reads.nf'
include { samtools_stats }    from './modules/samtools_stats.nf'
include { samtools_summary }  from './modules/samtools_summary.nf'
include { fetch_genome }      from './modules/fetch_genome.nf'
include { index_genome }      from './modules/index_genome.nf'
include { stage_local_genome } from './modules/stage_local_genome.nf'
include { analyze_read_stats } from './modules/analyze_read_stats.nf'

include { fastqc_raw }        from './modules/fastqc.nf'
include { fastqc_generic as fastqc_trim3 }    from './modules/fastqc.nf'
include { fastqc_generic as fastqc_clumpify } from './modules/fastqc.nf'
include { fastqc_generic as fastqc_trim5 }    from './modules/fastqc.nf'
include { fastqc_generic as fastqc_screen }   from './modules/fastqc.nf'
include { fastqc_generic as fastqc_repair }   from './modules/fastqc.nf'

// MultiQC processes with aliases
include { multiqc as multiqc_raw }      from './modules/multiqc.nf'
include { multiqc as multiqc_trim3 }    from './modules/multiqc.nf'
include { multiqc as multiqc_clumpify } from './modules/multiqc.nf'
include { multiqc as multiqc_trim5 }    from './modules/multiqc.nf'
include { multiqc as multiqc_screen }   from './modules/multiqc.nf'
include { multiqc as multiqc_repair }   from './modules/multiqc.nf'
include { multiqc as multiqc_mapping }   from './modules/multiqc.nf'