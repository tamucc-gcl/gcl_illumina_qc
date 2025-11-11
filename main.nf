// main.nf – GCL Illumina QC pipeline
// July 2025 – consolidated MultiQC per step
// Oct 2025 - implement local genome option
// Nov 2025 - add de novo assembly option

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

// Species ID parameters
params.run_species_id = true  // Enable/disable species identification
params.mito_reference = "databases/mito_gene_refs.fasta" // Mitochondrial gene reference provided with gcl_illumina_qc compiled from https://github.com/cmayer/MitoGeneExtractor
params.genetic_code = 2  //Mitochondrial genetic code: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
params.blast_db = null  // Path to local BLAST database (null = use NCBI nt)

// Assembly parameters
params.cutoff1 = 4          // Minimum reads per individual
params.cutoff2 = 4          // Minimum number of individuals  
params.cluster_similarity = 0.8
params.div_f = 0.5
params.div_K = 10
params.merge_r = 2
params.final_similarity = 0.9

//--------------------------------------------------------------------
// WORKFLOW DEFINITION
//--------------------------------------------------------------------

workflow {
    // Validate genome input parameters
    if (!params.genome && !params.accession && params.assembly_mode == "none") {
        log.info "No reference genome provided - will output cleaned FASTQ files only"
    } else if (!params.genome && !params.accession && params.assembly_mode == "denovo") {
        log.info "No reference genome provided - will perform de novo assembly"
    } else if (params.genome && params.accession) {
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
    // 3. QC PIPELINE STEPS - through repair
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
    
    //----------------------------------------------------------------
    // 4. GENOME PREPARATION AND MAPPING (branching logic)
    //----------------------------------------------------------------
    
    // Branch based on genome availability and assembly mode
    if (params.genome || params.accession) {
        // Option 1: Reference genome provided
        log.info "Using provided reference genome for mapping"
        
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
        
        // Map reads to reference
        map_reads( repair.out, genome_indexed )
        
        // Generate BAM statistics
        samtools_stats( map_reads.out )
        
        // MultiQC for mapping
        multiqc_mapping_out = multiqc_mapping(
            samtools_stats.out
                .map{ sid, stats, flagstats -> [stats, flagstats] }
                .flatten()
                .collect(),
            Channel.value('mapping')
        )
        
        // Create mapping summary
        samtools_summary(
            samtools_stats.out
                .map{ sid, stats, flagstats -> stats }
                .collect()
        )
        
        mapping_summary_ch = samtools_summary.out
        
        // No assembly stats for reference genome mode - create empty placeholder files inline
        assembly_stats_ch = Channel.fromPath("$projectDir/NO_ASSEMBLY", checkIfExists: false).ifEmpty { 
            file("NO_ASSEMBLY").text = "No assembly performed"
            return file("NO_ASSEMBLY")
        }
        filter_stats_ch = Channel.fromPath("$projectDir/NO_FILTER", checkIfExists: false).ifEmpty { 
            file("NO_FILTER").text = "No filtering performed"
            return file("NO_FILTER")
        }
        
    } else if (params.assembly_mode == "denovo") {
        // Option 2: De novo assembly workflow
        log.info "Performing de novo assembly from cleaned reads"
        
        // Run de novo assembly using repaired reads
        denovo_assembly( repair.out )
        
        // Capture assembly statistics
        assembly_stats_ch = denovo_assembly.out.assembly_stats
        filter_stats_ch = denovo_assembly.out.filter_stats
        
        // Use the de novo assembly as reference genome
        log.info "Indexing de novo assembly"
        prepare_genome_local( denovo_assembly.out.reference )
        genome_indexed = prepare_genome_local.out.genome
        
        // Map reads to de novo assembly
        map_reads( repair.out, genome_indexed )
        
        // Generate BAM statistics
        samtools_stats( map_reads.out )
        
        // MultiQC for mapping
        multiqc_mapping_out = multiqc_mapping(
            samtools_stats.out
                .map{ sid, stats, flagstats -> [stats, flagstats] }
                .flatten()
                .collect(),
            Channel.value('mapping_denovo')
        )
        
        // Create mapping summary
        samtools_summary(
            samtools_stats.out
                .map{ sid, stats, flagstats -> stats }
                .collect()
        )
        
        mapping_summary_ch = samtools_summary.out
        
    } else {
        // Option 3: No genome mode - just output cleaned reads
        log.info "No reference genome - outputting cleaned FASTQ files only"
        
        output_cleaned_reads( repair.out )
        
        // Create empty channels for mapping-related outputs
        multiqc_mapping_out = Channel.empty().mix(
            Channel.value("NO_MAPPING_HTML"),
            Channel.value("NO_MAPPING_STATS")
        )
        mapping_summary_ch = Channel.empty()  // Use empty channel instead of value
        
        // No assembly stats for no-genome mode - create empty placeholder files inline
        assembly_stats_ch = Channel.fromPath("$projectDir/NO_ASSEMBLY", checkIfExists: false).ifEmpty { 
            file("NO_ASSEMBLY").text = "No assembly performed"
            return file("NO_ASSEMBLY")
        }
        filter_stats_ch = Channel.fromPath("$projectDir/NO_FILTER", checkIfExists: false).ifEmpty { 
            file("NO_FILTER").text = "No filtering performed"
            return file("NO_FILTER")
        }
    }

    //----------------------------------------------------------------
    // 3b. SPECIES IDENTIFICATION (Optional)
    //----------------------------------------------------------------
    if (params.run_species_id) {
        log.info "Running species identification workflow"
        
        // Extract reads after 5' trimming for species ID
        species_id_input = fastp_trim_5.out
            .map{ sid, r1, r2, json, html -> tuple(sid, r1, r2) }
        
        // Run species identification subworkflow
        species_identification(species_id_input)
        
        // The results are available in:
        // - species_identification.out.species_report
        // - species_identification.out.species_consensus
        // - species_identification.out.species_stats
    }
    
    //----------------------------------------------------------------
    // 5. COLLECT STATS AND GENERATE REPORT
    //----------------------------------------------------------------
    
    // Collect all the general stats files based on mode
    if (params.assembly_mode != "denovo" && (!params.genome && !params.accession)) {
        // No mapping mode - collect stats without mapping
        all_stats = Channel.empty()
            .mix(
                multiqc_raw_out[1],
                multiqc_trim3_out[1],
                multiqc_clumpify_out[1],
                multiqc_trim5_out[1],
                multiqc_screen_out[1],
                multiqc_repair_out[1]
            )
            .collect()
    } else {
        // Include mapping stats
        all_stats = Channel.empty()
            .mix(
                multiqc_raw_out[1],
                multiqc_trim3_out[1],
                multiqc_clumpify_out[1],
                multiqc_trim5_out[1],
                multiqc_screen_out[1],
                multiqc_repair_out[1],
                multiqc_mapping_out[1],
                mapping_summary_ch
            )
            .collect()
    }
    
    // Run R analysis on collected stats
    analyze_read_stats(all_stats)
    
    //----------------------------------------------------------------
    // 6. GENERATE FINAL REPORT - WITH PLACEHOLDER PROCESS
    //----------------------------------------------------------------
    
    // Create a simple process to generate placeholder files when needed
    process create_assembly_placeholders {
        output:
            path "no_assembly_stats.txt"
            path "no_filter_stats.txt"
        
        script:
        """
        echo "No assembly performed" > no_assembly_stats.txt
        echo "No filtering performed" > no_filter_stats.txt
        """
    }
    
    // Prepare genome source information
    genome_source = params.genome 
        ? Channel.value("local:${params.genome}")
        : params.accession
        ? Channel.value("accession:${params.accession}")
        : params.assembly_mode == "denovo"
        ? Channel.value("denovo:assembled")
        : Channel.value("none:no_reference")
    
    // Collect all MultiQC HTML reports
    all_multiqc_reports = Channel.empty()
        .mix(
            multiqc_raw_out[0],
            multiqc_trim3_out[0],
            multiqc_clumpify_out[0],
            multiqc_trim5_out[0],
            multiqc_screen_out[0],
            multiqc_repair_out[0]
        )
        .collect()
    
    // Handle mapping summary for report generation
    // Use ifEmpty to provide a default value when no mapping was done
    mapping_summary_for_report = mapping_summary_ch.ifEmpty("NO_MAPPING")
    
    // Handle assembly stats - create placeholders if needed
    if (params.assembly_mode == "denovo") {
        final_assembly_stats = assembly_stats_ch
        final_filter_stats = filter_stats_ch
    } else {
        placeholder_outputs = create_assembly_placeholders()
        final_assembly_stats = placeholder_outputs[0]
        final_filter_stats = placeholder_outputs[1]
    }
    
    // Generate final report with assembly stats
    generate_report(
        analyze_read_stats.out[0],  // qc_summary_plot.png
        analyze_read_stats.out[1],  // read_counts_summary.txt
        all_multiqc_reports,
        genome_source,
        analyze_read_stats.out[5],  // initial_reads_histogram.png
        analyze_read_stats.out[6],  // mapped_reads_histogram.png
        mapping_summary_for_report,
        final_assembly_stats,        // Assembly statistics (actual or placeholder)
        final_filter_stats           // Filter statistics (actual or placeholder)
    )
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
// MODULE WORKFLOWS
//--------------------------------------------------------------------
include { species_identification } from './workflows/species_identification.nf'
include { denovo_assembly.nf } from './workflows/denovo_assembly.nf'

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
include { generate_report } from './modules/generate_report.nf'
include { output_cleaned_reads } from './modules/output_cleaned_reads.nf'

// de novo assembly modules
include { extract_unique_seqs } from './modules/extract_unique_seqs.nf'
include { filter_unique_seqs } from './modules/filter_unique_seqs.nf'
include { assemble_rainbow } from './modules/assemble_rainbow.nf'

// species id modules
include { get_mito_genes } from './modules/get_mito_genes.nf'
include { blast_mito_genes } from './modules/blast_mito_genes.nf'
include { summarize_species_id } from './modules/summarize_species_id.nf'

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