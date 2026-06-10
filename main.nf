// main.nf – GCL Illumina QC pipeline
// July 2025 – consolidated MultiQC per step
// Oct 2025 - implement local genome option
// Nov 2025 - add de novo assembly & species ID options
// Mar 2026 - add stacks sequencing type with uniform trim

nextflow.enable.dsl = 2

//--------------------------------------------------------------------
// USER PARAMETERS
//--------------------------------------------------------------------
params.outdir      = "illumina_qc"
params.publish_dir_mode = "copy"  // Options: "copy", "symlink", "move"

params.reads       = "data/fq_raw/*.{1,2}.fq.gz"    // paired‑end,  sampleID.1.fq.gz / .2.fq.gz
params.accession   = null                            // NCBI assembly accession (optional)
params.genome      = null                            // Path to local genome file (optional)
params.assembly_mode = "none"                     // Options: "none" (no assembly), "denovo" (assemble from reads)
params.decontam_conffile    = "${projectDir}/config_files/example_fastq-screen.conf"  // FastQ Screen config file
params.sequencing_type = "whole_genome"  // Options: "ddrad", "whole_genome", or "stacks"

// Uniform trim parameter - set automatically for stacks mode, or manually override
params.uniform_trim_length = 0  // 0 = disabled, >0 = trim all reads to this length (bp)

// Species ID parameters
params.run_species_id = false  // Enable/disable species identification
params.mito_reference = "${projectDir}/databases/mito_gene_refs.fasta" // Mitochondrial gene reference provided with gcl_illumina_qc compiled from https://github.com/cmayer/MitoGeneExtractor
params.genetic_code = 2  //Mitochondrial genetic code: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
params.blast_db = "/work/birdlab/databases/midori2_latest/CO1/midori2_latest"  // Path to local BLAST database
params.taxonomy_db = "/work/birdlab/databases/ncbi_taxonomy" // Path to local NCBI Taxonomy database

// Assembly parameters
params.cutoff1 = null       // Min reads per individual (null = auto-detect from data)
params.cutoff2 = null       // Min individuals (null = auto-detect from data)  
params.cluster_similarity = null // used when do_optimize = false; ignored during optimization
params.div_f = 0.5
params.div_K = 10
params.merge_r = 2
params.final_similarity = 0.9

// --- De novo cutoff/similarity optimization ---
params.do_optimize = false

// Grid to search. Knee values from diagnostics act as CEILINGS for c1/c2 unless
// pinned. Explicit --cutoff1/--cutoff2/--cluster_similarity PIN that dimension.
params.optimize_cluster_similarity = [0.85, 0.90, 0.95]
params.cutoff1_floor               = 2        // lowest cutoff1 (per-individual coverage) in the grid
params.cutoff2_floor               = 3        // lowest cutoff2 (n individuals); 2 = junk/bloat + slow CD-HIT, skipped by default

// Subset used to BUILD candidate references (assembly branch, intact individuals).
params.snp_sample_pct = 25               // fraction of samples used for STAGE-2 SNP-recovery scoring (NOT assembly; assembly uses all samples)
params.optimize_seed       = 42               // deterministic pseudo-rep split + stage-2 SNP-sample selection

// Pseudo-replicates for the concordance signal (stage 2). The N highest-depth
// individuals are auto-selected and split 50/50; halves are used ONLY for
// concordance and never published (the intact individual goes to production).
params.n_pseudo_reps = 6

// Cheap-signal: coverage-uniformity subset size (samples mapped to each candidate for CV/Gini)
params.cv_sample_n = 4

// Two-stage handoff: cheap signals rank all candidates; the top-N by provisional
// rank get the expensive bcftools step. N is auto-set by the largest gap in the
// provisional-rank-score curve, clamped to [min,max] as a compute safety valve.
params.stage2_min_candidates = 3
params.stage2_max_candidates = 10

// --- Optional biological locus-count anchor (signal 2) ---
// If ALL THREE are supplied, expected RAD locus count is estimated and proximity
// to it becomes one of the ranking signals. If ANY is null, signal 2 is skipped.
params.enzyme_pair       = null   // e.g. "SbfI-EcoRI" (informational; model uses recognition-site lengths below)
params.genome_size_est   = null   // estimated genome size in bp, e.g. 1.3e9
params.size_select_min   = null   // INSERT lower bound in bp (genomic DNA between cut sites; NOT the fragment-trace window, which includes adapters), e.g. 150
params.size_select_max   = null   // INSERT upper bound in bp, e.g. 221
// Recognition-site lengths (bp). Model auto-assigns the LARGER as the rare cutter.
// e.g. SbfI=8, EcoRI=6. Order does not matter.
params.enzyme1_site_len  = 6
params.enzyme2_site_len  = 4

//--------------------------------------------------------------------
// DERIVED PARAMETERS
//--------------------------------------------------------------------

// Stacks mode behaves like ddrad but with uniform trimming
// Resolve the effective sequencing type for tools that care about ddrad vs whole_genome
def effective_sequencing_type = params.sequencing_type == "stacks" ? "ddrad" : params.sequencing_type

// Resolve uniform trim length: explicit override takes priority, otherwise auto-set for stacks
def uniform_trim_length = params.uniform_trim_length > 0 
    ? params.uniform_trim_length 
    : (params.sequencing_type == "stacks" ? 140 : 0)

//--------------------------------------------------------------------
// WORKFLOW DEFINITION
//--------------------------------------------------------------------

workflow {
    // Validate sequencing type
    if (!(params.sequencing_type in ['ddrad', 'whole_genome', 'stacks'])) {
        error "Invalid sequencing_type '${params.sequencing_type}'. Must be 'ddrad', 'whole_genome', or 'stacks'"
    }

    // Log uniform trim info
    if (uniform_trim_length > 0) {
        log.info "Uniform trim enabled: all reads will be trimmed to ${uniform_trim_length} bp"
    }
    if (params.sequencing_type == "stacks") {
        log.info "Stacks mode: using ddRAD parameters with uniform read length (${uniform_trim_length} bp)"
    }

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
    
    // Step 2: Clumpify - use effective_sequencing_type so stacks behaves like ddrad
    clumpify( fastp_trim_3.out,
              Channel.value(effective_sequencing_type))
    
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
    
    // Step 3: 5' trimming with fastp - pass uniform trim length
    fastp_trim_5( 
        clumpify.out.map{ sid, r1, r2, stats -> tuple(sid, r1, r2) },
        Channel.value(uniform_trim_length)
    )
    
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
    // 3a. OUTPUT CLEANED READS (always runs)
    //----------------------------------------------------------------
    output_cleaned_reads( repair.out )

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
        map_reads( repair.out, genome_indexed.first() )
        
        // Generate BAM statistics
        samtools_stats( map_reads.out )
        
        // MultiQC for mapping
        multiqc_mapping_out = multiqc_mapping(
            samtools_stats.out[0]
                .map{ sid, stats, flagstats -> [stats, flagstats] }
                .flatten()
                .collect(),
            Channel.value('mapping')
        )
        
        // Create mapping summary
        samtools_summary(
            samtools_stats.out[0]
                .map{ sid, stats, flagstats -> stats }
                .collect(),
            samtools_stats.out[1]
                .map{ sid, soft_clip -> soft_clip }
                .collect(),
            samtools_stats.out[2]
                .map{ sid, aln_score -> aln_score }
                .collect()
        )
        
        mapping_summary_ch = samtools_summary.out
        
        // No assembly stats for reference genome mode - create empty placeholder files inline
        assembly_stats_ch = Channel.empty()
        filter_stats_ch = Channel.empty()
        cutoff1_plot_ch = Channel.empty()
        cutoff2_plot_ch = Channel.empty()
        sweep_plot_ch = Channel.empty()
        sweep_summary_ch = Channel.empty()
        
    } else if (params.assembly_mode == "denovo") {
        // Option 2: De novo assembly workflow
        log.info "Performing de novo assembly from cleaned reads"
        
        // Run de novo assembly using repaired reads
        denovo_assembly( repair.out )
        
        // Capture assembly statistics
        assembly_stats_ch = denovo_assembly.out.assembly_stats
        filter_stats_ch = denovo_assembly.out.filter_stats
        cutoff1_plot_ch = denovo_assembly.out.cutoff1_plot
        cutoff2_plot_ch = denovo_assembly.out.cutoff2_plot
        sweep_plot_ch = denovo_assembly.out.optimize_plot
        sweep_summary_ch = denovo_assembly.out.optimize_summary
        
        // Use the de novo assembly as reference genome
        log.info "Indexing de novo assembly"
        prepare_genome_local( denovo_assembly.out.reference )
        genome_indexed = prepare_genome_local.out.genome
        
        // Map reads to de novo assembly
        map_reads( repair.out, genome_indexed.first() )
        
        // Generate BAM statistics
        samtools_stats( map_reads.out )
        
        // MultiQC for mapping
        multiqc_mapping_out = multiqc_mapping(
            samtools_stats.out[0]
                .map{ sid, stats, flagstats -> [stats, flagstats] }
                .flatten()
                .collect(),
            Channel.value('mapping_denovo')
        )
        
        // Create mapping summary
        samtools_summary(
            samtools_stats.out[0]
                .map{ sid, stats, flagstats -> stats }
                .collect(),
            samtools_stats.out[1]
                .map{ sid, soft_clip -> soft_clip }
                .collect(),
            samtools_stats.out[2]
                .map{ sid, aln_score -> aln_score }
                .collect()
        )
        
        mapping_summary_ch = samtools_summary.out
        
    } else {
        // Option 3: No genome mode - just output cleaned reads
        log.info "No reference genome - outputting cleaned FASTQ files only"
        
        // Create empty channels for mapping-related outputs
        multiqc_mapping_out = Channel.empty().mix(
            Channel.value("NO_MAPPING_HTML"),
            Channel.value("NO_MAPPING_STATS")
        )
        mapping_summary_ch = Channel.empty()  // Use empty channel instead of value
        
        // No assembly stats for no-genome mode - create empty placeholder files inline
        assembly_stats_ch = Channel.empty()
        filter_stats_ch = Channel.empty()
        cutoff1_plot_ch = Channel.empty()
        cutoff2_plot_ch = Channel.empty()
        sweep_plot_ch = Channel.empty()
        sweep_summary_ch = Channel.empty()
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
        
        // Capture outputs for report
        species_blast_tsv = species_identification.out.combined_blast
        species_raw_pie = species_identification.out.raw_pie_chart
        species_summary_pie = species_identification.out.summary_pie_chart
        species_top_hits = species_identification.out.top_hits
        
    } else {
        // Create empty channels for species ID outputs
        species_blast_tsv = Channel.empty()
        species_raw_pie = Channel.empty()
        species_summary_pie = Channel.empty()
        species_top_hits = Channel.empty()
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
    process create_placeholders {
        output:
            path "no_assembly_stats.txt"
            path "no_filter_stats.txt"
            path "no_insert_size_violin.png"
            path "no_species_blast.txt"
            path "no_species_raw_pie.png"
            path "no_species_summary_pie.png"
            path "no_species_top_hits.csv"
            path "no_soft_clip_violin.png"    // [7] NEW
            path "no_aln_score_violin.png"    // [8] NEW
            path "no_cutoff1_curve.png"       // [9] NEW
            path "no_cutoff2_curve.png"       // [10] NEW
            path "no_sweep_plot.png"          // [11] NEW
            path "no_sweep_summary.tsv"       // [12] NEW
        
        script:
        """
        echo "No assembly performed" > no_assembly_stats.txt
        echo "No filtering performed" > no_filter_stats.txt
        echo "No mapping performed" > no_insert_size_violin.png
        echo "No species identification performed" > no_species_blast.txt
        echo "No species identification performed" > no_species_raw_pie.png
        echo "No species identification performed" > no_species_summary_pie.png
        echo "No species identification performed" > no_species_top_hits.csv
        echo "No mapping performed"             > no_soft_clip_violin.png
        echo "No mapping performed"             > no_aln_score_violin.png
        echo "No assembly performed"            > no_cutoff1_curve.png
        echo "No assembly performed"            > no_cutoff2_curve.png
        echo "No sweep performed"               > no_sweep_plot.png
        echo "No sweep performed"               > no_sweep_summary.tsv
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
    //mapping_summary_for_report = mapping_summary_ch.ifEmpty("NO_MAPPING")
    mapping_summary_for_report = mapping_summary_ch[0].ifEmpty {
            def f = file("${workDir}/NO_MAPPING_SUMMARY.txt")
            f.text = "No mapping performed"
            return f
        }

    // Create placeholders once for any non-actual outputs
    placeholder_outputs = create_placeholders()

    // Handle assembly stats
    if (params.assembly_mode == "denovo") {
        final_assembly_stats = assembly_stats_ch
        final_cutoff1_plot = cutoff1_plot_ch
        final_cutoff2_plot = cutoff2_plot_ch
        if (params.do_optimize) {
            final_filter_stats  = placeholder_outputs[1]
            final_sweep_plot    = sweep_plot_ch
            final_sweep_summary = sweep_summary_ch
        } else {
            final_filter_stats  = filter_stats_ch
            final_sweep_plot    = placeholder_outputs[11]
            final_sweep_summary = placeholder_outputs[12]
        }
    } else {
        final_assembly_stats = placeholder_outputs[0]
        final_filter_stats = placeholder_outputs[1]
        final_cutoff1_plot = placeholder_outputs[9]
        final_cutoff2_plot = placeholder_outputs[10]
        final_sweep_plot = placeholder_outputs[11]
        final_sweep_summary = placeholder_outputs[12]
    }
    
    // Handle insert size violin
    if (params.genome || params.accession || params.assembly_mode == "denovo") {
        final_insert_size_violin = mapping_summary_ch[1]
        final_soft_clip_violin     = samtools_summary.out[2]
        final_aln_score_violin     = samtools_summary.out[3]
    } else {
        final_insert_size_violin = placeholder_outputs[2]
        final_soft_clip_violin     = placeholder_outputs[7]   // new
        final_aln_score_violin     = placeholder_outputs[8]   // new
    }

    // Handle species ID outputs
    if (params.run_species_id) {
        final_species_blast = species_blast_tsv
        final_species_raw_pie = species_raw_pie
        final_species_summary_pie = species_summary_pie
        final_species_top_hits = species_top_hits
    } else {
        final_species_blast = placeholder_outputs[3]
        final_species_raw_pie = placeholder_outputs[4]
        final_species_summary_pie = placeholder_outputs[5]
        final_species_top_hits = placeholder_outputs[6]
    }
    
    // Generate final report with assembly stats and species ID
    generate_report(
        analyze_read_stats.out[0],  // qc_summary_plot.png
        analyze_read_stats.out[1],  // read_counts_summary.txt
        all_multiqc_reports,
        genome_source,
        analyze_read_stats.out[5],  // initial_reads_histogram.png
        analyze_read_stats.out[6],  // mapped_reads_histogram.png
        mapping_summary_for_report,
        final_insert_size_violin,    // insert_size_violin.png
        final_soft_clip_violin,      // NEW
        final_aln_score_violin,      // NEW
        final_assembly_stats,        // Assembly statistics (actual or placeholder)
        final_filter_stats,          // Filter statistics (actual or placeholder)
        final_species_blast,         // BLAST results TSV
        final_species_raw_pie,       // Raw BLAST pie chart
        final_species_summary_pie,   // Summary BLAST pie chart
        final_species_top_hits,       // Top BLAST hits CSV
        final_cutoff1_plot,          // NEW: de novo cutoff1 diagnostic plot
        final_cutoff2_plot,           // NEW: de novo cutoff2 diagnostic plot
        final_sweep_plot,            // NEW: cluster_similarity sweep comparison plot
        final_sweep_summary          // NEW: cluster_similarity sweep ranked table
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
include { denovo_assembly } from './workflows/denovo_assembly.nf'

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