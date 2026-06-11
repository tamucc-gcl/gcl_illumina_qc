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

// Assembly parameters — see the GRID MODEL block below. cutoff1/cutoff2/
// cluster_similarity/final_similarity/div_f/merge_r are all defined there (each
// accepts a scalar to pin or a list to sweep). div_K is fixed:
params.div_K = 10        // Rainbow div -K (fixed; not a sweep axis)

// --- De novo OPTIMIZATION (grid sweep + r80-elbow selection) ---
params.do_optimize = false

// GRID MODEL (pivot): six assembly axes, each accepts EITHER a scalar (FIXED) or a
// list (SWEPT). The candidate grid is the Cartesian product of the resolved axes.
// Set any axis to a single number to pin it; set it to a list to sweep it.
//
// WHY this set: the cutoff axes (c1,c2) are a monotone STRINGENCY dial — every
// quality signal we tried is monotone in assembly SIZE along them, so they have no
// interior optimum to "find". The real STACKS-M-analogue tradeoff (collapse
// paralogs vs. split alleles) lives in the Rainbow div/merge + final CD-HIT step.
// So by DEFAULT we PIN the cutoffs at principled values and SWEEP the assembly
// params, then pick the r80-vs-n_contigs ELBOW (diminishing-returns point).
//
// cutoff1 default = a wide SWEPT band. The cutoffs are the axes that move assembly
//   SIZE the most, so sweeping them is how the grid SPANS the r80-vs-size plateau
//   (a narrow/pinned grid gave only a ~4% size range with no interior elbow). The
//   r80 elbow is still the SELECTOR; the cutoff sweep just generates the curve it
//   needs. Pin to a scalar (e.g. --cutoff1 5) to fix it; null => include the NB
//   crossover value in the swept set as well.
// NOTE: any axis set to a SCALAR is pinned (not swept). Lists are swept.
params.cutoff1            = [2, 3, 4, 5, 6, 7, 8]   // SWEPT wide (size span); scalar pins; null adds NB value
params.cutoff2            = [3, 4]                  // SWEPT low band (keeps the r80-interesting regime)
params.cluster_similarity = 0.8                     // INITIAL CD-HIT (loose pre-grouping; pinned — minor knob)
params.final_similarity  = [0.90, 0.95]            // FINAL CD-HIT (per-taxon precision merge -> swept)
params.div_f             = [0.1, 0.2, 0.5]         // Rainbow div -f (allele-split freq -> swept)
params.merge_r           = 2                        // Rainbow merge -r (pinned; list to sweep)

// Grid-size guardrail: warn if the Cartesian product exceeds this. The DEFAULT grid
// (c1[2..8] x c2[3,4] x div_f[.1,.2,.5] x final_sim[.90,.95]) = 84 candidates; each
// is one full assembly + one SNP pass, so this is the expensive case BY DESIGN (it
// spans the size range the r80 elbow needs). Pin axes to scalars to shrink it.
params.max_grid_candidates = 100

// SNP-signal sampling + r80 selection
params.snp_sample_pct  = 75              // fraction of samples mapped for the r80/SNP pass (tunable; r80 stability scales with this)
params.optimize_seed   = 42              // deterministic pseudo-rep split + SNP-sample selection
params.n_pseudo_reps   = 6               // individuals split 50/50 for the (reported) concordance signal; never published
params.r80_threshold   = 0.8             // STACKS r80: locus genotyped in >= this fraction of SNP-subset samples

// r80-elbow selector: the elbow is where marginal r80 gain per added contig drops
// below this threshold (loci per 1000 added contigs). Below the elbow you are
// paying contigs (paralogs/junk) for ~no new broadly-shared loci. Smaller =>
// more permissive (larger reference); larger => more parsimonious (smaller).
params.r80_elbow_min_slope = 30          // loci gained per 1000 added contigs at the elbow

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
    // 1a. PRE-FLIGHT SUMMARY — print resolved inputs BEFORE heavy work submits,
    //     so a mis-parsed --reads glob or an oversized/mis-parsed grid is caught
    //     at launch rather than after dozens of assemblies are queued.
    //----------------------------------------------------------------
    // Count matched read FILES at parse time via a Groovy glob (does NOT touch/
    // consume raw_reads_pairs — a queue channel can only be consumed once, and the
    // pipeline needs it downstream). files() ALWAYS returns a list (file() does not
    // for a single match), so the count is robust. Pairs ~= files / 2.
    def pf_nfiles = files(params.reads).size()
    log.info "PRE-FLIGHT | --reads '${params.reads}' matched ${pf_nfiles} file(s) (~${(pf_nfiles/2) as int} sample pair(s))"
    if (pf_nfiles == 0)
        log.warn "PRE-FLIGHT | NO files matched --reads. Check the glob (needs a sample-ID wildcard + paired {F,R} or {1,2}), e.g. 'data/fq_raw/*.{F,R}.fq.gz'"
    else if (pf_nfiles % 2 != 0)
        log.warn "PRE-FLIGHT | --reads matched an ODD number of files (${pf_nfiles}); paired-end expects an even count. Check for an unpaired/stray file."

    if (params.do_optimize) {
        // Mirror optimize_denovo's axis normalizer to report the resolved grid up front.
        def pfList = { v ->
            if (v == null) return [null]
            if (v instanceof List) return v
            if (v instanceof CharSequence) {
                def s = v.toString().trim().replaceAll(/^\[|\]$/, '').trim()
                if (s == '') return [null]
                if (s.contains(',') || s.contains(' '))
                    return s.split(/[,\s]+/).findAll { it }.collect { it.trim() }
                return [s]
            }
            return [v]
        }
        def pf_c1   = pfList(params.cutoff1)
        def pf_c2   = pfList(params.cutoff2)
        def pf_isim = pfList(params.cluster_similarity)
        def pf_fsim = pfList(params.final_similarity)
        def pf_divf = pfList(params.div_f)
        def pf_mr   = pfList(params.merge_r)
        // c1=null contributes the (runtime) NB value = 1 value for counting purposes.
        def c1_n = pf_c1.collect { it == null ? 'NB' : it }.unique().size()
        def n_grid = c1_n * pf_c2.size() * pf_isim.size() * pf_fsim.size() * pf_divf.size() * pf_mr.size()
        def gmax = (params.max_grid_candidates ?: 100) as int
        log.info "PRE-FLIGHT | OPTIMIZE grid axes: cutoff1=${pf_c1} cutoff2=${pf_c2} init_sim=${pf_isim} final_sim=${pf_fsim} div_f=${pf_divf} merge_r=${pf_mr}"
        log.info "PRE-FLIGHT | grid = ${n_grid} candidate(s) (each = 1 assembly + 1 SNP pass); guardrail max_grid_candidates=${gmax}"
        if (n_grid > gmax)
            log.warn "PRE-FLIGHT | grid (${n_grid}) EXCEEDS max_grid_candidates (${gmax}) — pin axes to scalars or raise the guardrail."
        log.info "PRE-FLIGHT | SNP/r80 pass maps ${params.snp_sample_pct}% of samples per candidate (anchor: ${ (params.genome_size_est && params.size_select_min && params.size_select_max) ? 'ENABLED' : 'disabled' })"
    }

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