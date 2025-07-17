// modules/analyze_read_stats.nf
process analyze_read_stats {
    label 'r_analysis'
    tag "read_stats_analysis"
    
    publishDir "${params.outdir}/read_analysis", mode: 'copy'
    
    input:
        path(stats_files)
    
    output:
        path("qc_summary_plot.png")
        path("read_counts_summary.txt")
        path("model_summary.txt")
        path("analysis_log.txt")
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Setup logging
    sink("analysis_log.txt", split = TRUE)
    cat("=== QC Read Statistics Analysis ===\\n")
    cat("Starting at:", format(Sys.time()), "\\n\\n")
    
    # Load libraries
    suppressPackageStartupMessages({
        library(tidyverse)
        library(glmmTMB)
        library(emmeans)
        library(broom)
        library(scales)
    })
    
    # Helper function to read and process each file type
    process_file <- function(filepath) {
        cat("Processing:", basename(filepath), "\\n")
        
        # Determine stage from filename
        stage <- case_when(
            str_detect(filepath, "raw.*fastqc") ~ "raw",
            str_detect(filepath, "trim.*3") ~ "trim3", 
            str_detect(filepath, "clumpify") ~ "dedup",
            str_detect(filepath, "trim.*5") ~ "trim5",
            str_detect(filepath, "fastq.*screen") ~ "fqscrn",
            str_detect(filepath, "repair") ~ "repr",
            str_detect(filepath, "mapping_summary") ~ "map",
            TRUE ~ "unknown"
        )
        
        if (stage == "unknown") {
            cat("  Skipping unknown file type\\n")
            return(NULL)
        }
        
        # Read the file
        data <- read_delim(filepath, delim = "\\t", show_col_types = FALSE, na = c("", "NA"))
        
        # Process based on stage
        if (stage == "map") {
            # Mapping data: convert individual reads to read pairs
            result <- data %>%
                select(sample_id = Sample, n_reads = Mapped_Paired / 2) %>%
                mutate(stage = stage)
            cat("  Mapping data - converted to read pairs\\n")
            
        } else {
            # MultiQC data: use fastqc-total_sequences column
            fastqc_col <- "fastqc-total_sequences"
            if (!fastqc_col %in% names(data)) {
                cat("  Warning: No fastqc-total_sequences column found\\n")
                return(NULL)
            }
            
            result <- data %>%
                select(sample_id = Sample, n_reads = !!sym(fastqc_col)) %>%
                mutate(stage = stage)
            cat("  MultiQC data - using fastqc-total_sequences\\n")
        }
        
        cat("  Processed", nrow(result), "samples\\n")
        return(result)
    }
    
    # Process all files
    cat("Processing input files...\\n")
    file_list <- list.files(".", pattern = "\\.txt$", full.names = TRUE)
    file_list <- file_list[!str_detect(file_list, "analysis_log|debug_stats")]
    
    all_data <- map_dfr(file_list, process_file)
    
    if (is.null(all_data) || nrow(all_data) == 0) {
        cat("ERROR: No data processed\\n")
        quit(save = "no", status = 1)
    }
    
    # Clean and standardize the data
    cat("\\nCleaning data...\\n")
    clean_data <- all_data %>%
        # Filter for R1 files only (except mapping)
        filter(stage == "map" | str_ends(sample_id, ".1") | str_ends(sample_id, ".r1")) %>%
        # Convert read counts
        mutate(
            n_reads = case_when(
                stage == "map" ~ as.integer(round(n_reads)),  # Already converted to pairs
                TRUE ~ as.integer(round(n_reads * 1e6))       # Convert millions to actual count
            ),
            # Clean sample names
            sample_id = str_remove(sample_id, "_fp1.*") %>%
                       str_remove("\\.1$") %>%
                       str_remove("\\.r1$"),
            # Set stage as ordered factor
            stage = factor(stage, 
                          levels = c("raw", "trim3", "dedup", "trim5", "fqscrn", "repr", "map"), 
                          ordered = TRUE)
        ) %>%
        # Remove any rows with missing data
        filter(!is.na(sample_id), !is.na(n_reads), n_reads > 0)
    
    cat("Final data: ", nrow(clean_data), "rows,", 
        n_distinct(clean_data$sample_id), "samples,", 
        n_distinct(clean_data$stage), "stages\\n")
    
    # Save summary
    clean_data %>%
        arrange(sample_id, stage) %>%
        write_delim("read_counts_summary.txt", delim = "\\t")
    
    # Fit model
    cat("\\nFitting model...\\n")
    model <- glmmTMB(n_reads ~ stage + (1 | sample_id),
                     data = clean_data,
                     family = nbinom2(link = "log"))
    
    # Save model summary
    capture.output(summary(model)) %>%
        writeLines("model_summary.txt")
    
    # Generate plot data
    cat("Generating plot...\\n")
    em_results <- emmeans(model, ~ stage, type = "response")
    plot_data <- em_results %>%
        tidy(conf.int = TRUE) %>%
        mutate(stage = factor(stage, levels = levels(clean_data$stage)))
    
    # Create plot
    p <- ggplot(plot_data, aes(x = stage, y = response)) +
        # Individual sample trajectories
        geom_line(data = clean_data,
                 aes(y = n_reads, color = sample_id, group = sample_id),
                 alpha = 0.7, size = 0.5) +
        geom_point(data = clean_data,
                  aes(y = n_reads, color = sample_id),
                  size = 2, alpha = 0.8) +
        # Model estimates
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                     width = 0.2, size = 1) +
        geom_line(aes(group = 1), size = 1.5) +
        geom_point(size = 4) +
        # Formatting
        scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
        labs(
            title = "QC Pipeline Read Retention",
            subtitle = "Read pairs across processing stages",
            x = "Processing Stage",
            y = "Number of Read Pairs",
            color = "Sample"
        ) +
        theme_classic() +
        theme(
            plot.title = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    
    ggsave("qc_summary_plot.png", plot = p, width = 10, height = 8, dpi = 300)
    
    cat("\\nCompleted at:", format(Sys.time()), "\\n")
    sink()
    """
}