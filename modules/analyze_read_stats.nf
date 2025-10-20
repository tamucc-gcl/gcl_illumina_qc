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
        path("stage_comparison.txt")
        path("sample_trajectories.png")
        path("debug_all_data.txt")
    
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
        library(ggrepel)
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
            # Use backtick notation to avoid Nextflow variable interpretation
            result <- data %>%
                select(sample_id = `Sample`, `Mapped_Paired`, `Properly_Paired`) %>%
                pivot_longer(cols = ends_with('Paired),
                             names_to = 'stage',
                             values_to = 'n_reads') %>% 
                mutate(
                    n_reads = n_reads / 2,  # Convert to read pairs
                    stage = case_when(stage == "Mapped_Paired" ~ "map",
                                      stage == "Properly_Paired" ~ "pp")
                )
            cat("  Mapping data - converted to read pairs\\n")
            
        } else {
            # MultiQC data: use fastqc-total_sequences column
            fastqc_col <- "fastqc-total_sequences"
            
            # Check if column exists
            if (!fastqc_col %in% names(data)) {
                cat("  Warning: No fastqc-total_sequences column found\\n")
                cat("  Available columns:", paste(names(data), collapse = ", "), "\\n")
                return(NULL)
            }
            
            # Check for Sample column
            if (!"Sample" %in% names(data)) {
                cat("  Warning: No Sample column found\\n")
                return(NULL)
            }
            
            result <- data %>%
                select(sample_id = `Sample`, n_reads = all_of(fastqc_col)) %>%
                mutate(stage = stage)
            cat("  MultiQC data - using fastqc-total_sequences\\n")
        }
        
        cat("  Processed", nrow(result), "samples\\n")
        return(result)
    }
    
    # Process all files
    cat("\\nProcessing input files...\\n")
    file_list <- list.files(".", pattern = "\\\\.txt\$", full.names = TRUE)
    file_list <- file_list[!str_detect(file_list, "analysis_log|debug_stats|read_counts_summary|model_summary|stage_comparison|sample_trajectories")]
    
    cat("Found", length(file_list), "input files\\n")
    
    all_data <- map_dfr(file_list, process_file)
    write_delim(all_data, "debug_all_data.txt", delim = "\\t")

    if (is.null(all_data) || nrow(all_data) == 0) {
        cat("ERROR: No data processed\\n")
        quit(save = "no", status = 1)
    }
    
    # Clean and standardize the data
    cat("\\nCleaning data...\\n")
    clean_data <- all_data %>%
        # Filter for R1 files only (except mapping)
        filter(stage == "map" | str_ends(sample_id, ".1") | str_ends(sample_id, ".r1") | str_ends(sample_id, "_R1")) %>%
        # Convert read counts and clean sample names
        mutate(
            n_reads = case_when(
                stage %in% c("map", "pp") ~ as.integer(round(n_reads)),  # Already converted to pairs
                TRUE ~ as.integer(round(n_reads * 1e6))       # Convert millions to actual count
            )
        ) %>%
        # Clean sample names in separate step to avoid regex issues
        mutate(
            sample_id = str_split_i(sample_id, "_fp1", 1),
            sample_id = str_remove(sample_id, "\\\\.1\$"),
            sample_id = str_remove(sample_id, "\\\\.r1\$"),
            sample_id = str_remove(sample_id, "_R1\$")
        ) %>%
        mutate(
            # Set stage as ordered factor
            stage = factor(stage, 
                          levels = c("raw", "trim3", "dedup", "trim5", "fqscrn", "repr", "map", "pp"), 
                          ordered = TRUE)
        ) %>%
        # Remove any rows with missing data
        filter(!is.na(sample_id), !is.na(n_reads), n_reads > 0) %>%
        # Remove duplicates if any
        distinct(sample_id, stage, .keep_all = TRUE)
    
    cat("Final data: ", nrow(clean_data), "rows,", 
        n_distinct(clean_data\$sample_id), "samples,", 
        n_distinct(clean_data\$stage), "stages\\n")
    
    # Save summary
    clean_data %>%
        arrange(sample_id, stage) %>%
        write_delim("read_counts_summary.txt", delim = "\\t")
    
    # Calculate stage-to-stage retention rates
    cat("\\nCalculating retention rates...\\n")
    retention_rates <- clean_data %>%
        arrange(sample_id, stage) %>%
        group_by(sample_id) %>%
        mutate(
            prev_reads = lag(n_reads),
            retention_pct = ifelse(!is.na(prev_reads), n_reads / prev_reads * 100, NA),
            prev_stage = lag(stage)
        ) %>%
        filter(!is.na(retention_pct)) %>%
        mutate(transition = paste(prev_stage, "->", stage)) %>%
        group_by(transition) %>%
        summarise(
            mean_retention = mean(retention_pct),
            sd_retention = sd(retention_pct),
            min_retention = min(retention_pct),
            max_retention = max(retention_pct),
            n_samples = n()
        )
    
    write_delim(retention_rates, "stage_comparison.txt", delim = "\\t")
    
    # Fit model
    cat("\\nFitting negative binomial model...\\n")
    
    # Check if we have enough data
    if (n_distinct(clean_data\$sample_id) < 2) {
        cat("WARNING: Only one sample found, skipping mixed model\\n")
        # Simple model without random effects
        model <- glm(n_reads ~ stage, data = clean_data, family = negative.binomial(theta = 1))
    } else {
        model <- glmmTMB(n_reads ~ stage + (1 | sample_id),
                         data = clean_data,
                         family = nbinom2(link = "log"))
    }
    
    # Save model summary
    capture.output(summary(model)) %>%
        writeLines("model_summary.txt")
    
    # Generate plot data
    cat("\\nGenerating plots...\\n")
    
    # Get model predictions
    if (n_distinct(clean_data\$sample_id) < 2) {
        # For single sample, just use the data
        plot_data <- clean_data %>%
            group_by(stage) %>%
            summarise(
                response = mean(n_reads),
                conf.low = n_reads,
                conf.high = n_reads
            )
    } else {
        em_results <- emmeans(model, ~ stage, type = "response")
        plot_data <- em_results %>%
            tidy(conf.int = TRUE) %>%
            mutate(stage = factor(stage, levels = levels(clean_data\$stage)))
    }
    
    # Main summary plot
    p1 <- ggplot(plot_data, aes(x = stage, y = response)) +
        # Individual sample trajectories
        geom_line(data = clean_data,
                 aes(y = n_reads, color = sample_id, group = sample_id),
                 alpha = 0.5, linewidth = 0.5) +
        geom_point(data = clean_data,
                  aes(y = n_reads, color = sample_id),
                  size = 2, alpha = 0.7) +
        # Model estimates
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                     width = 0.2, linewidth = 1, color = "black") +
        geom_line(aes(group = 1), linewidth = 1.5, color = "black") +
        geom_point(size = 4, color = "black") +
        # Formatting
        scale_y_log10(
            labels = trans_format("log10", math_format(10^.x)),
            breaks = trans_breaks("log10", function(x) 10^x)
        ) +
        labs(
            title = "QC Pipeline Read Retention",
            subtitle = paste("Read pairs across processing stages (n =", n_distinct(clean_data\$sample_id), "samples)"),
            x = "Processing Stage",
            y = "Number of Read Pairs (log scale)",
            color = "Sample"
        ) +
        theme_classic() +
        theme(
            plot.title = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = if(n_distinct(clean_data\$sample_id) > 10) "none" else "right"
        )
    
    ggsave("qc_summary_plot.png", plot = p1, width = 10, height = 8, dpi = 300)
    
    # Individual sample trajectory plot with annotations
    p2 <- clean_data %>%
        group_by(sample_id) %>%
        mutate(
            final_reads = n_reads[stage == max(stage)],
            initial_reads = n_reads[stage == min(stage)],
            overall_retention = final_reads / initial_reads * 100
        ) %>%
        ggplot(aes(x = stage, y = n_reads, group = sample_id)) +
        geom_line(aes(color = overall_retention), linewidth = 1) +
        geom_point(aes(color = overall_retention), size = 2) +
        # Add sample labels at the end
        geom_text_repel(
            data = . %>% filter(stage == max(stage)),
            aes(label = paste0(sample_id, " (", round(overall_retention, 1), "%)")),
            hjust = -0.1, size = 3, direction = "y"
        ) +
        scale_y_log10(
            labels = trans_format("log10", math_format(10^.x)),
            breaks = trans_breaks("log10", function(x) 10^x)
        ) +
        scale_color_gradient2(
            low = "red", mid = "yellow", high = "green",
            midpoint = 50, limits = c(0, 100),
            name = "Overall\\nRetention %"
        ) +
        labs(
            title = "Individual Sample Trajectories",
            subtitle = "Colored by overall retention percentage",
            x = "Processing Stage",
            y = "Number of Read Pairs (log scale)"
        ) +
        theme_classic() +
        theme(
            plot.title = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        expand_limits(x = c(1, length(levels(clean_data\$stage)) + 1))
    
    ggsave("sample_trajectories.png", plot = p2, width = 12, height = 8, dpi = 300)
    
    # Print summary statistics
    cat("\\n=== Summary Statistics ===\\n")
    cat("Total samples analyzed:", n_distinct(clean_data\$sample_id), "\\n")
    cat("Processing stages:", paste(levels(clean_data\$stage), collapse = ", "), "\\n\\n")
    
    # Overall retention statistics
    overall_stats <- clean_data %>%
        group_by(sample_id) %>%
        filter(stage %in% c(min(stage), max(stage))) %>%
        summarise(
            initial = n_reads[stage == min(stage)],
            final = n_reads[stage == max(stage)],
            retention_pct = final / initial * 100
        )
    
    cat("Overall retention (first to last stage):\\n")
    cat("  Mean:", round(mean(overall_stats\$retention_pct), 2), "%\\n")
    cat("  SD:", round(sd(overall_stats\$retention_pct), 2), "%\\n")
    cat("  Range:", round(min(overall_stats\$retention_pct), 2), "-", 
        round(max(overall_stats\$retention_pct), 2), "%\\n")
    
    cat("\\nCompleted at:", format(Sys.time()), "\\n")
    sink()
    """
}