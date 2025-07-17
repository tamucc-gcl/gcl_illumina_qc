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
    
    # Capture all output for logging
    sink("analysis_log.txt", split = TRUE)
    
    cat("=== QC Read Statistics Analysis ===\\n")
    cat("Starting analysis at:", format(Sys.time()), "\\n\\n")
    
    #### Libraries ####
    cat("Loading required libraries...\\n")
    
    # Initialize flags
    use_multcomp <- FALSE
    
    suppressPackageStartupMessages({
        library(tidyverse)
        library(glmmTMB)
        library(emmeans)
        library(broom)
        
        # Re-attach dplyr after other packages to avoid MASS conflicts
        library(dplyr, warn.conflicts = FALSE)
        
        # Try to load optional packages
        if (!require(multcomp, quietly = TRUE)) {
            cat("multcomp not available, using basic comparisons\\n")
        } else {
            use_multcomp <- TRUE
        }
        
        if (!require(janitor, quietly = TRUE)) {
            cat("janitor not available, proceeding without it\\n")
        }
        
        if (!require(patchwork, quietly = TRUE)) {
            cat("patchwork not available, proceeding without it\\n")
        }
        
        if (!require(ggtext, quietly = TRUE)) {
            cat("ggtext not available, proceeding without it\\n")
        }
    })
    cat("Libraries loaded successfully\\n\\n")
    
    #### Data Processing ####
    cat("Processing input data files...\\n")
    
    # List all available files
    cat("Available files:\\n")
    available_files <- list.files(".", full.names = TRUE)
    available_files <- available_files[str_ends(available_files, ".txt")]
    cat(paste(available_files, collapse = "\\n"), "\\n\\n")
    
    # Process the data with explicit namespace calls
    data <- available_files %>%
        tibble(filepath = .) %>%
        dplyr::mutate(stage = case_when(
            str_detect(filepath, "raw.*fastqc") ~ "raw",
            str_detect(filepath, "trim.*3") ~ "trim3", 
            str_detect(filepath, "clumpify") ~ "dedup",
            str_detect(filepath, "trim.*5") ~ "trim5",
            str_detect(filepath, "fastq.*screen") ~ "fqscrn",
            str_detect(filepath, "repair") ~ "repr",
            str_detect(filepath, "mapping_summary") ~ "map",
            TRUE ~ "unknown"
        )) %>%
        # Print stage assignments for debugging
        {
            cat("Stage assignments:\\n")
            for(i in 1:nrow(.)) {
                current_row <- .[i,]
                cat("File:", current_row[["filepath"]], "-> Stage:", current_row[["stage"]], "\\n")
            }
            cat("\\n")
            .
        } %>%
        dplyr::filter(stage != "unknown") %>%
        dplyr::mutate(stage = factor(stage, 
                            levels = c("raw", "trim3", "dedup", "trim5", 
                                     "fqscrn", "repr", "map"), 
                            ordered = TRUE)) %>%
        # Process each file individually
        {
            all_data <- tibble()
            for(i in 1:nrow(.)) {
                current_row <- .[i,]
                current_file <- current_row[["filepath"]]
                current_stage <- current_row[["stage"]]
                
                cat("Processing file:", current_file, "as stage:", current_stage, "\\n")
                
                tryCatch({
                    # Read the file
                    file_data <- read_delim(current_file, delim = "\\t", 
                                          show_col_types = FALSE, 
                                          na = c("", "NA", "N/A"))
                    
                    # Print column names for debugging
                    cat("Columns in", current_file, ":", paste(names(file_data), collapse = ", "), "\\n")
                    
                    if(current_stage == "map") {
                        # For mapping summary - use explicit dplyr:: namespace
                        processed_data <- file_data %>%
                            dplyr::select(sample_id = Sample, n_reads = Mapped_Paired) %>%
                            dplyr::mutate(stage = current_stage)
                    } else {
                        # For MultiQC files, look for fastqc total sequences column
                        fastqc_col <- names(file_data)[str_detect(names(file_data), "fastqc.*total.*sequences")]
                        
                        if(length(fastqc_col) == 0) {
                            # If no fastqc column found, look for other sequence-related columns
                            fastqc_col <- names(file_data)[str_detect(names(file_data), "sequences|reads")]
                        }
                        
                        if(length(fastqc_col) == 0) {
                            cat("Warning: No suitable reads column found in", current_file, "\\n")
                            next
                        }
                        
                        cat("Using reads column:", fastqc_col[1], "for file:", current_file, "\\n")
                        
                        processed_data <- file_data %>%
                            dplyr::select(sample_id = Sample, n_reads = !!fastqc_col[1]) %>%
                            dplyr::mutate(stage = current_stage)
                    }
                    
                    # Add to combined data
                    all_data <- dplyr::bind_rows(all_data, processed_data)
                    
                }, error = function(e) {
                    cat("Error processing file:", current_file, "- Error:", conditionMessage(e), "\\n")
                })
            }
            all_data
        } %>%
        # Convert stage back to factor with explicit namespace
        dplyr::mutate(stage = factor(stage, 
                            levels = c("raw", "trim3", "dedup", "trim5", 
                                     "fqscrn", "repr", "map"), 
                            ordered = TRUE)) %>%
        # Remove rows with missing data
        dplyr::filter(!is.na(sample_id), !is.na(n_reads)) %>%
        # Filter for R1 files only (except mapping) and clean sample names
        dplyr::filter(stage == 'map' | str_ends(sample_id, ".1") | str_ends(sample_id, ".r1")) %>%
        dplyr::mutate(n_reads = if_else(stage != 'map', 1e6 * n_reads, n_reads) %>%
                       round(0) %>%
                       as.integer(),
               # Clean sample names
               sample_id = str_remove(sample_id, "_fp1.*") %>%
                          str_remove(fixed(".1")) %>%
                          str_remove(fixed(".r1")))
    
    cat("Data processing completed\\n")
    cat("Samples found:", paste(unique(data[["sample_id"]]), collapse = ", "), "\\n")
    cat("Stages found:", paste(levels(data[["stage"]]), collapse = ", "), "\\n")
    cat("Number of rows:", nrow(data), "\\n\\n")
    
    # Check if we have any data
    if(nrow(data) == 0) {
        cat("ERROR: No data found! Check file formats and column names.\\n")
        writeLines("No data found", "read_counts_summary.txt")
        writeLines("No data found", "model_summary.txt")
        
        png("qc_summary_plot.png", width = 800, height = 600)
        plot(1, 1, type = "n", main = "No Data Found", xlab = "", ylab = "")
        text(1, 1, "No data available for plotting", cex = 1.5)
        dev.off()
        
        quit(save = "no", status = 0)
    }
    
    # Save read counts summary
    data %>%
        dplyr::arrange(sample_id, stage) %>%
        write_delim("read_counts_summary.txt", delim = "\\t")
    
    #### Model Fitting ####
    cat("Fitting statistical model...\\n")
    
    tryCatch({
        model <- glmmTMB(n_reads ~ stage + (1 | sample_id),
                        data = data,
                        family = nbinom2(link = "log"),
                        control = glmmTMBControl(optimizer = "optim",
                                               optArgs = list(method = "BFGS")))
        
        cat("Model fitted successfully\\n")
        
        # Save model summary
        model_summary <- capture.output(summary(model))
        writeLines(model_summary, "model_summary.txt")
        
        #### Plot Generation ####
        cat("Generating summary plot...\\n")
        
        # Generate emmeans and handle grouping
        em_results <- emmeans(model, ~ stage, type = "response")
        
        if (use_multcomp) {
            plot_data <- em_results %>%
                cld(Letters = LETTERS, alpha = 0.05) %>%
                broom::tidy(conf.int = TRUE) %>%
                dplyr::mutate(.group = str_trim(.group),
                       stage = factor(stage, levels = levels(data[["stage"]])))
        } else {
            plot_data <- em_results %>%
                broom::tidy(conf.int = TRUE) %>%
                dplyr::mutate(.group = "A",
                       stage = factor(stage, levels = levels(data[["stage"]])))
        }
        
        # Create the plot
        plot <- plot_data %>%
            ggplot(aes(x = stage, y = response)) +
            geom_line(data = data,
                     aes(y = n_reads, colour = sample_id, group = sample_id),
                     position = position_dodge(0.25),
                     linewidth = 0.1, show.legend = FALSE) +
            geom_point(data = data,
                      aes(y = n_reads, colour = sample_id, group = sample_id),
                      position = position_dodge(0.25),
                      size = 2) +
            geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                         width = 0, linetype = 'dashed',
                         position = position_dodge(0.5),
                         show.legend = FALSE) +
            geom_errorbar(aes(ymin = response - std.error, 
                             ymax = response + std.error),
                         width = 0.25,
                         position = position_dodge(0.5),
                         show.legend = FALSE) +
            geom_line(aes(group = 1),
                     position = position_dodge(0.5),
                     show.legend = FALSE) +
            geom_point(size = 5, position = position_dodge(0.5)) +
            scale_y_continuous(labels = scales::trans_format("log10",
                                                           scales::math_format(10^.x)),
                             trans = scales::log10_trans(),
                             guide = "axis_logticks") +
            guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
            labs(x = 'Processing Stage', 
                 y = 'Number of Reads',
                 colour = 'Sample',
                 title = 'QC Pipeline Read Retention',
                 subtitle = 'Statistical analysis of read counts across processing stages') +
            theme_classic() +
            theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12))
        
        # Add significance annotations if multcomp is available
        if (use_multcomp) {
            sig_data <- plot_data %>%
                dplyr::filter(n_distinct(.group) > 1, .by = stage) %>%
                dplyr::distinct(stage) %>%
                dplyr::mutate(sig = '*')
            
            if (nrow(sig_data) > 0) {
                plot <- plot +
                    geom_text(data = sig_data,
                             aes(x = stage, y = Inf, label = sig),
                             inherit.aes = FALSE, size = 16, vjust = 1)
            }
        }
        
        ggsave('qc_summary_plot.png',
               plot = plot,
               width = 10, height = 8, dpi = 300)
        
        cat("Plot saved successfully\\n")
        
    }, error = function(e) {
        cat("Error in model fitting or plotting:\\n")
        cat(conditionMessage(e), "\\n")
        
        # Create a simple fallback plot if we have data
        if(nrow(data) > 0) {
            simple_plot <- data %>%
                ggplot(aes(x = stage, y = n_reads, colour = sample_id, group = sample_id)) +
                geom_line(linewidth = 1) +
                geom_point(size = 3) +
                scale_y_continuous(labels = scales::trans_format("log10",
                                                               scales::math_format(10^.x)),
                                 trans = scales::log10_trans(),
                                 guide = "axis_logticks") +
                labs(x = 'Processing Stage', 
                     y = 'Number of Reads',
                     colour = 'Sample',
                     title = 'QC Pipeline Read Retention (Simple Plot)',
                     subtitle = 'Read counts across processing stages') +
                theme_classic()
            
            ggsave('qc_summary_plot.png',
                   plot = simple_plot,
                   width = 10, height = 8, dpi = 300)
            
            cat("Fallback plot created\\n")
        } else {
            # Create empty plot
            png("qc_summary_plot.png", width = 800, height = 600)
            plot(1, 1, type = "n", main = "Error in Analysis", xlab = "", ylab = "")
            text(1, 1, "Error occurred during analysis", cex = 1.5)
            dev.off()
        }
        
        # Write error details to model summary
        writeLines(c("Model fitting failed:", conditionMessage(e)), "model_summary.txt")
    })
    
    cat("\\nAnalysis completed at:", format(Sys.time()), "\\n")
    sink()
    """
}