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
    
    # List all available files for debugging
    cat("Available files:\\n")
    available_files <- list.files(".", pattern = "\\.txt\\$", full.names = TRUE)
    cat(paste(available_files, collapse = "\\n"), "\\n\\n")
    
    # Process the data
    data <- available_files %>%
        tibble(file = .) %>%
        mutate(stage = case_when(
            str_detect(file, "raw.*fastqc") ~ "raw",
            str_detect(file, "trim.*3") ~ "trim3", 
            str_detect(file, "clumpify") ~ "dedup",
            str_detect(file, "trim.*5") ~ "trim5",
            str_detect(file, "fastq.*screen") ~ "fqscrn",
            str_detect(file, "repair") ~ "repr",
            str_detect(file, "mapping_summary") ~ "map",
            TRUE ~ "unknown"
        )) %>%
        filter(stage != "unknown") %>%
        mutate(stage = factor(stage, 
                            levels = c("raw", "trim3", "dedup", "trim5", 
                                     "fqscrn", "repr", "map"), 
                            ordered = TRUE)) %>%
        rowwise(stage) %>%
        reframe({
            if(stage == "map") {
                # For mapping summary, read different column
                read_delim(file, delim = "\\t", 
                         show_col_types = FALSE, 
                         na = c("", "NA", "N/A")) %>%
                select(sample_id = Sample,
                       n_reads = Mapped_Paired)
            } else {
                # For MultiQC files, read fastqc-total_sequences
                read_delim(file, delim = "\\t", 
                         show_col_types = FALSE, 
                         na = c("", "NA", "N/A")) %>%
                select(sample_id = Sample,
                       n_reads = matches('fastqc-total_sequences'))
            }
        }) %>%
        # Filter for R1 files only (except mapping) and clean sample names
        filter(stage == 'map' | str_detect(sample_id, '\\\\.(r)?1\\$')) %>%
        mutate(n_reads = if_else(stage != 'map', 1e6 * n_reads, n_reads) %>%
                       round(0) %>%
                       as.integer(),
               sample_id = str_remove(sample_id, '_fp1.*r1\\$') %>%
                          str_remove('\\\\.1\\$'))
    
    cat("Data processing completed\\n")
    cat("Samples found:", paste(unique(data\$sample_id), collapse = ", "), "\\n")
    cat("Stages found:", paste(levels(data\$stage), collapse = ", "), "\\n\\n")
    
    # Save read counts summary
    data %>%
        arrange(sample_id, stage) %>%
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
                tidy(conf.int = TRUE) %>%
                mutate(.group = str_trim(.group),
                       stage = factor(stage, levels = levels(data\$stage)))
        } else {
            plot_data <- em_results %>%
                tidy(conf.int = TRUE) %>%
                mutate(.group = "A",  # Default group when multcomp not available
                       stage = factor(stage, levels = levels(data\$stage)))
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
            scale_shape_discrete(labels = str_to_sentence) +
            guides(shape = guide_legend(override.aes = list(size = 5, alpha = 1)),
                   colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
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
                filter(n_distinct(.group) > 1, .by = stage) %>%
                distinct(stage) %>%
                mutate(sig = '*')
            
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
        
        # Create a simple fallback plot
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
        
        # Write error details to model summary
        writeLines(c("Model fitting failed:", conditionMessage(e)), "model_summary.txt")
    })
    
    cat("\\nAnalysis completed at:", format(Sys.time()), "\\n")
    sink()
    """
}