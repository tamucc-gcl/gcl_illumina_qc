#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(purrr)
  library(stringr)
  library(rlang)
  library(patchwork)
})

# Read the combined BLAST results from the pipeline
blast_results <- read_delim('blast_results.tsv', 
                            delim = "\t",
                            show_col_types = FALSE)

# Check if data is empty
if (nrow(blast_results) == 0) {
  cat("No BLAST results to analyze\n")
  # Create empty outputs
  write_csv(data.frame(), 'top_blast_hits.csv')
  dir.create('blast_posteriors', showWarnings = FALSE)
  write_csv(data.frame(), 'blast_posteriors/empty.csv')
  # Create empty plots
  png('blast_raw_pie.png')
  plot.new()
  text(0.5, 0.5, "No BLAST results", cex = 2)
  dev.off()
  png('blast_summary_pie.png')
  plot.new()
  text(0.5, 0.5, "No BLAST results", cex = 2)
  dev.off()
  quit(save = "no", status = 0)
}

taxonomy_ranks <- c("kingdom","phylum","class","order","family","genus","species")

#### Functions ####
safe_neglog10 <- function(x) {
  -log10(pmax(as.numeric(x), 1e-300))
}

# Score hits based on weighted quality metrics
# increasing lambda "sharpens" the softmax - basically make winner take more of the relative score
score_hits <- function(df,
                       w_bitscore = 1.0,
                       w_pident   = 0.25,
                       w_length   = 0.25,
                       w_evalue   = 0.25,
                       lambda     = 6) {
  # Per qseqid standardize the features, compute composite z, softmax -> omega
  df %>%
    mutate(.by = c(sample_id, qseqid),
           
           # zscore features
           z_bits = ifelse(sd(bitscore) > 0, (bitscore - mean(bitscore)) / sd(bitscore), 0),
           z_pid  = ifelse(sd(pident)   > 0, (pident  - mean(pident))  / sd(pident), 0),
           z_len  = ifelse(sd(length)   > 0, (length  - mean(length))  / sd(length), 0),
           z_eval = ifelse(sd(safe_neglog10(evalue)) > 0,
                           (safe_neglog10(evalue) - mean(safe_neglog10(evalue))) / sd(safe_neglog10(evalue)), 0)) %>%
    mutate(.by = c(sample_id, qseqid),
           .keep = 'unused',
           # composite score creation with weights
           z_comp = w_bitscore * z_bits + w_pident * z_pid + w_length * z_len + w_evalue * z_eval,
           
           # softmax within qseqid
           sm_exp = exp(lambda * (z_comp - max(z_comp, na.rm = TRUE))),
           omega_q = sm_exp / sum(sm_exp, na.rm = TRUE)) %>%
    select(-z_comp, -sm_exp)
}


# Get posterior probability of hits based on scores within each sample and if desired a prior based on the number of entries for that taxa existing in the database
# use eta to either eta < 0 downweight taxa overrepresented in the database of eta > 0 favor taxa more common in database
# prior alpha gives uniform pseudo-count to all identified taxa - generally ignore/keep at default. Small values give minor smoothing to "help" rare taxa. Large values make the database counts largely irrelevant
aggregate_rank <- function(df_scored, rank_col, priors_tbl = NULL, prior_eta = 0, prior_alpha0 = 0.5) {
  rank_sym <- sym(rank_col)
  
  # Sum weights per taxon at rank; also count hits and queries supporting
  agg <- df_scored %>%
    filter(!is.na(!!rank_sym), !!rank_sym != "") %>%
    summarize(.by = c(sample_id, !!rank_sym),
              score_sum   = sum(omega_q, na.rm = TRUE),
              n_hits      = n(),
              n_qseqids   = n_distinct(qseqid)) %>%
    rename(taxon = !!rank_sym)
  
  # Join priors if provided (columns expected: rank, taxon, N)
  if (!is.null(priors_tbl)) {
    agg <- agg %>%
      left_join(
        priors_tbl %>% filter(rank == rank_col) %>% select(taxon, N),
        by = "taxon"
      ) %>%
      mutate(N = coalesce(N, 0)) %>%
      mutate(prior = (N + prior_alpha0)^prior_eta)
  } else {
    agg <- agg %>% mutate(prior = 1)
  }
  
  # Posterior ∝ score_sum * prior, normalized per sample
  agg %>%
    mutate(.by = sample_id,
           posterior_raw = score_sum * prior,
           posterior = posterior_raw / sum(posterior_raw, na.rm = TRUE)) %>%
    arrange(sample_id, desc(posterior)) %>%
    select(-prior, -posterior_raw)
}

#Get best taxonomy for hits respecting higher level taxonomy
choose_hierarchical <- function(all_rank_tables, 
                                ranks = c("kingdom","phylum","class","order","family","genus","species"),
                                lineage_cols = c("kingdom","phylum","class","order","family","genus","species"),
                                original_hits) {
  # We ensure each lower rank chosen taxon occurs in at least one hit whose higher-rank lineage matches chosen ancestors.
  # Greedy top-down: at each rank, filter to taxa consistent with already chosen path; pick max posterior.
  chosen <- list()
  ancestors <- list()
  
  for (rk in ranks) {
    tbl <- all_rank_tables[[rk]]
    
    if (length(ancestors) > 0) {
      # Build a filter mask by checking presence in original hits with matching ancestry
      mask <- original_hits
      for (anc_rk in names(ancestors)) {
        mask <- mask %>% filter(.data[[anc_rk]] == ancestors[[anc_rk]])
      }
      valid_taxa <- unique(mask[[rk]])
      
      tbl <- tbl %>% filter(taxon %in% valid_taxa)
      # If filtering removed everything (rare), fall back to unfiltered table for that rank
      if (nrow(tbl) == 0) tbl <- all_rank_tables[[rk]]
    }
    
    top <- tbl %>%
      group_by(sample_id) %>%
      slice_max(order_by = posterior, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(sample_id, taxon, posterior, n_hits, n_qseqids)
    
    chosen[[rk]] <- top
    
    # update ancestors per sample
    ancestors[[rk]] <- NULL  # placeholder for loop scope
    # We store as a named list keyed by rank with scalar per sample below
  }
  
  # Stitch chosen per-rank into a single summary per sample
  summary <- reduce(
    lapply(names(chosen), function(rk) {
      df <- chosen[[rk]] %>%
        rename(!!rk := taxon, !!paste0(rk, "_likelihood") := posterior) %>%
        select(sample_id, all_of(rk), all_of(paste0(rk, "_likelihood")))
      df
    }),
    ~ left_join(.x, .y, by = "sample_id")
  )
  
  # Also provide the per-rank full posterior tables
  list(summary = summary, posteriors = all_rank_tables)
}

# Classify taxonomy
classify_taxa <- function(blast_df,
                          ranks = c("kingdom","phylum","class","order","family","genus","species"),
                          priors_tbl = NULL,
                          prior_eta = 0,
                          prior_alpha0 = 0.5,
                          w_bitscore = 1.0,
                          w_pident   = 0.25,
                          w_length   = 0.25,
                          w_evalue   = 0.25,
                          lambda     = 6) {
  
  stopifnot(all(ranks %in% colnames(blast_df)))
  
  # 1) Score per-query hits
  scored <- score_hits(
    blast_df,
    w_bitscore = w_bitscore,
    w_pident   = w_pident,
    w_length   = w_length,
    w_evalue   = w_evalue,
    lambda     = lambda
  )
  
  # 2) Per-rank aggregation and posterior
  rank_tables <- map(
    ranks,
    ~ aggregate_rank(scored, .x, priors_tbl = priors_tbl, prior_eta = prior_eta, prior_alpha0 = prior_alpha0)
  )
  names(rank_tables) <- ranks
  
  # 3) Choose a hierarchy-consistent path, return both summary and full posteriors
  choose_hierarchical(rank_tables, ranks = ranks, original_hits = blast_df)
}

#### Classify Taxa ####
cat("Classifying taxa from BLAST results...\n")
blast_classification <- classify_taxa(blast_results, prior_eta = -0.25)

#### Plots ####
cat("Generating pie charts...\n")

make_pie_chart <- function(df, max_cat = 5){
  plot_data <- df %>%
    mutate(taxa = fct_lump(taxa, n = max_cat, w = n)) %>%
    summarise(.by = taxa,
              n = sum(n)) %>% 
    mutate(taxa = fct_reorder(taxa, n, .desc = TRUE)) %>%
    mutate(taxa = if ("Other" %in% levels(taxa)) fct_relevel(taxa, "Other", after = Inf) else taxa) %>%
    arrange(taxa) %>%
    mutate(frac = n / sum(n),
           label = paste0(taxa, "\n", scales::percent(frac, 0.1)))
  
  label_cutoff <- 0.05
  
  plot_data %>%
    ggplot(aes(x = "", y = n, fill = taxa)) +
    geom_col(width = 1, color = "white",
             position = position_stack()) +
    coord_polar(theta = "y") +
    geom_text(aes(label = if_else(frac > label_cutoff, taxa, '')),
              position = position_stack(vjust = 0.5),
              size = 3, show.legend = FALSE) +
    theme_void() +
    theme(legend.position = "right")
}

# Raw pie species pie chart
raw_piecharts <- blast_results %>%
  select(sample_id, all_of(taxonomy_ranks)) %>%
  pivot_longer(cols = -sample_id,
               names_to = 'level',
               values_to = 'taxa') %>%
  filter(!is.na(taxa)) %>%
  summarise(.by = c(level, taxa),
            n = n()) %>%
  nest(data = -c(level)) %>%
  rowwise() %>% 
  mutate(plot = list(make_pie_chart(data) +
                       labs(fill = str_to_title(level)))) %>%
  pull(plot) %>%
  wrap_plots()

ggsave('blast_raw_pie.png',
       plot = raw_piecharts, 
       height = 10, width = 10)
cat("  Saved: blast_raw_pie.png\n")

# Summarized species pie chart
summary_piecharts <- blast_classification$summary %>%
  select(sample_id, all_of(taxonomy_ranks)) %>%
  pivot_longer(cols = -sample_id,
               names_to = 'level',
               values_to = 'taxa') %>%
  filter(!is.na(taxa)) %>%
  summarise(.by = c(level, taxa),
            n = n()) %>%
  nest(data = -c(level)) %>%
  rowwise() %>% 
  mutate(plot = list(make_pie_chart(data) +
                       labs(fill = str_to_title(level)))) %>%
  pull(plot) %>%
  wrap_plots()

ggsave('blast_summary_pie.png',
       plot = summary_piecharts, 
       height = 10, width = 10)
cat("  Saved: blast_summary_pie.png\n")

# Write results
write_csv(blast_classification$summary,
          'top_blast_hits.csv')
cat("  Saved: top_blast_hits.csv\n")

dir.create('blast_posteriors', showWarnings = FALSE)
walk2(blast_classification$posteriors,
      names(blast_classification$posteriors),
      ~write_csv(.x, 
                 file = str_c('blast_posteriors/blast_posteriors_', .y, '.csv')))
cat("  Saved posterior files to blast_posteriors/\n")

cat("\nSpecies identification analysis complete!\n")