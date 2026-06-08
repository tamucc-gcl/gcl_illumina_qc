#!/usr/bin/env Rscript
# score_assemblies.R
# Rank de novo cluster_similarity sweep candidates and pick the best.
#
# Usage: Rscript score_assemblies.R sweep_scores.tsv
# Input columns: sim  n_samples  map_rate  pp_rate  softclip_per_read  mean_AS
#
# Composite score (higher = better):
#   z(map_rate)*w_map + z(pp_rate)*w_pp - z(softclip_per_read)*w_sc + z(mean_AS)*w_as
# Properly-paired rate and mapping rate are rewarded; soft-clipping is penalized.
# Alignment score is included with a small weight. Tune the weights below.
#
# Outputs: best_sim.value, sweep_summary.tsv, sweep_comparison.png

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

# ---- weights (tune here) ----
w_map <- 1.0   # mapping rate
w_pp  <- 1.5   # properly-paired rate (primary signal of a good reference)
w_sc  <- 1.0   # soft-clipping penalty
w_as  <- 0.25  # mean alignment score

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]

scores <- read_tsv(infile, show_col_types = FALSE,
                   col_types = cols(sim = col_character(), .default = col_double()))

# z-score helper that is safe when all values are identical (sd == 0)
z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
}

ranked <- scores %>%
  mutate(
    z_map = z(map_rate),
    z_pp  = z(pp_rate),
    z_sc  = z(softclip_per_read),
    z_as  = z(mean_AS),
    composite = w_map * z_map + w_pp * z_pp - w_sc * z_sc + w_as * z_as
  ) %>%
  arrange(desc(composite))

write_tsv(ranked, "sweep_summary.tsv")

best <- ranked$sim[1]
writeLines(as.character(best), "best_sim.value")
cat("Ranked candidates:\n"); print(as.data.frame(ranked), row.names = FALSE)
cat(sprintf("\nSelected cluster_similarity = %s\n", best))

# ---- comparison plot: one panel per metric, best candidate highlighted ----
plot_df <- ranked %>%
  mutate(sim = factor(sim, levels = sim[order(as.numeric(sim))]),
         is_best = sim == best) %>%
  select(sim, is_best, map_rate, pp_rate, softclip_per_read, mean_AS, composite) %>%
  pivot_longer(c(map_rate, pp_rate, softclip_per_read, mean_AS, composite),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         map_rate = "Mapping rate (%)",
                         pp_rate = "Properly paired (%)",
                         softclip_per_read = "Soft-clip / read (lower better)",
                         mean_AS = "Mean alignment score",
                         composite = "Composite (selected = max)"))

p <- ggplot(plot_df, aes(x = sim, y = value, fill = is_best)) +
  geom_col() +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                    guide = "none") +
  labs(title = "Cluster-similarity sweep: assembly quality by candidate",
       subtitle = paste0("Selected cluster_similarity = ", best),
       x = "cluster_similarity", y = NULL) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave("sweep_comparison.png", p, width = 9, height = 6, dpi = 200)
cat("Saved sweep_comparison.png\n")
