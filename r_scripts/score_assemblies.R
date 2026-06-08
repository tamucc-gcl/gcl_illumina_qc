#!/usr/bin/env Rscript
# score_assemblies.R
# Rank de novo sweep candidates over (cutoff1, cutoff2, cluster_similarity) and
# pick the best by a composite of map-back quality.
#
# Usage: Rscript score_assemblies.R sweep_scores.tsv
# Columns: id c1 c2 sim n_samples map_rate pp_rate softclip_per_read mean_AS
#
# Composite (higher = better):
#   z(map_rate)*w_map + z(pp_rate)*w_pp - z(softclip_per_read)*w_sc + z(mean_AS)*w_as
#
# Outputs: best_id.value, sweep_summary.tsv, sweep_comparison.png

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

# ---- weights (tune here) ----
w_map <- 1.0
w_pp  <- 1.5   # properly-paired rate (primary signal of a good reference)
w_sc  <- 1.0   # soft-clipping penalty
w_as  <- 0.25

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]

scores <- read_tsv(infile, show_col_types = FALSE,
                   col_types = cols(id = col_character(), sim = col_character(),
                                    .default = col_double()))

z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
}

ranked <- scores %>%
  mutate(z_map = z(map_rate), z_pp = z(pp_rate),
         z_sc = z(softclip_per_read), z_as = z(mean_AS),
         composite = w_map*z_map + w_pp*z_pp - w_sc*z_sc + w_as*z_as) %>%
  arrange(desc(composite))

write_tsv(ranked, "sweep_summary.tsv")

best <- ranked$id[1]
writeLines(as.character(best), "best_id.value")
cat("Top candidates:\n"); print(head(as.data.frame(ranked), 10), row.names = FALSE)
cat(sprintf("\nSelected candidate = %s\n", best))

# ---- heatmap: composite over (cutoff1 x cutoff2), faceted by cluster_similarity ----
plot_df <- ranked %>%
  mutate(c1 = factor(c1, levels = sort(unique(c1))),
         c2 = factor(c2, levels = sort(unique(c2))),
         sim_lab = paste0("sim = ", sim),
         is_best = id == best)

p <- ggplot(plot_df, aes(x = c1, y = c2, fill = composite)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_tile(data = subset(plot_df, is_best),
            color = "black", linewidth = 1.1, fill = NA) +
  geom_text(aes(label = sprintf("%.1f", map_rate)), size = 2.6) +
  facet_wrap(~ sim_lab) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", name = "composite") +
  labs(title = "De novo sweep: composite score by cutoffs and cluster_similarity",
       subtitle = paste0("Best = ", best, " (tile labels = mapping rate %); black border = selected"),
       x = "cutoff1 (min reads / individual)",
       y = "cutoff2 (min individuals)") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9))

ggsave("sweep_comparison.png", p, width = 10, height = 7, dpi = 200)
cat("Saved sweep_comparison.png\n")
