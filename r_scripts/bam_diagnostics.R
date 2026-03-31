library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(forcats)

# Helper: read one of the count-value files produced by samtools_stats
# Format: "  <count> <value>" (output of `sort -n | uniq -c`)
read_count_value <- function(file, value_col) {
  read_table(file,
             col_names = c("count", value_col),
             col_types = "ii",
             comment = "") |>
    filter(!is.na(count), !is.na(.data[[value_col]]))
}

# Expand count-value pairs into one-row-per-read (via uncount),
# then attach sample_id
load_files <- function(pattern, value_col, max_val = Inf) {
  files <- list.files(pattern = pattern, full.names = FALSE)
  if (length(files) == 0L) return(tibble())

  lapply(files, function(f) {
    sample <- str_remove(f, paste0("_", pattern |> str_remove("\\*") |> str_remove("\\..*$"), ".*$"))
    # Simpler: strip the fixed suffix
    sample <- str_remove(f, str_c("_", value_col, "_stats\\.txt$"))
    read_count_value(f, value_col) |>
      filter(.data[[value_col]] <= max_val) |>
      tidyr::uncount(count) |>
      mutate(sample_id = sample)
  }) |>
    bind_rows()
}

# ── Soft clipping ──────────────────────────────────────────────────────────────
soft_clip_data <- list.files(pattern = "soft_clipping_stats\\.txt$") |>
  tibble(file = _) |>
  mutate(sample_id = str_remove(file, "_soft_clipping_stats\\.txt$")) |>
  rowwise(sample_id) |>
  reframe({
    read_table(file, col_names = c("count", "soft_clip_bases"), col_types = "ii") |>
      filter(!is.na(count), !is.na(soft_clip_bases))
  }) |>
  tidyr::uncount(count) |>
  mutate(sample_id = fct_reorder(sample_id, soft_clip_bases))

if (nrow(soft_clip_data) > 0) {
  n_samples <- n_distinct(soft_clip_data$sample_id)

  calc_height <- function(n, base = 1.5, per_sample = 0.35, max_h = 20)
    min(base + n * per_sample, max_h)

  p_sc <- ggplot(soft_clip_data, aes(y = sample_id, x = soft_clip_bases)) +
    geom_violin(fill = "steelblue", alpha = 0.7) +
    scale_x_continuous(labels = scales::comma_format()) +
    labs(x = "Soft-clipped bases per read",
         y = NULL,
         title = "Soft Clipping Distribution") +
    theme_classic() +
    theme(panel.background = element_rect(colour = "black"))

  ggsave("soft_clipping_violin.png", plot = p_sc,
         height = calc_height(n_samples), width = 5)
  cat("Saved soft_clipping_violin.png\n")
} else {
  cat("No soft clipping data found — writing placeholder\n")
  png("soft_clipping_violin.png", width = 400, height = 300)
  plot.new()
  text(0.5, 0.5, "No soft clipping data", cex = 1.5, col = "grey50")
  dev.off()
}

# ── Alignment scores ───────────────────────────────────────────────────────────
aln_score_data <- list.files(pattern = "alignment_score_stats\\.txt$") |>
  tibble(file = _) |>
  mutate(sample_id = str_remove(file, "_alignment_score_stats\\.txt$")) |>
  rowwise(sample_id) |>
  reframe({
    read_table(file, col_names = c("count", "alignment_score"), col_types = "ii") |>
      filter(!is.na(count), !is.na(alignment_score))
  }) |>
  tidyr::uncount(count) |>
  mutate(sample_id = fct_reorder(sample_id, alignment_score))

if (nrow(aln_score_data) > 0) {
  n_samples <- n_distinct(aln_score_data$sample_id)

  p_as <- ggplot(aln_score_data, aes(y = sample_id, x = alignment_score)) +
    geom_violin(fill = "darkorange", alpha = 0.7) +
    scale_x_continuous(labels = scales::comma_format()) +
    labs(x = "Alignment score (AS tag)",
         y = NULL,
         title = "Alignment Score Distribution") +
    theme_classic() +
    theme(panel.background = element_rect(colour = "black"))

  ggsave("alignment_score_violin.png", plot = p_as,
         height = calc_height(n_samples), width = 5)
  cat("Saved alignment_score_violin.png\n")
} else {
  cat("No alignment score data found — writing placeholder\n")
  png("alignment_score_violin.png", width = 400, height = 300)
  plot.new()
  text(0.5, 0.5, "No alignment score data", cex = 1.5, col = "grey50")
  dev.off()
}