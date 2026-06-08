#!/usr/bin/env Rscript
# assembly_diagnostics.R
# Auto-detect a de novo assembly cutoff (cutoff1 or cutoff2) from a frequency
# table and plot the "sequences retained vs. threshold" curve.
#
# Usage:
#   Rscript assembly_diagnostics.R cutoff1 coverage_freq.txt
#   Rscript assembly_diagnostics.R cutoff2 indiv_freq.txt
#
# Input freq file: two whitespace-separated columns produced by `... | sort -n | uniq -c`
#   col1 = number of records (how many unique sequences had that value)
#   col2 = value (within-individual coverage for cutoff1; #individuals for cutoff2)
#
# Outputs:
#   <mode>.value      single integer (the auto-selected cutoff)
#   <mode>_curve.png  the diagnostic curve with the chosen cutoff marked

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

args      <- commandArgs(trailingOnly = TRUE)
mode      <- args[1]            # "cutoff1" or "cutoff2"
freq_file <- args[2]

stopifnot(mode %in% c("cutoff1", "cutoff2"))

value_file <- paste0(mode, ".value")
plot_file  <- paste0(mode, "_curve.png")

write_default <- function(val = 4L, msg = "No data") {
  writeLines(as.character(as.integer(val)), value_file)
  png(plot_file, width = 800, height = 600)
  plot.new(); text(0.5, 0.5, msg, cex = 1.5, col = "grey40")
  dev.off()
  cat(sprintf("[%s] %s -> default = %d\n", mode, msg, as.integer(val)))
  quit(save = "no", status = 0)
}

# ---- Read frequency table -----------------------------------------------------
freq <- tryCatch(
  read_table(freq_file, col_names = c("records", "value"), col_types = "ii"),
  error = function(e) NULL
)
if (is.null(freq) || nrow(freq) == 0) write_default(msg = "No data for cutoff curve")
freq <- freq %>% filter(!is.na(records), !is.na(value), records > 0)
if (nrow(freq) == 0) write_default(msg = "No data for cutoff curve")

# ---- Build "number of sequences with value >= x" curve ------------------------
# cutoff1 (coverage) is capped at 40; cutoff2 (n individuals) runs to its max.
max_observed <- max(freq$value)
max_x <- if (mode == "cutoff1") min(max_observed, 40L) else max_observed
xs <- seq.int(2L, max(2L, max_x))

curve <- tibble(
  x = xs,
  n = vapply(xs, function(t) sum(freq$records[freq$value >= t]), numeric(1))
) %>%
  filter(n > 0)

if (nrow(curve) < 3) write_default(val = max(2L, min(curve$x)),
                                   msg = "Too few points for knee detection")

# ---- Kneedle-style knee detection ---------------------------------------------
# Curves are monotonically decreasing and convex (steep drop -> flat tail).
# Normalize to [0,1], draw the chord between endpoints; the knee is the point of
# maximum vertical gap between the chord and the curve.
# change to use log10 to scale y better for post-clumpify results
find_knee <- function(x, y) {
  y <- log10(pmax(y, 1))          # ← compress the orders-of-magnitude drop
  xn <- (x - min(x)) / (max(x) - min(x))
  yn <- (y - min(y)) / (max(y) - min(y))
  chord <- yn[1] + (yn[length(yn)] - yn[1]) * (xn - xn[1]) / (xn[length(xn)] - xn[1])
  x[which.max(chord - yn)]
}

knee <- as.integer(round(find_knee(curve$x, curve$n)))
# cap kept low for deduplicated (Clumpify'd) input so cutoffs don't over-prune
if (mode == "cutoff1") knee <- max(2L, min(knee, 5L))
if (mode == "cutoff2") knee <- max(2L, min(knee, 4L))
writeLines(as.character(knee), value_file)

# ---- Plot ---------------------------------------------------------------------
title <- if (mode == "cutoff1")
  "Cutoff 1: unique sequences vs. within-individual coverage" else
  "Cutoff 2: unique sequences vs. number of individuals"
xlab  <- if (mode == "cutoff1")
  "Minimum within-individual coverage (reads per individual)" else
  "Minimum number of individuals"

p <- ggplot(curve, aes(x, n)) +
  geom_line(linewidth = 0.8, color = "grey30") +
  geom_point(size = 1.6, color = "grey30") +
  geom_vline(xintercept = knee, color = "red", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = knee, y = max(curve$n),
           label = paste0("  auto = ", knee), color = "red", hjust = 0, vjust = 1) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = title, x = xlab, y = "Unique sequences retained") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave(plot_file, plot = p, width = 7, height = 5, dpi = 200)
cat(sprintf("[%s] auto-selected = %d (curve range %d-%d)\n",
            mode, knee, min(curve$x), max(curve$x)))
