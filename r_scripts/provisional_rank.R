#!/usr/bin/env Rscript
# provisional_rank.R  (CHUNK 3b)
# Weight-free rank aggregation of the CHEAP signals over all candidates, to
# produce a provisional ranking. The top candidates by this ranking are handed to
# the expensive stage-2 (bcftools) step (gap-based top-N selection = chunk 4; for
# now we emit the full ranked table and a simple top-N survivors list).
#
# Signals (each becomes a per-candidate RANK; lower rank = better):
#   1a inflection : across-grid contig-count vs total-length elbow. We rank by
#        distance from the inflection point of the n_contigs-vs-cutoff curve,
#        approximated as the point of maximum curvature on sorted n_contigs.
#   1b redundancy : rank by redundancy fraction (lower = better).
#   2 anchor      : if expected_loci supplied (>0), rank by |n_contigs - expected|.
#                   If not supplied, this signal is DROPPED from the aggregation.
#   5b NB cutoff1 : rank by |c1 - nb_cutoff1| (proximity on the c1 axis).
#
# Aggregation: MEAN of available signal ranks (weight-free). Ties broken by 1b
# redundancy then by id for determinism. This is intentionally NOT a weighted
# composite — the whole point of the redesign is to avoid free weights.
#
# Usage:
#   Rscript provisional_rank.R cheap_all.tsv nb_cutoff1.value EXPECTED_LOCI MIN_SURV
#     cheap_all.tsv  : concatenated per-candidate rows (no header):
#                      id c1 c2 sim n_contigs total_len mean_len redundancy
#     nb_cutoff1     : file with the NB-derived cutoff1 integer
#     EXPECTED_LOCI  : integer expected RAD loci, or "NA"/0 to disable signal 2
#     MIN_SURV       : how many survivors to emit (stage-2 input; chunk 4 will
#                      replace this fixed N with the gap-based rule)
#
# Outputs:
#   provisional_rank.tsv  : full table with per-signal ranks + aggregate
#   survivors.txt         : one candidate id per line (top MIN_SURV)

suppressPackageStartupMessages({ library(readr); library(dplyr) })

args <- commandArgs(trailingOnly = TRUE)
cheap_file <- args[1]
nb_file    <- args[2]
expected   <- suppressWarnings(as.numeric(args[3]))
min_surv   <- suppressWarnings(as.integer(args[4])); if (is.na(min_surv)) min_surv <- 3L

nb_cut <- tryCatch(as.numeric(readLines(nb_file)[1]), error = function(e) NA_real_)

d <- read_tsv(cheap_file,
              col_names = c("id","c1","c2","sim","n_contigs","total_len","mean_len","redundancy"),
              col_types = "ciidiiid")

if (nrow(d) == 0) {
  writeLines("id\tc1\tc2\tsim\tn_contigs\tredundancy\trank_inflect\trank_redun\trank_anchor\trank_nb\tagg_rank",
             "provisional_rank.tsv")
  writeLines(character(0), "survivors.txt")
  quit(save = "no", status = 0)
}

# ---- signal 1a: inflection of n_contigs across the grid ----
# Sort candidates by n_contigs; the inflection (elbow) is the point of maximum
# second-difference (curvature) on the sorted curve. Rank by distance from it.
ord <- order(d$n_contigs)
nc_sorted <- d$n_contigs[ord]
if (length(nc_sorted) >= 3) {
  d2 <- c(0, diff(diff(nc_sorted)), 0)            # discrete curvature
  elbow_idx <- which.max(abs(d2))
  elbow_val <- nc_sorted[elbow_idx]
} else {
  elbow_val <- median(d$n_contigs)
}
d$dist_inflect <- abs(d$n_contigs - elbow_val)
d$rank_inflect <- rank(d$dist_inflect, ties.method = "average")

# ---- signal 1b: redundancy (lower better) ----
d$rank_redun <- rank(d$redundancy, ties.method = "average")

# ---- signal 2: anchor proximity (optional) ----
use_anchor <- is.finite(expected) && expected > 0
if (use_anchor) {
  d$dist_anchor <- abs(d$n_contigs - expected)
  d$rank_anchor <- rank(d$dist_anchor, ties.method = "average")
} else {
  d$dist_anchor <- NA_real_
  d$rank_anchor <- NA_real_
}

# ---- signal 5b: NB cutoff1 proximity on the c1 axis ----
if (is.finite(nb_cut)) {
  d$dist_nb <- abs(d$c1 - nb_cut)
  d$rank_nb <- rank(d$dist_nb, ties.method = "average")
} else {
  d$dist_nb <- NA_real_
  d$rank_nb <- NA_real_
}

# ---- weight-free aggregation: mean of AVAILABLE signal ranks ----
rank_cols <- c("rank_inflect","rank_redun","rank_anchor","rank_nb")
d$agg_rank <- rowMeans(d[, rank_cols], na.rm = TRUE)

ranked <- d %>%
  arrange(agg_rank, redundancy, id) %>%
  mutate(provisional_order = row_number())

write_tsv(
  ranked %>% select(id, c1, c2, sim, n_contigs, redundancy,
                    rank_inflect, rank_redun, rank_anchor, rank_nb,
                    agg_rank, provisional_order),
  "provisional_rank.tsv"
)

# survivors = top min_surv ids (chunk 4 replaces with gap-based N)
n_take <- max(1L, min(min_surv, nrow(ranked)))
writeLines(ranked$id[seq_len(n_take)], "survivors.txt")

cat(sprintf("provisional_rank: %d candidates, anchor=%s, nb_cutoff1=%s, survivors=%d\n",
            nrow(d), ifelse(use_anchor, as.character(expected), "disabled"),
            ifelse(is.finite(nb_cut), as.character(nb_cut), "NA"), n_take))
cat("Top survivors:\n"); cat(paste0("  ", ranked$id[seq_len(n_take)]), sep = "\n"); cat("\n")
