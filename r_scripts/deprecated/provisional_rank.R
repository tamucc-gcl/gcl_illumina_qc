#!/usr/bin/env Rscript
# provisional_rank.R  (CHUNK 3c)
# Weight-free rank aggregation of the CHEAP signals over all candidates.
#
# RANKING signals (each -> per-candidate rank, lower = better):
#   5b NB cutoff1 : rank by |c1 - nb_cutoff1| (proximity on the c1 axis)
#   2  cov CV     : rank by coverage-uniformity CV (lower = more uniform = better)
#   2  cov Gini   : rank by coverage-uniformity Gini (lower = more uniform = better)
#   2  anchor     : if expected_loci > 0, rank by |n_contigs - expected| (else dropped)
# REPORT-ONLY (computed, NOT in aggregate):
#   1a inflection : distance from across-grid contig-count elbow (duplicates NB)
#
# Aggregation: MEAN of available signal ranks (weight-free). Ties broken by CV
# then id for determinism.
#
# Usage:
#   Rscript provisional_rank.R cheap_all.tsv cv_all.tsv nb_cutoff1.value EXPECTED_LOCI MIN_SURV
#     cheap_all.tsv : id c1 c2 sim n_contigs total_len mean_len   (no header)
#     cv_all.tsv    : id c1 c2 sim cv gini                        (no header; NA allowed)
#     nb_cutoff1    : file with NB-derived cutoff1 integer
#     EXPECTED_LOCI : integer, or "NA"/0 to disable anchor
#     MIN_SURV      : survivors to emit (chunk 4 replaces with gap-based rule)
#
# Outputs: provisional_rank.tsv, survivors.txt

suppressPackageStartupMessages({ library(readr); library(dplyr) })

args <- commandArgs(trailingOnly = TRUE)
cheap_file <- args[1]
cv_file    <- args[2]
nb_file    <- args[3]
expected   <- suppressWarnings(as.numeric(args[4]))
min_surv   <- suppressWarnings(as.integer(args[5])); if (is.na(min_surv)) min_surv <- 3L

nb_cut <- tryCatch(as.numeric(readLines(nb_file)[1]), error = function(e) NA_real_)

cheap <- read_tsv(cheap_file,
            col_names = c("id","c1","c2","sim","n_contigs","total_len","mean_len"),
            col_types = "ciidiii")

cv <- tryCatch(
  read_tsv(cv_file, col_names = c("id","c1","c2","sim","frac_hi","gini","cv"),
           col_types = "ciidddd", na = c("NA","")),
  error = function(e) NULL
)

if (nrow(cheap) == 0) {
  writeLines("id\tc1\tc2\tsim\tn_contigs\tfrac_hi\tgini\tcv\tdist_inflect\trank_nb\trank_frac_hi\trank_gini\trank_cv\trank_anchor\tagg_rank\tprovisional_order",
             "provisional_rank.tsv")
  writeLines(character(0), "survivors.txt"); quit(save="no", status=0)
}

# join coverage-uniformity stats onto the candidate table by id
if (!is.null(cv) && nrow(cv) > 0) {
  d <- left_join(cheap, cv %>% select(id, frac_hi, gini, cv), by = "id")
} else {
  d <- cheap %>% mutate(frac_hi = NA_real_, gini = NA_real_, cv = NA_real_)
}

# ---- 1a inflection (REPORT-ONLY): distance from across-grid contig-count elbow ----
ord <- order(d$n_contigs); nc_sorted <- d$n_contigs[ord]
if (length(nc_sorted) >= 3) {
  d2 <- c(0, diff(diff(nc_sorted)), 0)
  elbow_val <- nc_sorted[which.max(abs(d2))]
} else elbow_val <- median(d$n_contigs)
d$dist_inflect <- abs(d$n_contigs - elbow_val)   # reported, not ranked

# ---- 5b NB proximity on the c1 axis ----
d$rank_nb <- if (is.finite(nb_cut)) rank(abs(d$c1 - nb_cut), ties.method="average") else NA_real_

# ---- 2 coverage uniformity (lower better): frac_hi PRIMARY, gini + cv secondary ----
# (5a: per-base depth, median-normalized. frac_hi = paralog-pileup fraction.)
d$rank_frac_hi <- if (any(is.finite(d$frac_hi))) rank(d$frac_hi, ties.method="average", na.last="keep") else NA_real_
d$rank_gini    <- if (any(is.finite(d$gini)))    rank(d$gini,    ties.method="average", na.last="keep") else NA_real_
d$rank_cv      <- if (any(is.finite(d$cv)))      rank(d$cv,      ties.method="average", na.last="keep") else NA_real_

# ---- 2 anchor proximity (optional) ----
use_anchor <- is.finite(expected) && expected > 0
d$rank_anchor <- if (use_anchor) rank(abs(d$n_contigs - expected), ties.method="average") else NA_real_

# ---- weight-free aggregation: mean of AVAILABLE ranks (inflection EXCLUDED) ----
# NOTE (5a): emitting frac_hi/gini/cv all as columns to inspect; only frac_hi (the
# size-robust primary) enters the aggregate by default. gini/cv are reported for
# comparison and can be promoted after we see real-data behavior.
rank_cols <- c("rank_nb","rank_frac_hi","rank_anchor")
d$agg_rank <- rowMeans(d[, rank_cols], na.rm = TRUE)

ranked <- d %>%
  arrange(agg_rank, frac_hi, id) %>%
  mutate(provisional_order = row_number())

write_tsv(
  ranked %>% select(id, c1, c2, sim, n_contigs, frac_hi, gini, cv, dist_inflect,
                    rank_nb, rank_frac_hi, rank_gini, rank_cv, rank_anchor,
                    agg_rank, provisional_order),
  "provisional_rank.tsv"
)

n_take <- max(1L, min(min_surv, nrow(ranked)))
writeLines(ranked$id[seq_len(n_take)], "survivors.txt")

cat(sprintf("provisional_rank: %d candidates | nb_cutoff1=%s | anchor=%s | CV signal=%s | survivors=%d\n",
            nrow(d), ifelse(is.finite(nb_cut), as.character(nb_cut), "NA"),
            ifelse(use_anchor, as.character(expected), "disabled"),
            ifelse(any(is.finite(d$frac_hi)), "active", "MISSING"), n_take))
cat("Top survivors:\n"); cat(paste0("  ", ranked$id[seq_len(n_take)]), sep="\n"); cat("\n")
