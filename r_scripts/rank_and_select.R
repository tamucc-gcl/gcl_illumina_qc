#!/usr/bin/env Rscript
# rank_and_select.R  (CHUNK 5b — single 1-pass aggregation, no gate)
# Weight-free rank aggregation over all candidate signals -> winner.
#
# RANKING signals (each -> per-candidate rank; the code orients each so LOWER
# rank = BETTER, then averages available ranks):
#   NB cutoff1      : rank by |c1 - nb_cutoff1|            (proximity, lower better)
#   anchor          : rank by |n_contigs - expected_loci|  (if expected>0; lower better)
#   concordance     : rank by concordance DESC             (higher better -> primary)
#   r80_loci        : rank by r80 DESC                     (higher better)
#   snps_per_locus  : rank by density DESC                 (higher better, secondary)
# REPORT-ONLY (computed, excluded from aggregate):
#   inflection      : distance from across-grid contig-count elbow (duplicates NB)
#
# Aggregation: MEAN of available signal ranks (weight-free). The winner is the
# minimum aggregate rank; ties broken by concordance DESC then id.
#
# Usage:
#   Rscript rank_and_select.R cheap_all.tsv snp_all.tsv nb_cutoff1.value EXPECTED_LOCI
#     cheap_all.tsv : id c1 c2 sim n_contigs total_len mean_len      (no header)
#     snp_all.tsv   : id c1 c2 sim concordance r80_loci snps_per_locus n_snps n_called_contigs
#     nb_cutoff1    : file with NB-derived cutoff1 integer
#     EXPECTED_LOCI : integer expected RAD loci, or "NA"/0 to disable anchor
#
# Outputs: final_rank.tsv, best_id.value, optimize_plot.png

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(tidyr) })

args <- commandArgs(trailingOnly = TRUE)
cheap_file <- args[1]; snp_file <- args[2]; nb_file <- args[3]
expected <- suppressWarnings(as.numeric(args[4]))

nb_cut <- tryCatch(as.numeric(readLines(nb_file)[1]), error = function(e) NA_real_)

cheap <- read_tsv(cheap_file,
            col_names = c("id","c1","c2","sim","n_contigs","total_len","mean_len"),
            col_types = "ciidiii")

snp <- tryCatch(
  read_tsv(snp_file,
    col_names = c("id","c1","c2","sim","concordance","r80_loci","snps_per_locus","n_snps","n_called_contigs"),
    col_types = "ciidddddi", na = c("NA","")),
  error = function(e) NULL
)

if (nrow(cheap) == 0) {
  writeLines("id\tc1\tc2\tsim\tn_contigs\tconcordance\tr80_loci\tsnps_per_locus\tdist_inflect\trank_nb\trank_anchor\trank_conc\trank_r80\trank_snpdens\tagg_rank\tfinal_order", "final_rank.tsv")
  writeLines("NA", "best_id.value"); quit(save="no", status=0)
}

# join SNP signals onto candidate table by id
if (!is.null(snp) && nrow(snp) > 0) {
  d <- left_join(cheap, snp %>% select(id, concordance, r80_loci, snps_per_locus), by = "id")
} else {
  d <- cheap %>% mutate(concordance = NA_real_, r80_loci = NA_real_, snps_per_locus = NA_real_)
}

# ---- inflection (REPORT-ONLY) ----
ord <- order(d$n_contigs); nc_sorted <- d$n_contigs[ord]
if (length(nc_sorted) >= 3) {
  d2 <- c(0, diff(diff(nc_sorted)), 0); elbow_val <- nc_sorted[which.max(abs(d2))]
} else elbow_val <- median(d$n_contigs)
d$dist_inflect <- abs(d$n_contigs - elbow_val)

# ---- ranks (lower rank = better in every case) ----
# proximity signals: rank ascending on |.|
d$rank_nb     <- if (is.finite(nb_cut)) rank(abs(d$c1 - nb_cut), ties.method="average") else NA_real_
use_anchor <- is.finite(expected) && expected > 0
d$rank_anchor <- if (use_anchor) rank(abs(d$n_contigs - expected), ties.method="average") else NA_real_
# quality signals: higher is better -> rank on NEGATED value so lower rank = better
d$rank_conc    <- if (any(is.finite(d$concordance)))    rank(-d$concordance,    ties.method="average", na.last="keep") else NA_real_
d$rank_r80     <- if (any(is.finite(d$r80_loci)))       rank(-d$r80_loci,       ties.method="average", na.last="keep") else NA_real_
d$rank_snpdens <- if (any(is.finite(d$snps_per_locus))) rank(-d$snps_per_locus, ties.method="average", na.last="keep") else NA_real_

# ---- weight-free aggregation: mean of AVAILABLE ranks (inflection excluded) ----
rank_cols <- c("rank_nb","rank_anchor","rank_conc","rank_r80","rank_snpdens")
d$agg_rank <- rowMeans(d[, rank_cols], na.rm = TRUE)

ranked <- d %>%
  arrange(agg_rank, desc(concordance), id) %>%
  mutate(final_order = row_number())

write_tsv(
  ranked %>% select(id, c1, c2, sim, n_contigs, concordance, r80_loci, snps_per_locus,
                    dist_inflect, rank_nb, rank_anchor, rank_conc, rank_r80, rank_snpdens,
                    agg_rank, final_order),
  "final_rank.tsv"
)

best_id <- ranked$id[1]
writeLines(as.character(best_id), "best_id.value")

# ---- plot: signals for the top candidates (faceted, for the report) ----
tryCatch({
  topn <- head(ranked, 15)
  pl <- topn %>%
    select(id, concordance, r80_loci, snps_per_locus, agg_rank) %>%
    pivot_longer(-c(id, agg_rank), names_to = "signal", values_to = "value") %>%
    mutate(id = factor(id, levels = rev(topn$id)))
  p <- ggplot(pl, aes(value, id)) +
    geom_point(size = 2, color = "steelblue") +
    facet_wrap(~ signal, scales = "free_x") +
    labs(title = "De novo optimization: top candidates by signal",
         subtitle = paste0("Winner: ", best_id, " (lowest mean rank). Concordance is primary."),
         x = NULL, y = NULL) +
    theme_bw() + theme(plot.title = element_text(face="bold"))
  ggsave("optimize_plot.png", p, width = 9, height = 6, dpi = 150)
}, error = function(e) {
  png("optimize_plot.png", width = 700, height = 500)
  plot.new(); text(0.5, 0.5, paste("plot failed:", conditionMessage(e)), col="grey40")
  dev.off()
})

cat(sprintf("rank_and_select: %d candidates | nb=%s | anchor=%s | conc=%s r80=%s snpdens=%s | winner=%s\n",
    nrow(d), ifelse(is.finite(nb_cut), as.character(nb_cut), "NA"),
    ifelse(use_anchor, as.character(expected), "disabled"),
    ifelse(any(is.finite(d$concordance)), "active","MISSING"),
    ifelse(any(is.finite(d$r80_loci)), "active","MISSING"),
    ifelse(any(is.finite(d$snps_per_locus)), "active","MISSING"), best_id))
cat("Top 5:\n"); cat(paste0("  ", head(ranked$id,5)), sep="\n"); cat("\n")
