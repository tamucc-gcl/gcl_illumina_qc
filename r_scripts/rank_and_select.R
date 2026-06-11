#!/usr/bin/env Rscript
# rank_and_select.R  (PIVOT: r80-vs-n_contigs ELBOW selector)
# ============================================================================
# Selection principle (replaces the old weight-free rank aggregation):
#   The (c1,c2,init_sim,div_f,merge_r,final_sim) grid produces a SIZE CONTINUUM.
#   Every quality signal we tried (coverage CV/Gini, concordance, snps/locus) is
#   MONOTONE in assembly size, so none has an interior optimum. The ONE curve with
#   a real turning point is r80 (STACKS-style broadly-shared polymorphic loci) vs
#   n_contigs: r80 rises steeply then PLATEAUS — past the elbow you add contigs
#   (paralogs/junk) for ~no new broadly-shared loci.
#
#   WINNER = the elbow of r80(n_contigs): Kneedle on the r80 upper-envelope (the
#   point bowing furthest above the endpoints chord). Among candidates at/below the
#   elbow contig count, pick MAX r80, ties broken by FEWEST contigs (parsimony:
#   the smallest reference capturing the bulk of achievable broadly-shared loci).
#
#   concordance / anchor / NB / snps_per_locus are REPORTED for context but do NOT
#   drive selection (concordance is size-monotone; the r80 elbow is the real,
#   size-aware optimum). This is parameter-agnostic: it works no matter which axes
#   the user swept, because it operates on the OUTCOME (size<->r80), not any input.
#
# Usage:
#   Rscript rank_and_select.R cheap_all.tsv snp_all.tsv nb_cutoff1.value EXPECTED_LOCI MIN_SLOPE
#     cheap_all.tsv : id c1 c2 sim n_contigs total_len mean_len            (no header)
#     snp_all.tsv   : id c1 c2 sim concordance r80_loci snps_per_locus n_snps n_called_contigs
#     nb_cutoff1    : NB-derived cutoff1 integer (reported)
#     EXPECTED_LOCI : expected RAD loci, or "NA"/0 (reported as anchor)
#     MIN_SLOPE     : optional override — if >0, elbow = first contig level where
#                     marginal r80/1000-contigs (on the envelope) drops below it.
#                     Default behavior (<=0 or unset) = parameter-free Kneedle.
#
# Outputs: final_rank.tsv, best_id.value, optimize_plot.png (elbow curve),
#          optimize_params_plot.png (r80 faceted vs each tuned parameter)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(tidyr) })

args <- commandArgs(trailingOnly = TRUE)
cheap_file <- args[1]; snp_file <- args[2]; nb_file <- args[3]
expected   <- suppressWarnings(as.numeric(args[4]))
min_slope  <- suppressWarnings(as.numeric(args[5])); if (is.na(min_slope)) min_slope <- 0

nb_cut <- tryCatch(as.numeric(readLines(nb_file)[1]), error = function(e) NA_real_)

cheap <- read_tsv(cheap_file,
            col_names = c("id","c1","c2","isim","divf","mr","fsim",
                          "n_contigs","total_len","mean_len"),
            col_types = "ciiddidiii")

snp <- tryCatch(
  read_tsv(snp_file,
    col_names = c("id","c1","c2","isim","divf","mr","fsim",
                  "concordance","r80_loci","snps_per_locus","n_snps","n_called_contigs"),
    col_types = "ciiddidddddi", na = c("NA","")),
  error = function(e) NULL
)

empty_out <- function() {
  writeLines("id\tc1\tc2\tisim\tdivf\tmr\tfsim\tn_contigs\tconcordance\tr80_loci\tsnps_per_locus\tis_winner",
             "final_rank.tsv")
  writeLines("NA", "best_id.value")
  png("optimize_plot.png", width=700, height=500); plot.new(); text(.5,.5,"no candidates"); dev.off()
  png("optimize_params_plot.png", width=700, height=500); plot.new(); text(.5,.5,"no candidates"); dev.off()
  quit(save="no", status=0)
}
if (nrow(cheap) == 0) empty_out()

# join SNP signals onto candidate table
if (!is.null(snp) && nrow(snp) > 0) {
  d <- left_join(cheap, snp %>% select(id, concordance, r80_loci, snps_per_locus), by = "id")
} else {
  d <- cheap %>% mutate(concordance = NA_real_, r80_loci = NA_real_, snps_per_locus = NA_real_)
}

if (!any(is.finite(d$r80_loci))) {
  # No r80 signal -> cannot run the elbow; fall back to fewest contigs and warn.
  d <- d %>% arrange(n_contigs)
  best_id <- d$id[1]
  message("rank_and_select: WARNING no r80 signal; falling back to fewest-contigs pick.")
} else {
  # ---- order by size, build r80 UPPER ENVELOPE (running max) ----
  d <- d %>% arrange(n_contigs)
  env <- cummax(ifelse(is.na(d$r80_loci), -Inf, d$r80_loci))
  env[!is.finite(env)] <- NA_real_
  d$r80_envelope <- env

  x <- d$n_contigs; y <- d$r80_envelope
  ok <- is.finite(x) & is.finite(y)
  xk <- x[ok]; yk <- y[ok]

  elbow_nc <- NA_real_
  if (length(xk) >= 3 && diff(range(xk)) > 0 && diff(range(yk)) > 0) {
    if (min_slope > 0) {
      # marginal-slope override: first point where r80 gain per 1000 contigs < min_slope
      slope <- c(Inf, diff(yk) / (diff(xk)/1000))
      hit <- which(slope < min_slope)
      elbow_nc <- if (length(hit) > 0) xk[max(1, hit[1]-1)] else max(xk)
    } else {
      # parameter-free Kneedle: knee = max vertical distance ABOVE endpoints chord
      xn <- (xk - min(xk)) / (max(xk) - min(xk))
      yn <- (yk - min(yk)) / (max(yk) - min(yk))
      chord <- yn[1] + (yn[length(yn)] - yn[1]) * (xn - xn[1]) / (xn[length(xn)] - xn[1])
      elbow_nc <- xk[ which.max(yn - chord) ]
    }
  } else {
    elbow_nc <- if (length(xk)) min(xk) else min(d$n_contigs)
  }

  # ---- WINNER: among candidates at/below the elbow contig count, MAX r80;
  #       ties -> FEWEST contigs (parsimony) ----
  elig <- d %>% filter(n_contigs <= elbow_nc, is.finite(r80_loci))
  if (nrow(elig) == 0) elig <- d %>% filter(is.finite(r80_loci))
  winner_row <- elig %>% arrange(desc(r80_loci), n_contigs) %>% slice(1)
  best_id  <- winner_row$id
  attr(best_id, "elbow_nc") <- elbow_nc
}

d$is_winner <- d$id == best_id
writeLines(as.character(best_id), "best_id.value")

write_tsv(
  d %>% arrange(n_contigs) %>%
    select(id, c1, c2, isim, divf, mr, fsim, n_contigs,
           concordance, r80_loci, snps_per_locus, is_winner),
  "final_rank.tsv"
)

# ---------------- PLOT 1: the r80-vs-n_contigs elbow ----------------
# Encode the swept-parameter combination on each point. Identify which axes vary;
# map the (up to two) most-varied non-cutoff axes to color + shape so the user can
# see which parameter combo sits where on the size<->r80 curve. Falls back to a
# single color when the grid is large (too many combos to distinguish).
tryCatch({
  elbow_nc <- attr(best_id, "elbow_nc")
  wr <- d %>% filter(is_winner)

  axis_cols <- c("c1","c2","isim","divf","mr","fsim")
  varied <- axis_cols[ sapply(axis_cols, function(c) dplyr::n_distinct(d[[c]]) > 1) ]

  # choose a color axis and a shape axis from the varied ones (prefer divf/fsim for
  # color since cutoffs already drive the x-axis/size); keep it legible.
  pick <- function(prefer, pool) { p <- intersect(prefer, pool); if (length(p)) p[1] else if (length(pool)) pool[1] else NA }
  color_axis <- pick(c("divf","fsim","isim","mr"), varied)
  shape_axis <- pick(c("fsim","c2","c1","mr"), setdiff(varied, color_axis))

  too_many <- nrow(d) > 40   # if huge grid, don't overload aesthetics
  p1 <- ggplot(d, aes(n_contigs, r80_loci))
  if ("r80_envelope" %in% names(d))
    p1 <- p1 + geom_line(aes(y = r80_envelope), color = "grey60", linewidth = 0.4)
  if (is.finite(elbow_nc))
    p1 <- p1 + geom_vline(xintercept = elbow_nc, linetype = "dashed", color = "red")

  if (!too_many && !is.na(color_axis) && !is.na(shape_axis)) {
    p1 <- p1 + geom_point(aes(color = factor(.data[[color_axis]]),
                              shape = factor(.data[[shape_axis]])), size = 3, alpha = 0.9) +
      labs(color = color_axis, shape = shape_axis)
  } else if (!too_many && !is.na(color_axis)) {
    p1 <- p1 + geom_point(aes(color = factor(.data[[color_axis]])), size = 3, alpha = 0.9) +
      labs(color = color_axis)
  } else {
    p1 <- p1 + geom_point(color = "steelblue", size = 2.4, alpha = 0.85)
  }

  p1 <- p1 +
    geom_point(data = wr, color = "black", shape = 21, fill = "gold", size = 4.5, stroke = 1.2) +
    annotate("text", x = elbow_nc, y = min(d$r80_loci, na.rm=TRUE),
             label = paste0("  elbow ~ ", format(round(elbow_nc), big.mark=",")),
             color = "red", hjust = 0, vjust = 0) +
    labs(title = "r80 broadly-shared loci vs assembly size",
         subtitle = paste0("Winner (gold): ", best_id, " — fewest contigs at the r80 plateau"),
         x = "n_contigs (assembly size)",
         y = "r80 loci (polymorphic & genotyped in >=80% of samples)") +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    theme_bw() + theme(plot.title = element_text(face="bold"))
  ggsave("optimize_plot.png", p1, width = 9.5, height = 6, dpi = 150)
}, error = function(e) {
  png("optimize_plot.png", width=700, height=500); plot.new()
  text(.5,.5,paste("elbow plot failed:", conditionMessage(e)), col="grey40"); dev.off()
})

# ---------------- PLOT 2: r80 faceted vs each SWEPT parameter ----------------
# Facet over EVERY axis that varies (c1,c2,isim,divf,mr,fsim), so the user sees
# which knob actually moves r80. Color by n_contigs to expose size-confounding.
tryCatch({
  axis_cols <- c("c1","c2","isim","divf","mr","fsim")
  param_long <- d %>%
    select(id, r80_loci, n_contigs, all_of(axis_cols)) %>%
    pivot_longer(all_of(axis_cols), names_to = "parameter", values_to = "value") %>%
    group_by(parameter) %>% filter(dplyr::n_distinct(value) > 1) %>% ungroup() %>%
    mutate(parameter = factor(parameter, levels = axis_cols))
  if (nrow(param_long) > 0) {
    p2 <- ggplot(param_long, aes(factor(value), r80_loci)) +
      geom_jitter(aes(color = n_contigs), width = 0.12, height = 0, size = 2, alpha = 0.8) +
      facet_wrap(~ parameter, scales = "free_x") +
      scale_color_viridis_c(labels = scales::comma) +
      labs(title = "r80 vs each swept parameter",
           subtitle = "Which knob moves broadly-shared loci? (color = assembly size; if r80 tracks color, it's size-confounded)",
           x = NULL, y = "r80 loci", color = "n_contigs") +
      theme_bw() + theme(plot.title = element_text(face="bold"))
    ggsave("optimize_params_plot.png", p2, width = 10, height = 6, dpi = 150)
  } else {
    png("optimize_params_plot.png", width=700, height=400); plot.new()
    text(.5,.5,"No parameter was swept (all axes fixed)", col="grey40"); dev.off()
  }
}, error = function(e) {
  png("optimize_params_plot.png", width=700, height=400); plot.new()
  text(.5,.5,paste("param plot failed:", conditionMessage(e)), col="grey40"); dev.off()
})

cat(sprintf("rank_and_select(ELBOW): %d candidates | nb_cutoff1=%s | anchor=%s | winner=%s (n_contigs=%s, r80=%s)\n",
    nrow(d), ifelse(is.finite(nb_cut), as.character(nb_cut), "NA"),
    ifelse(is.finite(expected) && expected>0, as.character(expected), "disabled"),
    best_id,
    format(d$n_contigs[d$is_winner]), format(d$r80_loci[d$is_winner])))
