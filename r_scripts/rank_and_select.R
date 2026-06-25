#!/usr/bin/env Rscript
# rank_and_select.R  (curve-model selector: fitted r80(n_contigs) + diminishing returns)
# ============================================================================
# Selection principle:
#   The (c1,c2,init_sim,div_f,merge_r,final_sim) grid produces a SIZE CONTINUUM.
#   r80 (STACKS-style broadly-shared polymorphic loci) vs n_contigs rises steeply
#   then LEVELS OFF (and may gently decline once the assembly fills with paralogs/
#   junk/fragmented loci that are not broadly genotypeable). We want the
#   DIMINISHING-RETURNS point: the smallest assembly that already captures the bulk
#   of the achievable broadly-shared loci.
#
#   Earlier versions ran Kneedle on the RAW r80 envelope, which latched onto the
#   bottom-left corner (the most aggressive/worst assembly) when the curve was steep.
#   This version instead FITS A SMOOTH MODEL to the r80-vs-n_contigs envelope and
#   finds the leveling-off point ON THE FIT, which is robust to point-to-point noise
#   and behaves correctly for BOTH plateau-only and peaked curves:
#
#     1. Fit a monotone-saturating asymptotic model  r80 = Asym + (R0-Asym)*exp(-rate*nc)
#        (nls/SSasymp; robust fallbacks: loess, then raw envelope).
#     2. LEVELING-OFF point = smallest n_contigs where the fitted MARGINAL slope has
#        fallen below `knee_frac` of the curve's initial slope (default 0.10). This is
#        the diminishing-returns knee; it does NOT chase the peak (a saturating fit
#        plateaus, so points past the knee that decline sit BELOW the curve and are
#        not selected), and it generalizes to non-peaked runs.
#     3. WINNER = among candidates with n_contigs <= leveling-off point, the MAX r80;
#        ties -> FEWEST contigs (parsimony). If the curve genuinely peaks before the
#        knee, that peak candidate is selected (it is the max r80 in the eligible set).
#
#   The fitted curve + the leveling-off line + the winner are drawn on optimize_plot.png.
#   concordance / anchor / NB / snps_per_locus are REPORTED for context only.
#
# Usage:
#   Rscript rank_and_select.R cheap_all.tsv snp_all.tsv nb_cutoff1.value EXPECTED_LOCI KNEE_FRAC
#     KNEE_FRAC : marginal-slope fraction for the leveling-off point (default 0.10).
#                 Smaller => stricter "fully flat" (larger assembly); larger => earlier
#                 (smaller assembly). 0.10 = "90% of the rise is behind us".
#
# Outputs: final_rank.tsv, best_id.value, optimize_plot.png, optimize_params_plot.png

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(tidyr) })

args <- commandArgs(trailingOnly = TRUE)
cheap_file <- args[1]; snp_file <- args[2]; nb_file <- args[3]
expected   <- suppressWarnings(as.numeric(args[4]))
knee_frac  <- suppressWarnings(as.numeric(args[5])); if (is.na(knee_frac) || knee_frac <= 0) knee_frac <- 0.10

nb_cut <- tryCatch(as.numeric(readLines(nb_file)[1]), error = function(e) NA_real_)

cheap <- read_tsv(cheap_file,
            col_names = c("id","c1","c2","isim","divf","mr","fsim","n_contigs"),
            col_types = "ciiddidi")

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

if (!is.null(snp) && nrow(snp) > 0) {
  d <- left_join(cheap, snp %>% select(id, concordance, r80_loci, snps_per_locus), by = "id")
} else {
  d <- cheap %>% mutate(concordance = NA_real_, r80_loci = NA_real_, snps_per_locus = NA_real_)
}

# Curve-fit objects kept for plotting (NULL unless a fit succeeds)
fit_grid   <- NULL          # data.frame(n_contigs, r80_fit) for the drawn curve
level_nc   <- NA_real_      # leveling-off n_contigs
fit_method <- "none"

if (!any(is.finite(d$r80_loci))) {
  d <- d %>% arrange(n_contigs)
  best_id <- d$id[1]
  message("rank_and_select: WARNING no r80 signal; falling back to fewest-contigs pick.")
} else {
  # ---- order by size, build r80 UPPER ENVELOPE (running max) ----
  d <- d %>% arrange(n_contigs)
  env <- cummax(ifelse(is.na(d$r80_loci), -Inf, d$r80_loci))
  env[!is.finite(env)] <- NA_real_
  d$r80_envelope <- env

  e  <- d %>% filter(is.finite(n_contigs), is.finite(r80_envelope))
  xk <- e$n_contigs; yk <- e$r80_envelope

  # dense x grid for prediction / knee search
  xs <- if (length(unique(xk)) >= 2) seq(min(xk), max(xk), length.out = 400) else xk

  pred <- NULL
  if (length(unique(xk)) >= 3 && diff(range(xk)) > 0 && diff(range(yk)) > 0) {
    # ---- 1) primary: asymptotic (monotone-saturating) via SSasymp ----
    pred <- tryCatch({
      fm <- nls(yk ~ SSasymp(xk, Asym, R0, lrc))
      p  <- predict(fm, newdata = data.frame(xk = xs))
      fit_method <<- "asymptotic (SSasymp)"
      p
    }, error = function(e1) {
      # ---- 2) fallback: loess (span tuned to point count) ----
      tryCatch({
        sp <- if (length(xk) >= 6) 0.9 else 1.2
        lo <- loess(yk ~ xk, span = sp, degree = 1,
                    control = loess.control(surface = "direct"))
        p  <- predict(lo, newdata = data.frame(xk = xs))
        fit_method <<- "loess"
        p
      }, error = function(e2) NULL)
    })
  }

  if (!is.null(pred) && all(is.finite(pred))) {
    fit_grid <- data.frame(n_contigs = xs, r80_fit = pred)
    # ---- 2) leveling-off = first x where fitted marginal slope < knee_frac * initial slope
    dy <- c(NA, diff(pred) / diff(xs))
    init_slope <- dy[2]
    if (is.finite(init_slope) && init_slope > 0) {
      below <- which(is.finite(dy) & dy < knee_frac * init_slope)
      level_nc <- if (length(below) > 0) xs[below[1]] else max(xs)
    } else {
      # non-increasing fit (already flat / declining): take the fitted-max location
      level_nc <- xs[ which.max(pred) ]
    }
  } else {
    # ---- 3) raw fallback: parameter-free Kneedle on the envelope (legacy) ----
    fit_method <- "kneedle (raw, no fit converged)"
    if (length(xk) >= 3) {
      xn <- (xk - min(xk)) / (max(xk) - min(xk))
      yn <- (yk - min(yk)) / (max(yk) - min(yk))
      chord <- yn[1] + (yn[length(yn)] - yn[1]) * (xn - xn[1]) / (xn[length(xn)] - xn[1])
      level_nc <- xk[ which.max(yn - chord) ]
    } else level_nc <- min(xk)
  }

  # ---- WINNER: snap the leveling-off point to the NEAREST actual candidate size
  #       (so the candidate that DEFINES the knee is itself eligible — a strict
  #       n_contigs < level_nc test would exclude it and pick one size too small),
  #       then take MAX r80 among candidates at/below that size; ties -> FEWEST contigs.
  level_snapped <- level_nc
  if (is.finite(level_nc)) {
    sizes <- d$n_contigs[is.finite(d$n_contigs)]
    if (length(sizes)) level_snapped <- sizes[ which.min(abs(sizes - level_nc)) ]
  }
  elig <- d %>% filter(n_contigs <= level_snapped + 1e-6, is.finite(r80_loci))
  if (nrow(elig) == 0) elig <- d %>% filter(is.finite(r80_loci))
  winner_row <- elig %>% arrange(desc(r80_loci), n_contigs) %>% slice(1)
  best_id  <- winner_row$id

  # ---- CONFIRMATION: guard against the fitted knee UNDERSHOOTING the true optimum.
  # The asymptotic fit runs on the running-max ENVELOPE, so a genuine DECLINE past a
  # peak is masked (cummax flattens it); combined with snapping the knee to the
  # nearest candidate, the knee can land LEFT of a materially-better candidate,
  # excluding it. Because the final pick is data-driven, we cross-check the
  # knee-winner against the RAW r80 maximum: if the knee-winner's r80 is more than
  # CONFIRM_TOL below the best observed r80, the knee undershot -> OVERRIDE to the
  # SMALLEST candidate whose r80 is within CONFIRM_TOL of the best (parsimonious
  # near-best; still does NOT chase size — among near-best it takes the fewest
  # contigs). Otherwise the knee-winner is confirmed and stands (we do NOT move to a
  # smaller near-best, to respect the leveling-off criterion).
  CONFIRM_TOL <- 0.02
  r80_max <- max(d$r80_loci, na.rm = TRUE)
  win_r80 <- d$r80_loci[d$id == best_id][1]
  near_best <- d %>% filter(is.finite(r80_loci), r80_loci >= (1 - CONFIRM_TOL) * r80_max) %>%
                 arrange(n_contigs)
  smallest_near_best <- near_best$id[1]
  if (is.finite(win_r80) && is.finite(r80_max) && win_r80 < (1 - CONFIRM_TOL) * r80_max
      && !is.na(smallest_near_best)) {
    message(sprintf(
      "rank_and_select: CONFIRM FAILED — knee-winner %s (r80=%s) is >%.0f%% below best observed r80=%s; OVERRIDING to smallest near-best %s (r80=%s, n_contigs=%s).",
      best_id, format(win_r80), CONFIRM_TOL*100, format(r80_max),
      smallest_near_best, format(d$r80_loci[d$id==smallest_near_best]),
      format(d$n_contigs[d$id==smallest_near_best])))
    best_id <- smallest_near_best
  } else {
    message(sprintf(
      "rank_and_select: CONFIRMED knee-winner %s (r80=%s) within %.0f%% of best observed r80=%s.",
      best_id, format(win_r80), CONFIRM_TOL*100, format(r80_max)))
  }
}

d$is_winner <- d$id == best_id
writeLines(as.character(best_id), "best_id.value")

write_tsv(
  d %>% arrange(n_contigs) %>%
    select(id, c1, c2, isim, divf, mr, fsim, n_contigs,
           concordance, r80_loci, snps_per_locus, is_winner),
  "final_rank.tsv"
)

# ---------------- PLOT 1: r80 vs n_contigs with the FITTED curve ----------------
tryCatch({
  wr <- d %>% filter(is_winner)
  axis_cols <- c("c1","c2","isim","divf","mr","fsim")
  varied <- axis_cols[ sapply(axis_cols, function(c) dplyr::n_distinct(d[[c]]) > 1) ]
  pick <- function(prefer, pool){ p<-intersect(prefer,pool); if(length(p)) p[1] else if(length(pool)) pool[1] else NA }
  color_axis <- pick(c("divf","fsim","isim","mr"), varied)
  shape_axis <- pick(c("fsim","c2","c1","mr"), setdiff(varied, color_axis))
  too_many <- nrow(d) > 60

  p1 <- ggplot(d, aes(n_contigs, r80_loci))
  # fitted model curve (the key addition)
  if (!is.null(fit_grid))
    p1 <- p1 + geom_line(data = fit_grid, aes(n_contigs, r80_fit),
                         color = "grey35", linewidth = 0.9)
  # leveling-off line (red dashed)
  if (is.finite(level_nc))
    p1 <- p1 + geom_vline(xintercept = level_nc, linetype = "dashed", color = "red")
  # biological anchor: expected ddRAD locus count (blue dotted), if enabled
  if (is.finite(expected) && expected > 0)
    p1 <- p1 + geom_vline(xintercept = expected, linetype = "dotted",
                          color = "steelblue", linewidth = 0.8)

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
    annotate("text", x = level_nc, y = min(d$r80_loci, na.rm=TRUE),
             label = paste0("  levels off ~ ", format(round(level_nc), big.mark=",")),
             color = "red", hjust = 0, vjust = 0)
  if (is.finite(expected) && expected > 0)
    p1 <- p1 + annotate("text", x = expected, y = max(d$r80_loci, na.rm=TRUE),
                        label = paste0("expected ~ ", format(round(expected), big.mark=","), "  "),
                        color = "steelblue", hjust = 1, vjust = 1)
  p1 <- p1 +
    labs(title = "r80 broadly-shared loci vs assembly size",
         subtitle = paste0("Winner (gold): ", best_id,
                           "  |  fit: ", fit_method,
                           "  |  leveling-off at ", round(knee_frac*100), "% of initial slope"),
         x = "n_contigs (assembly size)",
         y = "r80 loci (polymorphic & genotyped in >=80% of samples)") +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    theme_bw() + theme(plot.title = element_text(face="bold"))
  ggsave("optimize_plot.png", p1, width = 9.5, height = 6, dpi = 150)
}, error = function(e) {
  png("optimize_plot.png", width=700, height=500); plot.new()
  text(.5,.5,paste("fit plot failed:", conditionMessage(e)), col="grey40"); dev.off()
})

# ---------------- PLOT 2: r80 faceted vs each SWEPT parameter ----------------
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

cat(sprintf("rank_and_select(CURVE-FIT): %d candidates | fit=%s | nb_cutoff1=%s | anchor=%s | levels_off~%s | winner=%s (n_contigs=%s, r80=%s)\n",
    nrow(d), fit_method,
    ifelse(is.finite(nb_cut), as.character(nb_cut), "NA"),
    ifelse(is.finite(expected) && expected>0, as.character(expected), "disabled"),
    ifelse(is.finite(level_nc), format(round(level_nc)), "NA"),
    best_id,
    format(d$n_contigs[d$is_winner]), format(d$r80_loci[d$is_winner])))
