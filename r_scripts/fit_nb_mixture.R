#!/usr/bin/env Rscript
# fit_nb_mixture.R
# Fit a 2-component negative-binomial mixture to the pooled within-individual
# coverage frequency distribution and return a principled cutoff1 = the coverage
# at which the high-mean ("true locus") component overtakes the low-mean
# ("error/noise") component in posterior probability.
#
# Usage: Rscript fit_nb_mixture.R coverage_freq.txt knee_fallback.value
#   coverage_freq.txt : two cols "<n_records> <coverage>" (from sort|uniq -c)
#   knee_fallback     : file with the geometric-knee cutoff1 (fallback value)
#
# Outputs:
#   nb_cutoff1.value    integer cutoff1
#   nb_mixture_fit.txt  human-readable fit summary
#
# Robustness: any failure (too few points, non-convergence, degenerate fit)
# falls back to the knee value so the pipeline never breaks on this signal.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
freq_file <- args[1]
knee_file <- args[2]

fallback <- tryCatch(as.integer(readLines(knee_file)[1]), error = function(e) 4L)
if (is.na(fallback)) fallback <- 4L

PLOT_FILE <- "nb_mixture_plot.png"

# Emit a placeholder plot when we cannot/do not fit (keeps the output present so
# the channel cardinality and report wiring stay valid).
placeholder_plot <- function(msg) {
  png(PLOT_FILE, width = 800, height = 600)
  plot.new(); text(0.5, 0.5, msg, cex = 1.4, col = "grey40")
  dev.off()
}

# Full overlay plot: observed coverage histogram + the two fitted NB component
# densities (scaled to counts) + the crossover cutoff1 line.
fit_plot <- function(freq_df, mu1, th1, mu2, th2, p1, cutoff) {
  total <- sum(freq_df$n)
  xs <- seq(min(freq_df$cov), max(freq_df$cov))
  comp1 <- p1       * dnbinom(xs, size = th1, mu = mu1) * total
  comp2 <- (1 - p1) * dnbinom(xs, size = th2, mu = mu2) * total
  dens <- bind_rows(
    tibble(cov = xs, count = comp1, component = sprintf("error (mu=%.1f)", mu1)),
    tibble(cov = xs, count = comp2, component = sprintf("locus (mu=%.1f)", mu2))
  )
  # On a log y-axis, values < 1 (incl. component densities in the tail and the
  # implicit 0 baseline of geom_col) go negative/-Inf and render as a degenerate
  # axis (the "1 0 0 0" labels) or vanish. Floor everything at 1 and use explicit
  # log breaks with log-aware labels so the spike AND the components are visible.
  freq_plot <- freq_df %>% mutate(n = pmax(n, 1))
  dens <- dens %>% filter(is.finite(count)) %>% mutate(count = pmax(count, 1))
  ymax <- max(freq_plot$n, dens$count, na.rm = TRUE)

  p <- ggplot() +
    geom_col(data = freq_plot, aes(cov, n), fill = "grey80", color = NA, width = 1) +
    geom_line(data = dens, aes(cov, count, color = component), linewidth = 1) +
    geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", linewidth = 0.9) +
    annotate("text", x = cutoff, y = ymax,
             label = paste0("  cutoff1 = ", cutoff), color = "red", hjust = 0, vjust = 1) +
    scale_color_manual(values = c("steelblue", "darkorange"), name = NULL) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    annotation_logticks(sides = "l") +
    coord_cartesian(xlim = c(0, min(max(freq_df$cov), 60)), ylim = c(1, ymax)) +
    labs(title = "NB-mixture fit on within-individual coverage",
         subtitle = "Grey = observed (log scale); lines = fitted NB components; red = posterior crossover (cutoff1)",
         x = "Within-individual coverage (reads per unique sequence)",
         y = "Number of unique sequences (log scale)") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "top")
  ggsave(PLOT_FILE, p, width = 8, height = 6, dpi = 200)
}

write_out <- function(val, msg, freq_df = NULL, fit = NULL) {
  val <- as.integer(round(val))
  if (is.na(val) || val < 2L) val <- 2L
  writeLines(as.character(val), "nb_cutoff1.value")
  writeLines(msg, "nb_mixture_fit.txt")
  cat(msg, "\n")
  # plot: full overlay if we have a fit, else a labelled placeholder
  if (!is.null(freq_df) && !is.null(fit)) {
    tryCatch(
      fit_plot(freq_df, fit$mu1, fit$th1, fit$mu2, fit$th2, fit$p1, val),
      error = function(e) placeholder_plot(paste("NB plot failed:", conditionMessage(e)))
    )
  } else {
    placeholder_plot(sub("\n.*", "", msg))
  }
  quit(save = "no", status = 0)
}

# ---- read frequency table -> expand to a coverage vector (weighted) ----
freq <- tryCatch(
  read_table(freq_file, col_names = c("n", "cov"), col_types = "ii"),
  error = function(e) NULL
)
if (is.null(freq) || nrow(freq) == 0)
  write_out(fallback, paste0("NB-mixture: no coverage data; fallback=", fallback))

freq <- freq %>% filter(!is.na(n), !is.na(cov), n > 0, cov > 0)
# Cap coverage to a sane ceiling so a few ultra-high-coverage repeats don't
# dominate the fit (matches the cutoff1<=40-ish regime).
freq <- freq %>% filter(cov <= 200)
if (nrow(freq) < 4)
  write_out(fallback, paste0("NB-mixture: too few points; fallback=", fallback))

# Method-of-moments style 2-component NB mixture via EM on the binned counts.
# We avoid heavy deps (no mixtools/flexmix) — implement a compact EM on the
# weighted histogram. NB parameterized by mean mu and size theta (dispersion).
cov <- freq$cov
wt  <- as.numeric(freq$n)

dnb <- function(x, mu, theta) dnbinom(x, size = theta, mu = mu)

# init: split at median coverage into low/high components
med <- max(2, round(weighted.mean(cov, wt)))
lo  <- cov <= med
mu1 <- max(1, weighted.mean(cov[lo],  wt[lo]));  if (is.nan(mu1)) mu1 <- 1
mu2 <- max(mu1 + 1, weighted.mean(cov[!lo], wt[!lo])); if (is.nan(mu2)) mu2 <- mu1 + 2
th1 <- 2; th2 <- 2
p1  <- 0.5

ok <- TRUE
for (iter in 1:200) {
  d1 <- p1 * dnb(cov, mu1, th1)
  d2 <- (1 - p1) * dnb(cov, mu2, th2)
  den <- d1 + d2
  den[den <= 0 | !is.finite(den)] <- .Machine$double.eps
  r1 <- d1 / den                 # responsibility for component 1 (low)
  r1[!is.finite(r1)] <- 0.5

  # weighted M-step for means and mixing
  w1 <- wt * r1; w2 <- wt * (1 - r1)
  sw1 <- sum(w1); sw2 <- sum(w2)
  if (sw1 <= 0 || sw2 <= 0) { ok <- FALSE; break }
  mu1n <- sum(w1 * cov) / sw1
  mu2n <- sum(w2 * cov) / sw2
  p1n  <- sw1 / (sw1 + sw2)

  # crude dispersion update via moment match (variance/mean relationship)
  v1 <- sum(w1 * (cov - mu1n)^2) / sw1
  v2 <- sum(w2 * (cov - mu2n)^2) / sw2
  th1 <- ifelse(v1 > mu1n, mu1n^2 / (v1 - mu1n), 50)
  th2 <- ifelse(v2 > mu2n, mu2n^2 / (v2 - mu2n), 50)
  th1 <- max(0.05, min(th1, 1e4)); th2 <- max(0.05, min(th2, 1e4))

  if (abs(mu1n - mu1) < 1e-4 && abs(mu2n - mu2) < 1e-4 && abs(p1n - p1) < 1e-5) {
    mu1 <- mu1n; mu2 <- mu2n; p1 <- p1n; break
  }
  mu1 <- mu1n; mu2 <- mu2n; p1 <- p1n
}

if (!ok || !is.finite(mu1) || !is.finite(mu2))
  write_out(fallback, paste0("NB-mixture: EM failed; fallback=", fallback))

# Ensure component 1 is the LOW-mean (error) component
if (mu1 > mu2) { tmp<-mu1;mu1<-mu2;mu2<-tmp; tmp<-th1;th1<-th2;th2<-tmp; p1<-1-p1 }

# cutoff1 = smallest coverage where the HIGH component posterior >= low component
xs <- seq(2, max(cov))
post_hi <- ((1 - p1) * dnb(xs, mu2, th2)) /
           (p1 * dnb(xs, mu1, th1) + (1 - p1) * dnb(xs, mu2, th2))
post_hi[!is.finite(post_hi)] <- 0
cross <- xs[which(post_hi >= 0.5)]
nb_cut <- if (length(cross) > 0) min(cross) else fallback

# sanity clamp: cutoff1 in a reasonable RAD regime
if (!is.finite(nb_cut) || nb_cut < 2) nb_cut <- fallback
nb_cut <- min(nb_cut, 40L)

msg <- sprintf(paste0("NB-mixture cutoff1 fit\n",
  "  low component:  mu=%.2f theta=%.2f (mix=%.2f)\n",
  "  high component: mu=%.2f theta=%.2f (mix=%.2f)\n",
  "  posterior crossover (cutoff1) = %d\n",
  "  knee fallback = %d\n",
  "  input: %d coverage bins, %s total unique sequences, max count = %s\n"),
  mu1, th1, p1, mu2, th2, 1 - p1, as.integer(nb_cut), fallback,
  nrow(freq), format(sum(freq$n), big.mark=","), format(max(freq$n), big.mark=","))

write_out(nb_cut, msg,
          freq_df = freq,
          fit = list(mu1 = mu1, th1 = th1, mu2 = mu2, th2 = th2, p1 = p1))
