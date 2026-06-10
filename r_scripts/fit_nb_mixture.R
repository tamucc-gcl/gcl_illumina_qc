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
  library(readr); library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
freq_file <- args[1]
knee_file <- args[2]

fallback <- tryCatch(as.integer(readLines(knee_file)[1]), error = function(e) 4L)
if (is.na(fallback)) fallback <- 4L

write_out <- function(val, msg) {
  val <- as.integer(round(val))
  if (is.na(val) || val < 2L) val <- 2L
  writeLines(as.character(val), "nb_cutoff1.value")
  writeLines(msg, "nb_mixture_fit.txt")
  cat(msg, "\n")
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
  "  knee fallback = %d\n"),
  mu1, th1, p1, mu2, th2, 1 - p1, as.integer(nb_cut), fallback)

write_out(nb_cut, msg)
