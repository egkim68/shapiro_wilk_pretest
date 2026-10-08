# REAL-DATA DEMONSTRATION: the Shapiro-Wilk pre-test on warpbreaks before and after missing data.
# Base R only; uses the built-in `warpbreaks` data (number of warp breaks per loom, wool A vs wool B,
# 27 looms each), a right-skewed count. Writes demo_warpbreaks.csv to <run>/analysis (created if missing).
# Deletion is repeated many times at a 30% missing rate under three mechanisms, one rule for both wools:
#   MCAR              every value equally likely to go missing
#   tail lost (MNAR)  larger counts more likely to go missing (logistic in the standardized value)
#   low end lost      smaller counts more likely to go missing
# For each deletion the per-group Shapiro-Wilk pre-test is run on what is left, and the test it would
# choose (Welch t if it passes, Mann-Whitney if not) is recorded.
# Use:  source("demo_warpbreaks.R")   or   Rscript demo_warpbreaks.R

OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer1")
# the table goes into the analysis folder of the current run, so it is part of the upload package
DEMO_DIR <- file.path(OUTPUT_ROOT, paste0("pretest_", Sys.getenv("SIM_MODE", unset = "medium")), "analysis")
N_DELETIONS <- 2000
RATE <- 0.30
STRENGTH <- 1.5
SEED <- 2041

calibrated_missing <- function(eta, rate) {
  intercept <- uniroot(function(a) mean(plogis(a + eta)) - rate, c(-max(eta) - 40, -min(eta) + 40), tol = 1e-8)$root
  runif(length(eta)) < plogis(intercept + eta)
}
sw_group <- function(a, b) min(shapiro.test(a)$p.value, shapiro.test(b)$p.value)

run_demo <- function() {
  dir.create(DEMO_DIR, recursive = TRUE, showWarnings = FALSE)
  a0 <- warpbreaks$breaks[warpbreaks$wool == "A"]; b0 <- warpbreaks$breaks[warpbreaks$wool == "B"]
  skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
  cat(sprintf("warpbreaks: wool A n=%d (mean %.1f, skewness %.2f), wool B n=%d (mean %.1f, skewness %.2f)\n",
              length(a0), mean(a0), skew(a0), length(b0), mean(b0), skew(b0)))
  full <- c(sw_A = shapiro.test(a0)$p.value, sw_B = shapiro.test(b0)$p.value,
            welch = t.test(a0, b0)$p.value, mw = suppressWarnings(wilcox.test(a0, b0, exact = FALSE)$p.value))
  cat(sprintf("Complete data: Shapiro-Wilk p = %.4f (A), %.4f (B) -> per-group pre-test %s; Welch p = %.3f, Mann-Whitney p = %.3f\n",
              full[1], full[2], if (min(full[1:2]) < 0.05) "FAILS (Mann-Whitney chosen)" else "passes (t test chosen)", full[3], full[4]))
  set.seed(SEED)
  mechs <- list(MCAR = function(y) runif(length(y)) < RATE,
                tail_lost = function(y) calibrated_missing(STRENGTH * (y - mean(y)) / sd(y), RATE),
                low_end_lost = function(y) calibrated_missing(-STRENGTH * (y - mean(y)) / sd(y), RATE))
  rows <- list()
  for (m in names(mechs)) {
    res <- t(replicate(N_DELETIONS, {
      g <- mechs[[m]](c(a0, b0)); a <- a0[!g[seq_along(a0)]]; b <- b0[!g[-seq_along(a0)]]
      if (length(a) < 4 || length(b) < 4) c(NA, NA, NA) else {
        p <- sw_group(a, b)
        c(reject = p < 0.05, skew_a = skew(a), skew_b = skew(b)) }
    }))
    res <- res[complete.cases(res), , drop = FALSE]
    rows[[m]] <- data.frame(mechanism = m, deletions = nrow(res), sw_rejects = mean(res[, 1]),
                            t_test_chosen = 1 - mean(res[, 1]), mean_skew_A = mean(res[, 2]), mean_skew_B = mean(res[, 3]))
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out <- rbind(data.frame(mechanism = "complete data", deletions = 0, sw_rejects = as.numeric(min(full[1:2]) < 0.05),
                          t_test_chosen = as.numeric(min(full[1:2]) >= 0.05), mean_skew_A = skew(a0), mean_skew_B = skew(b0)), out)
  write.csv(out, file.path(DEMO_DIR, "demo_warpbreaks.csv"), row.names = FALSE)
  cat(sprintf("\nAfter %d deletions of %.0f%% of the values in each wool:\n", N_DELETIONS, 100 * RATE))
  for (i in seq_len(nrow(out))) cat(sprintf("  %-14s pre-test rejects %.2f, t test chosen %.2f, mean skewness A %.2f, B %.2f\n",
    out$mechanism[i], out$sw_rejects[i], out$t_test_chosen[i], out$mean_skew_A[i], out$mean_skew_B[i]))
  cat("\nWritten to", normalizePath(file.path(DEMO_DIR, "demo_warpbreaks.csv")), "\n")
  invisible(out)
}
demo_results <- run_demo()
