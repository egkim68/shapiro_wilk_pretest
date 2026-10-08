# CHECK THE FROZEN HYPOTHESES (PROTOCOL_pretest_v1.md) ON THE CONFIRMATORY RUN. Base R only; changes no result.
# Reads the tables that summarize_pretest.R wrote and prints, for H1 to H5, every criterion with its value and a
# verdict. Writes <run>/analysis/confirmation_pretest_<mode>.csv.
# Use:  Rscript check_pretest_hypotheses.R                (the run of SIM_MODE, default full)
#       Rscript check_pretest_hypotheses.R <run folder>   (a run folder, or a flat folder of an upload package)

args <- commandArgs(trailingOnly = TRUE)
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer1")
MODE <- Sys.getenv("SIM_MODE", unset = "full")
RUN_DIR <- if (length(args) >= 1) args[1] else file.path(OUTPUT_ROOT, paste0("pretest_", MODE))
TWO_STAGE <- c("Two-stage: per group, Student t or MW", "Two-stage: per group, Welch t or MW",
               "Two-stage: residuals, Student t or MW", "Two-stage: residuals, Welch t or MW")

read_table <- function(dir, name) {
  f <- c(file.path(dir, "analysis", name), file.path(dir, name), list.files(dir, pattern = paste0("__", name, "$"), full.names = TRUE))
  f <- f[file.exists(f)]
  if (!length(f)) stop("Table ", name, " not found in ", dir, ". Run summarize_pretest.R first.", call. = FALSE)
  read.csv(f[1], stringsAsFactors = FALSE)
}

check <- function(run_dir) {
  sw <- read_table(run_dir, "pretest_sw_rates.csv"); er <- read_table(run_dir, "pretest_error.csv")
  mk <- read_table(run_dir, "pretest_masking.csv"); cd <- read_table(run_dir, "pretest_conditional.csv")
  dg <- read_table(run_dir, "pretest_diagnostic.csv")
  chg <- function(shape, mech) sw$group_change[sw$shape == shape & sw$mechanism == mech]
  harm <- function(g, p) er$harmful_share[er$group == g & er$procedure == p]
  merr <- function(g, p) er$mean_error[er$group == g & er$procedure == p]
  rise <- function(mech) { r <- mk[mk$mechanism == mech, ]; r$sw_pass_observed - r$sw_pass_complete }
  cond <- function(shape, which) cd$conditional_error[cd$label == paste0("Cond: Welch t | per group ", which) & cd$shape == shape]
  auc <- function(subset) dg$auc_per_group_test[dg$subset == subset]
  rows <- list()
  add <- function(h, criterion, value, ok) rows[[length(rows) + 1]] <<- data.frame(hypothesis = h, criterion = criterion,
                                                                                    value = round(value, 4), met = isTRUE(ok), stringsAsFactors = FALSE)   # a missing value counts as not met
  # H1 masking: upper-tail loss lowers detection of skew by at least 10 points, and by 5 points more than MCAR
  for (s in c("mild_skew", "skewed")) {
    v <- chg(s, "MNAR_large_lost"); add("H1", paste0(s, ": change under upper-tail loss <= -0.10"), v, v <= -0.10)
    d <- v - chg(s, "MCAR"); add("H1", paste0(s, ": change under upper-tail loss minus change under MCAR <= -0.05"), d, d <= -0.05)
  }
  # H2 random loss: MCAR and group-dependent loss lower detection of every non-normal shape by at least 3 points;
  # normal populations move by at most 4 points under every mechanism
  for (s in c("mild_skew", "skewed", "heavy")) for (m in c("MCAR", "MAR_group")) {
    v <- chg(s, m); add("H2", paste0(s, ", ", m, ": change <= -0.03"), v, v <= -0.03) }
  for (m in c("MCAR", "MAR_group", "MNAR_large_lost", "MNAR_small_lost")) {
    v <- chg("normal", m); add("H2", paste0("normal, ", m, ": |change| <= 0.04"), v, abs(v) <= 0.04) }
  # H3 no benefit: Welch t alone has a lower harmful share and (identical distributions) a lower mean error than every two-stage procedure
  for (g in c("identical distributions", "all")) for (p in TWO_STAGE) {
    d <- harm(g, "Welch t") - harm(g, p); add("H3", paste0(g, ": harmful share of Welch t minus ", p, " < 0"), d, d < 0) }
  for (p in TWO_STAGE) {
    d <- merr("identical distributions", "Welch t") - merr("identical distributions", p)
    add("H3", paste0("identical distributions: mean error of Welch t minus ", p, " < 0"), d, d < 0) }
  # H4 hidden skew: upper-tail loss raises the pass rate by at least 15 points, 10 points more than random loss;
  # on skewed populations Welch t errs more after a pass than after a fail, by at least 0.02
  v <- rise("MNAR_large_lost"); add("H4", "rise in pass rate under upper-tail loss >= 0.15", v, v >= 0.15)
  for (m in c("MCAR", "MAR_group")) { d <- v - rise(m); add("H4", paste0("rise under upper-tail loss minus rise under ", m, " >= 0.10"), d, d >= 0.10) }
  for (s in c("mild_skew", "skewed")) { d <- cond(s, "pass") - cond(s, "fail"); add("H4", paste0(s, ": Welch t error after pass minus after fail >= 0.02"), d, d >= 0.02) }
  # H5 no warning: held-out AUC of the observed rejection rate for outcome-driven missingness
  v <- auc("all designs"); add("H5", "all designs: held-out AUC within 0.35 to 0.65", v, v > 0.35 && v < 0.65)
  v <- auc("normal populations, no outliers"); add("H5", "normal populations, no outliers: held-out AUC < 0.70", v, v < 0.70)
  out <- do.call(rbind, rows)
  verdict <- tapply(out$met, out$hypothesis, all)
  cat(sprintf("Frozen hypotheses (PROTOCOL_pretest_v1.md) checked on %s\n\n", normalizePath(run_dir, winslash = "/")))
  for (h in names(verdict)) {
    cat(sprintf("%s  %s\n", h, if (verdict[[h]]) "SUPPORTED (every criterion met)" else "NOT SUPPORTED"))
    x <- out[out$hypothesis == h, ]
    for (i in seq_len(nrow(x))) cat(sprintf("    [%s] %-95s %8.4f\n", if (x$met[i]) "met" else "NOT", x$criterion[i], x$value[i]))
  }
  dest <- if (dir.exists(file.path(run_dir, "analysis"))) file.path(run_dir, "analysis") else run_dir
  mode <- if (grepl("full", run_dir)) "full" else if (grepl("medium", run_dir)) "medium" else MODE
  write.csv(out, file.path(dest, paste0("confirmation_pretest_", mode, ".csv")), row.names = FALSE)
  invisible(out)
}
confirmation <- check(RUN_DIR)
