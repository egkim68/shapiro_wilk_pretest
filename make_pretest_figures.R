# FIGURES FOR THE PAPER (PNG). Base R only. Runs nothing and changes no result: it reads the tables that
# summarize_pretest.R wrote to <run>/analysis and draws three images into <run>/paper_figures (created if
# missing). Rerun any time, for example after the confirmatory run.
#   Figure_1_sw_complete_vs_observed.png   Shapiro-Wilk rejection rate, complete -> observed data
#   Figure_2_hidden_skew.png               (A) pre-test pass rate on skewed data, (B) t test error after pass or fail
#   Figure_3_harmful_share.png             harmful share: Welch t alone vs the two-stage procedure, identical distributions
#                                          vs equal means with unequal spreads (computed from the counts)
# Use:  Rscript make_pretest_figures.R                (the run of SIM_MODE, default medium)
#       Rscript make_pretest_figures.R <run folder>   (a run folder, or a flat folder of an upload package)

args <- commandArgs(trailingOnly = TRUE)
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer1")
MODE <- Sys.getenv("SIM_MODE", unset = "medium")
RUN_DIR <- if (length(args) >= 1) args[1] else file.path(OUTPUT_ROOT, paste0("pretest_", MODE))
FIG_WIDTH_IN <- 6.5
FIG_DPI <- 300
COL <- c(complete = "#999999", observed = "#0072B2", welch = "#0072B2", two_stage = "#D55E00", pass = "#D55E00", fail = "#009E73")
SHAPE_NAMES <- c(normal = "Normal", mild_skew = "Mild skew", skewed = "Strong skew", heavy = "Heavy-tailed")
MECH_NAMES <- c(MCAR = "MCAR", MAR_group = "Group-dependent", MNAR_large_lost = "Upper tail lost", MNAR_small_lost = "Lower tail lost")

read_table <- function(dir, name) {
  f <- c(file.path(dir, "analysis", name), file.path(dir, name), list.files(dir, pattern = paste0("__", name, "$"), full.names = TRUE))
  f <- f[file.exists(f)]
  if (!length(f)) stop("Table ", name, " not found in ", dir, ". Run summarize_pretest.R first.", call. = FALSE)
  read.csv(f[1], stringsAsFactors = FALSE)
}
open_png <- function(dir, name, height_in) {
  f <- file.path(dir, paste0(name, ".png"))
  png(f, width = FIG_WIDTH_IN, height = height_in, units = "in", res = FIG_DPI, pointsize = 11)
  f
}

make_figures <- function(run_dir) {
  out <- file.path(run_dir, "paper_figures"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out)) stop("Cannot create ", out, call. = FALSE)
  made <- character(0)

  # ---- Figure 1: complete -> observed Shapiro-Wilk rejection rate, by shape and mechanism --------------
  sw <- read_table(run_dir, "pretest_sw_rates.csv")
  sw$shape <- factor(sw$shape, levels = names(SHAPE_NAMES)); sw$mechanism <- factor(sw$mechanism, levels = names(MECH_NAMES))
  sw <- sw[order(sw$shape, sw$mechanism), ]
  f <- open_png(out, "Figure_1_sw_complete_vs_observed", 4.6); made <- c(made, f)
  par(mar = c(3.4, 13.5, 0.8, 1), mgp = c(2.1, 0.6, 0), tcl = -0.3, las = 1)
  y <- rev(seq_len(nrow(sw)) + (as.integer(sw$shape) - 1) * 0.8)
  plot(NA, xlim = c(0, 1), ylim = range(y) + c(-0.6, 0.6), axes = FALSE, xlab = "Share of replications in which Shapiro-Wilk rejects normality", ylab = "")
  abline(v = seq(0, 1, 0.2), col = "gray92"); abline(v = 0.05, lty = 3, col = "gray50")
  arrows(sw$group_complete, y, sw$group_observed, y, length = 0.05, col = "gray40", lwd = 1.2)
  points(sw$group_complete, y, pch = 21, bg = COL[["complete"]], col = "black", cex = 1.1)
  points(sw$group_observed, y, pch = 21, bg = COL[["observed"]], col = "black", cex = 1.1)
  axis(1, at = seq(0, 1, 0.2)); axis(2, at = y, labels = MECH_NAMES[as.character(sw$mechanism)], tick = FALSE, cex.axis = 1)
  for (s in levels(sw$shape)) { yy <- y[sw$shape == s]; if (length(yy)) mtext(SHAPE_NAMES[[s]], side = 2, at = mean(yy), line = 8.3, font = 2, cex = 0.95, adj = 1) }
  legend("topright", c("Complete data", "Observed data"), pt.bg = COL[c("complete", "observed")], pch = 21, bty = "n", cex = 1)
  box(); dev.off()

  # ---- Figure 2: hidden skew changes the test, and the test chosen after a pass errs more ---------------
  mk <- read_table(run_dir, "pretest_masking.csv"); cd <- read_table(run_dir, "pretest_conditional.csv")
  mk$mechanism <- factor(mk$mechanism, levels = names(MECH_NAMES)); mk <- mk[order(mk$mechanism), ]
  f <- open_png(out, "Figure_2_hidden_skew", 3.9); made <- c(made, f)
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 1.6, 0.6), mgp = c(2.4, 0.6, 0), tcl = -0.3, las = 1)
  m <- rbind(mk$sw_pass_complete, mk$sw_pass_observed)
  bp <- barplot(m, beside = TRUE, col = COL[c("complete", "observed")], border = "white", ylim = c(0, max(m) * 1.3),
                names.arg = rep("", ncol(m)), ylab = "Share of replications passing")
  mtext(c("MCAR", "Group-\ndependent", "Upper tail\nlost", "Lower tail\nlost")[seq_len(ncol(m))], side = 1, at = colMeans(bp), line = 0.4, padj = 1, cex = 0.85)
  legend("topleft", c("Complete data", "Observed data"), fill = COL[c("complete", "observed")], border = "white", bty = "n", cex = 0.9)
  mtext("A", side = 3, adj = 0, line = 0.3, font = 2)
  w <- cd[grepl("^Cond: Welch t", cd$label), ]; w$shape <- factor(w$shape, levels = names(SHAPE_NAMES))
  pass <- w[grepl("pass$", w$label), ]; fail <- w[grepl("fail$", w$label), ]
  pass <- pass[order(pass$shape), ]; fail <- fail[order(fail$shape), ]
  m2 <- rbind(fail$conditional_error, pass$conditional_error)
  bp2 <- barplot(m2, beside = TRUE, col = COL[c("fail", "pass")], border = "white", ylim = c(0, max(m2) * 1.3),
          names.arg = rep("", ncol(m2)), ylab = expression("Rejection rate of Welch " * italic(t)))   # APA: italic t
  mtext(sub("-", "-\n", sub(" ", "\n", SHAPE_NAMES[as.character(pass$shape)])), side = 1, at = colMeans(bp2), line = 0.4, padj = 1, cex = 0.85)
  abline(h = 0.05, lty = 3, col = "gray40")
  legend("topleft", c("After a failed pre-test", "After a passed pre-test"), fill = COL[c("fail", "pass")], border = "white", bty = "n", cex = 0.9)
  mtext("B", side = 3, adj = 0, line = 0.3, font = 2)
  dev.off()

  # ---- Figure 3: harmful share, Welch t alone vs pre-test then Welch t or Mann-Whitney ------------------
  e <- read_table(run_dir, "pretest_error.csv")
  get <- function(g, p, col = "harmful_share") e[[col]][e$group == g & e$procedure == p]
  # unequal spreads = all designs with equal means minus those with identical distributions (exact, from the counts)
  share <- function(g, p) {
    if (g != "unequal spreads" || "unequal spreads" %in% e$group) return(get(g, p))   # written by summarize_pretest.R 0.9
    na <- get("all", p, "designs"); ni <- get("identical distributions", p, "designs")
    (get("all", p) * na - get("identical distributions", p) * ni) / (na - ni) }
  groups <- c("identical distributions", "unequal spreads"); gnames <- c("Identical distributions", "Equal means, unequal spreads")
  m3 <- sapply(groups, function(g) c(share(g, "Complete: Welch t"), share(g, "Complete: per group, Welch t or MW"),
                                     share(g, "Welch t"), share(g, "Two-stage: per group, Welch t or MW")))
  f <- open_png(out, "Figure_3_harmful_share", 3.6); made <- c(made, f)
  par(mar = c(3.2, 4, 0.8, 0.6), mgp = c(2.4, 0.6, 0), tcl = -0.3, las = 1)
  light <- function(x) { v <- col2rgb(x) * 0.4 + 255 * 0.6; rgb(v[1], v[2], v[3], maxColorValue = 255) }
  cols <- c(light(COL[["welch"]]), light(COL[["two_stage"]]), COL[["welch"]], COL[["two_stage"]])
  bp3 <- barplot(m3, beside = TRUE, col = cols, border = "white", ylim = c(0, max(m3) * 1.6), names.arg = gnames, cex.names = 1,
          ylab = "Share of harmful designs")
  text(bp3, m3, sprintf("%.1f%%", 100 * m3), pos = 3, cex = 0.8, offset = 0.3)
  legend("top", expression("Welch " * italic(t) * ", complete data", "Welch " * italic(t) * ", observed data", "Pre-test, complete data", "Pre-test, observed data"),
         fill = cols[c(1, 3, 2, 4)], border = "white", bty = "n", cex = 0.9, ncol = 2)
  dev.off()

  cat(sprintf("%d figures written to %s\n", length(made), normalizePath(out))); for (f in made) cat("  ", basename(f), "\n")
  invisible(made)
}
paper_figures <- make_figures(RUN_DIR)
