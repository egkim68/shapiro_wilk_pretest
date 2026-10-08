# SUMMARIZE THE PRE-TEST STUDY. Base R only. Changes no result and runs no simulation: rerun any time.
# Reads explore_pretest_train.csv and explore_pretest_test.csv of one run (the mode in SIM_MODE, default
# medium) and writes these tables to <run>/analysis, printing a short version of each:
#   pretest_design_summary.csv     designs generated/analysed, failed replications, median planned sample size
#   pretest_sw_rates.csv       Q1  Shapiro-Wilk rejection on the complete data and on the observed data,
#                                  by population shape and missingness mechanism (all designs)
#   pretest_error.csv          Q2  mean Type I error and share of harmful designs of each test and each
#                                  two-stage procedure (designs with equal means): all, identical
#                                  distributions, unequal spreads, and by shape and mechanism
#   pretest_conditional.csv    Q2  error of the t tests when the per-group pre-test passed or failed
#   pretest_masking.csv        Q2  skewed populations: does hidden skew lead to the t test, and at what cost
#   pretest_power.csv          Q2  mean power of each test and procedure (designs with a true difference)
#   pretest_diagnostic.csv     Q3  can an observed Shapiro-Wilk rejection rate tell selective (MNAR)
#                                  missingness from random missingness? held-out AUC, normal populations
# Graphs are drawn by a separate script, so they can be redrawn later without rerunning anything.
# Use:  Rscript summarize_pretest.R            or   Rscript summarize_pretest.R <run folder>

args <- commandArgs(trailingOnly = TRUE)
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer1")
MODE <- Sys.getenv("SIM_MODE", unset = "medium")
RUN_DIR <- if (length(args) >= 1) args[1] else file.path(OUTPUT_ROOT, paste0("pretest_", MODE))
HARM_LEVEL <- as.numeric(Sys.getenv("SIM_HARM_LEVEL", unset = "0.075"))
MIN_SUCCESS <- 20

summarize_pretest <- function(run_dir) {
  f <- file.path(run_dir, c("explore_pretest_train.csv", "explore_pretest_test.csv"))
  if (!all(file.exists(f))) stop("No finished pre-test run in ", run_dir, ". Run explore_pretest.R first.", call. = FALSE)
  d <- rbind(read.csv(f[1], stringsAsFactors = FALSE), read.csv(f[2], stringsAsFactors = FALSE))
  d$key <- paste(d$set, d$config_id)
  out_dir <- file.path(run_dir, "analysis"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  save <- function(x, name) { write.csv(x, file.path(out_dir, name), row.names = FALSE); x }
  val <- function(label, col = "rejection_rate", data = d) { z <- data[data$label == label & data$successful >= MIN_SUCCESS, ]; setNames(z[[col]], z$key) }
  designs <- d[!duplicated(d$key), ]
  cat(sprintf("Pre-test study: %d designs (%d with equal means) in %s; harmful above %.3f\n",
              nrow(designs), sum(designs$p_effect == 0), normalizePath(run_dir), HARM_LEVEL))

  # ---- design summary: counts, failed replications, planned sample size (reported in the Methods) ------
  w <- d[d$label == "Welch t", ]
  analysed <- w$key[w$successful >= MIN_SUCCESS]; an <- designs[designs$key %in% analysed, ]
  ds <- data.frame(designs_generated = nrow(designs), designs_analysed = nrow(an), designs_set_aside = nrow(designs) - nrow(an),
                   equal_means = sum(an$p_effect == 0), identical_distributions = sum(an$p_effect == 0 & an$p_sd_ratio == 1),
                   unequal_spreads = sum(an$p_effect == 0 & an$p_sd_ratio != 1), true_difference = sum(an$p_effect != 0),
                   failed_replication_share = sum(w$failed_replications) / sum(w$failed_replications + w$successful),
                   median_planned_n = median(designs$p_n_total), designs_without_outliers = sum(an$p_outlier == 0))
  ds <- save(ds, "pretest_design_summary.csv")
  cat(sprintf("Designs: %d generated, %d analysed (%d set aside); equal means %d (identical %d, unequal spreads %d); true difference %d\n",
              ds$designs_generated, ds$designs_analysed, ds$designs_set_aside, ds$equal_means, ds$identical_distributions, ds$unequal_spreads, ds$true_difference))
  cat(sprintf("Failed replications %.1f%%; median planned sample size %.0f; designs without outliers %d\n",
              100 * ds$failed_replication_share, ds$median_planned_n, ds$designs_without_outliers))

  # ---- Q1: what missingness does to the Shapiro-Wilk verdict --------------------------------------
  # MNAR is split by direction: larger values lost (direction +1) or smaller values lost (-1)
  mech <- ifelse(designs$p_mechanism != "MNAR", designs$p_mechanism, ifelse(designs$p_direction > 0, "MNAR_large_lost", "MNAR_small_lost"))
  sw <- data.frame(key = designs$key, shape = designs$p_shape, mechanism = mech,
                   outliers = designs$p_outlier > 0, strength = designs$p_strength, stringsAsFactors = FALSE)
  for (l in c("SW per group (complete)", "SW per group", "SW residuals (complete)", "SW residuals")) sw[[l]] <- val(l)[sw$key]
  q1 <- do.call(rbind, lapply(split(sw[!sw$outliers, ], list(sw$shape[!sw$outliers], sw$mechanism[!sw$outliers]), drop = TRUE), function(x)
    data.frame(shape = x$shape[1], mechanism = x$mechanism[1], designs = nrow(x),
               group_complete = mean(x[["SW per group (complete)"]], na.rm = TRUE), group_observed = mean(x[["SW per group"]], na.rm = TRUE),
               resid_complete = mean(x[["SW residuals (complete)"]], na.rm = TRUE), resid_observed = mean(x[["SW residuals"]], na.rm = TRUE))))
  q1$group_change <- q1$group_observed - q1$group_complete; q1$resid_change <- q1$resid_observed - q1$resid_complete
  for (l in c("SW decision flipped", "SW flipped to pass", "SW flipped to fail")) sw[[l]] <- val(l)[sw$key]
  q1$flipped <- sapply(seq_len(nrow(q1)), function(i) mean(sw[["SW decision flipped"]][!sw$outliers & sw$shape == q1$shape[i] & sw$mechanism == q1$mechanism[i]], na.rm = TRUE))
  q1$flipped_to_pass <- sapply(seq_len(nrow(q1)), function(i) mean(sw[["SW flipped to pass"]][!sw$outliers & sw$shape == q1$shape[i] & sw$mechanism == q1$mechanism[i]], na.rm = TRUE))
  q1 <- save(q1[order(q1$shape, q1$mechanism), ], "pretest_sw_rates.csv")
  cat("\nQ1  Shapiro-Wilk rejection rate, complete data -> observed data (designs without outliers)\n")
  for (i in seq_len(nrow(q1))) cat(sprintf("  %-10s %-16s n=%-4d per group %.2f -> %.2f (%+.2f) | residuals %.2f -> %.2f (%+.2f) | decision flipped %.2f (to pass %.2f)\n",
    q1$shape[i], q1$mechanism[i], q1$designs[i], q1$group_complete[i], q1$group_observed[i], q1$group_change[i],
    q1$resid_complete[i], q1$resid_observed[i], q1$resid_change[i], q1$flipped[i], q1$flipped_to_pass[i]))
  fl <- sw[!sw$outliers, ]; cat(sprintf("  all designs without outliers: the per-group decision flipped in %.1f%% of replications (to pass %.1f%%, to fail %.1f%%)\n",
    100 * mean(fl[["SW decision flipped"]], na.rm = TRUE), 100 * mean(fl[["SW flipped to pass"]], na.rm = TRUE), 100 * mean(fl[["SW flipped to fail"]], na.rm = TRUE)))

  # ---- Q2: error, conditional error and power ------------------------------------------------------
  procs <- c("Student t", "Welch t", "Mann-Whitney", "Two-stage: per group, Student t or MW", "Two-stage: per group, Welch t or MW",
             "Two-stage: residuals, Student t or MW", "Two-stage: residuals, Welch t or MW")
  procs <- c(procs, "Complete: Welch t", "Complete: per group, Welch t or MW")
  nul <- d[d$p_effect == 0 & d$successful >= MIN_SUCCESS, ]
  # identical distributions: equal means AND equal spreads, so the t tests and Mann-Whitney share the same null
  same <- nul[nul$p_sd_ratio == 1, ]
  # unequal spreads: equal means, different spreads (the t tests address a true null, Mann-Whitney does not)
  diff_sd <- nul[nul$p_sd_ratio != 1, ]
  err_row <- function(x, group) do.call(rbind, lapply(procs, function(l) { r <- x$rejection_rate[x$label == l]
    data.frame(group = group, procedure = l, designs = length(r), mean_error = mean(r), harmful_share = mean(r > HARM_LEVEL),
               above_0.10 = mean(r > 0.10), stringsAsFactors = FALSE) }))
  q2 <- rbind(err_row(nul, "all"), err_row(same, "identical distributions"), err_row(diff_sd, "unequal spreads"),
              do.call(rbind, lapply(split(same, same$p_shape), function(x) err_row(x, paste("identical,", x$p_shape[1])))),
              do.call(rbind, lapply(split(nul, nul$p_shape), function(x) err_row(x, x$p_shape[1]))),
              do.call(rbind, lapply(split(nul, nul$p_mechanism), function(x) err_row(x, x$p_mechanism[1]))))
  q2 <- save(q2, "pretest_error.csv")
  cat(sprintf("\nQ2  Type I error with equal means (Mann-Whitney: rejection rate), all designs; harmful = above %.3f\n", HARM_LEVEL))
  for (g in c("all", "identical distributions", "unequal spreads")) { a <- q2[q2$group == g, ]; cat(sprintf(" [%s, %d designs]\n", g, a$designs[1]))
    for (i in seq_len(nrow(a))) cat(sprintf("  %-40s mean %.3f, harmful %4.1f%%, above 0.10 %4.1f%%\n", a$procedure[i], a$mean_error[i], 100 * a$harmful_share[i], 100 * a$above_0.10[i])) }

  # hidden skew -> the pre-test passes -> a t test is used on skewed data: by mechanism, identical distributions
  sk <- same[same$p_shape %in% c("skewed", "mild_skew"), ]
  sk$mech <- ifelse(sk$p_mechanism != "MNAR", sk$p_mechanism, ifelse(sk$p_direction > 0, "MNAR_large_lost", "MNAR_small_lost"))
  q2m <- do.call(rbind, lapply(split(sk, sk$mech), function(x) { k <- unique(x$key); g <- function(l) mean(x$rejection_rate[x$label == l])
    data.frame(mechanism = x$mech[1], designs = length(k), sw_pass_complete = 1 - g("SW per group (complete)"), sw_pass_observed = 1 - g("SW per group"),
               welch_complete = g("Complete: Welch t"), welch_observed = g("Welch t"),
               two_stage_complete = g("Complete: per group, Welch t or MW"), two_stage_observed = g("Two-stage: per group, Welch t or MW"),
               stringsAsFactors = FALSE) }))
  q2m <- save(q2m, "pretest_masking.csv")
  cat("\nQ2  skewed populations with identical distributions: pre-test pass rate and error, complete -> observed\n")
  for (i in seq_len(nrow(q2m))) cat(sprintf("  %-16s n=%-3d SW passes %.2f -> %.2f | Welch %.3f -> %.3f | two-stage %.3f -> %.3f\n", q2m$mechanism[i], q2m$designs[i],
    q2m$sw_pass_complete[i], q2m$sw_pass_observed[i], q2m$welch_complete[i], q2m$welch_observed[i], q2m$two_stage_complete[i], q2m$two_stage_observed[i]))

  cond <- d[d$p_effect == 0 & grepl("^Cond: ", d$label) & d$successful > 0, ]
  q2c <- do.call(rbind, lapply(split(cond, list(cond$label, cond$p_shape), drop = TRUE), function(x)
    data.frame(label = x$label[1], shape = x$p_shape[1], designs = nrow(x), replications = sum(x$successful),
               conditional_error = sum(x$rejection_rate * x$successful) / sum(x$successful), stringsAsFactors = FALSE)))
  q2c <- save(q2c[order(q2c$label, q2c$shape), ], "pretest_conditional.csv")
  cat("\nQ2  error of the t tests conditional on the per-group pre-test (equal means; pooled over replications)\n")
  for (i in seq_len(nrow(q2c))) cat(sprintf("  %-36s %-10s %.3f  (%d replications)\n", q2c$label[i], q2c$shape[i], q2c$conditional_error[i], q2c$replications[i]))

  alt <- d[d$p_effect != 0 & d$successful >= MIN_SUCCESS, ]
  q2p <- save(do.call(rbind, lapply(procs, function(l) data.frame(procedure = l, designs = sum(alt$label == l),
    mean_power = mean(alt$rejection_rate[alt$label == l]), stringsAsFactors = FALSE))), "pretest_power.csv")
  cat("\nQ2  power (raw, not size-adjusted), designs with a true difference\n")
  for (i in seq_len(nrow(q2p))) cat(sprintf("  %-40s %.3f\n", q2p$procedure[i], q2p$mean_power[i]))

  # ---- Q3: is a Shapiro-Wilk rejection a warning of selective missingness? ------------------------
  auc <- function(score, y) { ok <- !is.na(score); score <- score[ok]; y <- y[ok]; r <- rank(score); n1 <- sum(y); n0 <- sum(!y)
    if (!n1 || !n0) NA_real_ else (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
  q3 <- do.call(rbind, lapply(list(list("normal populations, no outliers", sw$shape == "normal" & !sw$outliers),
                                   list("all shapes, no outliers", !sw$outliers), list("all designs", rep(TRUE, nrow(sw)))), function(s) {
    x <- sw[s[[2]], ]; y <- grepl("^MNAR", x$mechanism); te <- grepl("^test", x$key)
    data.frame(subset = s[[1]], designs = nrow(x), mnar = sum(y),
               auc_per_group_all = auc(x[["SW per group"]], y), auc_per_group_test = auc(x[["SW per group"]][te], y[te]),
               auc_residuals_all = auc(x[["SW residuals"]], y), auc_residuals_test = auc(x[["SW residuals"]][te], y[te]),
               stringsAsFactors = FALSE) }))
  q3 <- save(q3, "pretest_diagnostic.csv")
  cat("\nQ3  does the observed Shapiro-Wilk rejection rate separate MNAR from MCAR/MAR designs? (AUC; 0.5 = no)\n")
  for (i in seq_len(nrow(q3))) cat(sprintf("  %-32s n=%-4d MNAR=%-4d per group %.2f (held-out %.2f) | residuals %.2f (held-out %.2f)\n",
    q3$subset[i], q3$designs[i], q3$mnar[i], q3$auc_per_group_all[i], q3$auc_per_group_test[i], q3$auc_residuals_all[i], q3$auc_residuals_test[i]))
  cat("\nTables written to", normalizePath(out_dir), "\n")
  invisible(list(sw = q1, error = q2, conditional = q2c, power = q2p, diagnostic = q3))
}
pretest_summary <- summarize_pretest(RUN_DIR)
