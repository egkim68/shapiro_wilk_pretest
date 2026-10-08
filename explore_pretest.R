# EXPLORER: the Shapiro-Wilk pre-test when missing data shape what is left
# Base R only. Writes CSV files only. Run: source("explore_pretest.R") in R (pasting also works).
# Finished configurations are checkpointed, so an interrupted run resumes when you run it again.
#
# QUESTION: researchers often test normality with Shapiro-Wilk and then choose a t test or the
# Mann-Whitney test. Missing data change the sample the pre-test sees. Selective missingness can
# make normal data look non-normal (larger or smaller values go missing), and random loss can hide
# real non-normality (smaller groups, less power). This explorer measures
#   (1) how often Shapiro-Wilk rejects on the observed data, next to the same test on the complete
#       data before deletion (so the change caused by missingness is visible design by design);
#   (2) the Type I error and power of the tests alone and of four two-stage procedures that follow
#       the pre-test (per group or on pooled residuals; Student t or Welch t when it passes,
#       Mann-Whitney when it fails);
#   (3) the error of the t tests conditional on the pre-test passing or failing.
#
# METHOD: every configuration is DRAWN AT RANDOM from wide ranges (see sample_parameters).
# TRAIN configurations are used to find rules; a separate TEST set of fresh configurations checks
# them (analyze_explorer.R). Two groups; parameters (p_):
#   p_n_total 10-200 (log-uniform) split by p_ratio_n = n1/n2 (0.5-2); at least 5 planned per group
#   p_sd_ratio = SD2/SD1: 1 in 50% of designs, otherwise 0.33-3 (log-uniform)
#   p_mechanism (25% / 25% / 50%, so that selective loss on skewed data is well covered):
#     MCAR (same rate in both groups), MAR_group (each group its own rate), MNAR (the
#     probability of deletion is logistic in the raw value, strength p_strength, sign p_direction,
#     one rule for both groups with its intercept calibrated on the pooled sample so the average rate
#     equals p_rate1)
#   p_rate1, p_rate2: missing rates 0-0.5 (equal under MCAR and MNAR)
#   p_shape: normal, mild skew (lognormal sigma 0.4), skewed (lognormal sigma 0.8), heavy-tailed (t, 3 df),
#     each standardized to mean 0 and SD 1
#   p_outlier: none in 60% of designs, otherwise 1% to 5% of cells shifted by p_outlier_size (4-8) of
#     the group's own SD
#   p_effect: 0 in 60% of designs (equal population means), otherwise a shift of 0.2-1.0 group-1 SDs
# Under equal means with unequal spread or skew, the Mann-Whitney test addresses a different null
# (stochastic equality), so its rejections there are rejection rates, not Type I errors.
# A replication with fewer than 4 observations left in a group counts as a failure.
#
# LABELS (one row per configuration and label; rejection_rate = share of replications with p < alpha):
#   "SW per group (complete)", "SW residuals (complete)"   pre-test on the data BEFORE deletion
#   "SW per group", "SW residuals"                         pre-test on the OBSERVED data
#       per group = non-normal if either group's Shapiro-Wilk p < alpha (the common practice);
#       residuals = one Shapiro-Wilk test on the pooled within-group residuals
#   "SW decision flipped", "SW flipped to pass", "SW flipped to fail"   share of replications in which the
#       per-group pre-test decides differently on the observed data than on the complete data
#   "Student t", "Welch t", "Mann-Whitney"                 the tests alone
#   "Two-stage: <which pre-test>, <t test> or MW"          the pre-test decides: t test if it passes, MW if not
#   "Complete: Welch t", "Complete: per group, Welch t or MW"   the same on the data BEFORE deletion
#   "Cond: <t test> | <which pre-test> pass/fail"          the t test only in replications where the
#       pre-test passed (or failed); rejection_rate is then CONDITIONAL on that outcome, and
#       `successful` counts the replications in that condition
# SW and Cond labels are descriptive and are excluded from the Type I error models in analyze_explorer.R.

# "pilot": quick check; "medium": the exploratory phase (400 configurations, 100 replications);
# "full": the confirmatory phase (2,000 configurations, 1,000 replications; about 5 hours, resumable). Set the environment variable
# SIM_MODE (the launchers do this) or change the text below.
MODE <- Sys.getenv("SIM_MODE", unset = "medium")
N_TRAIN <- switch(MODE, pilot = 60L, medium = 300L, full = 1500L)
N_TEST <- switch(MODE, pilot = 20L, medium = 100L, full = 500L)
N_REPLICATIONS <- switch(MODE, pilot = 60L, medium = 100L, full = 1000L)   # full: 1,000 replications, SE about 0.008 near 0.075
# The confirmatory FULL run uses a different, declared seed (exploratory seed + 730000).
BASE_SEED <- if (MODE == "full") 2041L + 730000L else 2041L
ALPHA_LEVEL <- 0.05
PRETEST_ALPHA <- 0.05
RESUME <- TRUE
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer1")
OUTPUT_DIR <- file.path(OUTPUT_ROOT, paste0("pretest_", MODE))
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")


# ---- building blocks -----------------------------------------------------------------------------
calibrated_missing <- function(eta, rate) {
  n <- length(eta)
  if (rate <= 0) return(rep(FALSE, n))
  intercept <- uniroot(function(a) mean(plogis(a + eta)) - rate,
                       c(-max(eta) - 40, -min(eta) + 40), tol = 1e-8)$root
  runif(n) < plogis(intercept + eta)
}
contaminate <- function(x, prob, size) {
  if (prob <= 0) return(x)
  hit <- runif(length(x)) < prob
  size <- rep_len(size, length(x))
  x[hit] <- x[hit] + size[hit] * sample(c(-1, 1), sum(hit), replace = TRUE)
  x
}
std_lognormal <- function(n, s) (rlnorm(n, 0, s) - exp(s^2 / 2)) / sqrt((exp(s^2) - 1) * exp(s^2))
draw_shape <- function(n, shape) {
  switch(shape, normal = rnorm(n), mild_skew = std_lognormal(n, 0.4), skewed = std_lognormal(n, 0.8),
         heavy = rt(n, 3) / sqrt(3), stop("Unknown shape"))
}
POP_SKEW <- c(normal = 0, mild_skew = (exp(0.16) + 2) * sqrt(exp(0.16) - 1),
              skewed = (exp(0.64) + 2) * sqrt(exp(0.64) - 1), heavy = 0)
sw_p <- function(x) if (length(x) >= 3 && length(x) <= 5000 && diff(range(x)) > 0) shapiro.test(x)$p.value else NA_real_
row_out <- function(label, p_value, estimate, n_used)
  data.frame(label = label, p_value = p_value, estimate = estimate, n_used = n_used,
             status = "Success", stringsAsFactors = FALSE)
mean_safe <- function(x) if (length(x)) mean(x) else NA_real_


# ---- configuration sampler and design ------------------------------------------------------------
TESTS <- c("Student t", "Welch t", "Mann-Whitney")
TWO_STAGE <- c("Two-stage: per group, Student t or MW", "Two-stage: per group, Welch t or MW",
               "Two-stage: residuals, Student t or MW", "Two-stage: residuals, Welch t or MW")
SW_LABELS <- c("SW per group (complete)", "SW residuals (complete)", "SW per group", "SW residuals",
               "SW decision flipped", "SW flipped to pass", "SW flipped to fail")
# the same Welch test and two-stage procedure on the complete data, so the part added by missingness is visible
COMPLETE_LABELS <- c("Complete: Welch t", "Complete: per group, Welch t or MW")
ANALYSIS_LABELS <- c(SW_LABELS, TESTS, TWO_STAGE, COMPLETE_LABELS)   # Cond labels appear only when their condition occurs
expected_labels <- function(p) ANALYSIS_LABELS

sample_parameters <- function(n) {
  mech <- sample(c("MCAR", "MAR_group", "MNAR"), n, replace = TRUE, prob = c(0.25, 0.25, 0.5))
  r1 <- round(runif(n, 0, 0.5), 3); r2 <- round(runif(n, 0, 0.5), 3)
  data.frame(
    p_n_total = round(exp(runif(n, log(10), log(200)))),
    p_ratio_n = round(exp(runif(n, log(0.5), log(2))), 3),
    p_sd_ratio = round(ifelse(runif(n) < 0.5, 1, exp(runif(n, log(1 / 3), log(3)))), 3),
    p_mechanism = mech,
    p_rate1 = r1,
    p_rate2 = ifelse(mech == "MAR_group", r2, r1),
    p_strength = round(runif(n, 0.5, 2.5), 2),
    p_direction = sample(c(-1, 1), n, replace = TRUE),
    p_outlier = round(ifelse(runif(n) < 0.6, 0, runif(n, 0.01, 0.05)), 3),
    p_outlier_size = round(runif(n, 4, 8), 1),
    p_shape = sample(c("normal", "mild_skew", "skewed", "heavy"), n, replace = TRUE),
    p_effect = ifelse(runif(n) < 0.6, 0, round(runif(n, 0.2, 1.0), 2)),
    stringsAsFactors = FALSE)
}
make_design <- function(p) {
  n1 <- max(5L, as.integer(round(p$p_n_total * p$p_ratio_n / (1 + p$p_ratio_n))))
  n2 <- max(5L, as.integer(p$p_n_total - n1))
  list(n1 = n1, n2 = n2)
}
describe <- function(p, design) {
  e1 <- design$n1 * (1 - p$p_rate1); e2 <- design$n2 * (1 - p$p_rate2)
  data.frame(d_n1 = design$n1, d_n2 = design$n2, d_expected_n1 = e1, d_expected_n2 = e2,
             d_min_expected_n = min(e1, e2), d_log_sd_ratio = log(p$p_sd_ratio),
             d_pop_skew = unname(POP_SKEW[p$p_shape]))
}

rep_function <- function(p, design) {
  y1 <- contaminate(draw_shape(design$n1, p$p_shape), p$p_outlier, p$p_outlier_size)
  y2 <- contaminate(p$p_effect + p$p_sd_ratio * draw_shape(design$n2, p$p_shape),
                    p$p_outlier, p$p_outlier_size * p$p_sd_ratio)
  # MNAR: one deletion rule for both groups, its intercept calibrated once on the pooled sample (as in the
  # ANOVA toolkit). Calibrating separately within each group would force each group to lose exactly its
  # share, which makes the kept values dependent and inflates the t tests' error even with identical groups.
  gone <- function(y, rate) if (rate <= 0) rep(FALSE, length(y)) else runif(length(y)) < rate
  # the pre-test on the complete data, before any value is deleted
  sw_group_full <- min(sw_p(y1), sw_p(y2))
  sw_resid_full <- sw_p(c(y1 - mean(y1), y2 - mean(y2)))
  p_w_full <- t.test(y1, y2)$p.value
  p_mw_full <- suppressWarnings(wilcox.test(y1, y2, exact = FALSE)$p.value)
  d_full <- mean(y2) - mean(y1)
  if (p$p_mechanism == "MNAR") {
    g <- if (p$p_rate1 <= 0) rep(FALSE, length(y1) + length(y2)) else calibrated_missing(p$p_direction * p$p_strength * c(y1, y2), p$p_rate1)
    g1 <- g[seq_along(y1)]; g2 <- g[-seq_along(y1)]
  } else { g1 <- gone(y1, p$p_rate1); g2 <- gone(y2, p$p_rate2) }
  a <- y1[!g1]; b <- y2[!g2]
  if (length(a) < 4 || length(b) < 4) stop("Too few observations")
  sw_group <- min(sw_p(a), sw_p(b))
  sw_resid <- sw_p(c(a - mean(a), b - mean(b)))
  if (is.na(sw_group) || is.na(sw_resid)) stop("Shapiro-Wilk not computable")
  d <- mean(b) - mean(a); used <- length(a) + length(b)
  p_t <- t.test(a, b, var.equal = TRUE)$p.value
  p_w <- t.test(a, b)$p.value
  p_mw <- suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value)
  pass_g <- sw_group >= PRETEST_ALPHA; pass_r <- sw_resid >= PRETEST_ALPHA
  out <- rbind(
    row_out("SW per group (complete)", sw_group_full, NA_real_, length(y1) + length(y2)),
    row_out("SW residuals (complete)", sw_resid_full, NA_real_, length(y1) + length(y2)),
    row_out("SW per group", sw_group, NA_real_, used),
    row_out("SW residuals", sw_resid, NA_real_, used),
    # decision flips: the per-group pre-test decides differently on the observed data than on the complete
    # data (coded so that rejection_rate = share of replications with that flip)
    row_out("SW decision flipped", if (pass_g != (sw_group_full >= PRETEST_ALPHA)) 0 else 1, NA_real_, used),
    row_out("SW flipped to pass", if (pass_g && sw_group_full < PRETEST_ALPHA) 0 else 1, NA_real_, used),
    row_out("SW flipped to fail", if (!pass_g && sw_group_full >= PRETEST_ALPHA) 0 else 1, NA_real_, used),
    row_out("Student t", p_t, d, used), row_out("Welch t", p_w, d, used), row_out("Mann-Whitney", p_mw, d, used),
    row_out(TWO_STAGE[1], if (pass_g) p_t else p_mw, d, used), row_out(TWO_STAGE[2], if (pass_g) p_w else p_mw, d, used),
    row_out(TWO_STAGE[3], if (pass_r) p_t else p_mw, d, used), row_out(TWO_STAGE[4], if (pass_r) p_w else p_mw, d, used),
    row_out(COMPLETE_LABELS[1], p_w_full, d_full, length(y1) + length(y2)),
    row_out(COMPLETE_LABELS[2], if (sw_group_full >= PRETEST_ALPHA) p_w_full else p_mw_full, d_full, length(y1) + length(y2)),
    row_out(sprintf("Cond: Student t | per group %s", if (pass_g) "pass" else "fail"), p_t, d, used),
    row_out(sprintf("Cond: Welch t | per group %s", if (pass_g) "pass" else "fail"), p_w, d, used))
  out
}


# ---- running and summarising (shared with the other explorers) -------------------------------------
aggregate_reps <- function(rows, expected = ANALYSIS_LABELS) {
  failed <- sum(rows$label == "FAILED")
  rows <- rows[rows$label != "FAILED", , drop = FALSE]
  pieces <- lapply(split(rows, rows$label), function(d) {
    ok <- d$status == "Success" & !is.na(d$p_value); s <- d[ok, ]; m <- nrow(s)
    rate <- mean_safe(s$p_value < ALPHA_LEVEL)
    data.frame(label = d$label[1], successful = m,
               failure_rate = 1 - m / N_REPLICATIONS, rejection_rate = rate,
               rejection_mcse = if (m > 0) sqrt(rate * (1 - rate) / m) else NA_real_,
               mean_estimate = mean_safe(s$estimate[!is.na(s$estimate)]), mean_n_used = mean_safe(s$n_used),
               stringsAsFactors = FALSE)
  })
  for (label in setdiff(expected, names(pieces)))
    pieces[[label]] <- data.frame(label = label, successful = 0L, failure_rate = 1,
      rejection_rate = NA_real_, rejection_mcse = NA_real_, mean_estimate = NA_real_,
      mean_n_used = NA_real_, stringsAsFactors = FALSE)
  out <- do.call(rbind, pieces); rownames(out) <- NULL
  out$failed_replications <- failed
  out
}
sample_configs <- function(n, seed, offset) {
  set.seed(seed)
  configs <- sample_parameters(n)
  cbind(data.frame(config_id = seq_len(n), config_seed = BASE_SEED + offset + 100000000L + seq_len(n)), configs)
}
run_set <- function(set_name, configs, offset) {
  results <- vector("list", nrow(configs))
  for (i in seq_len(nrow(configs))) {
    p <- configs[i, ]
    checkpoint <- file.path(CHECKPOINT_DIR, sprintf("%s_%05d.rds", set_name, i))
    if (RESUME && file.exists(checkpoint)) res <- readRDS(checkpoint) else {
      design <- make_design(p)
      rows <- lapply(seq_len(N_REPLICATIONS), function(b) {
        set.seed(BASE_SEED + offset + p$config_id * 1000L + b)
        tryCatch(rep_function(p, design), error = function(e) data.frame(label = "FAILED",
          p_value = NA_real_, estimate = NA_real_, n_used = NA_real_,
          status = conditionMessage(e), stringsAsFactors = FALSE))
      })
      res <- aggregate_reps(do.call(rbind, rows), expected_labels(p))
      desc <- describe(p, design)
      res <- cbind(p[rep(1L, nrow(res)), , drop = FALSE], desc[rep(1L, nrow(res)), , drop = FALSE], res)
      rownames(res) <- NULL
      temporary <- paste0(checkpoint, ".tmp")
      saveRDS(res, temporary)
      if (file.exists(checkpoint)) unlink(checkpoint)
      if (!file.rename(temporary, checkpoint)) stop("Checkpoint save failed")
    }
    results[[i]] <- res
    if (i %% 25 == 0 || i == nrow(configs)) cat(sprintf("[%s] configuration %d/%d\n", set_name, i, nrow(configs)))
  }
  out <- do.call(rbind, results); rownames(out) <- NULL
  out$set <- set_name
  out
}

CODE_VERSION <- "0.5.0"
run_explorer <- function() {
  stopifnot(MODE %in% c("pilot", "medium", "full"))
  cat(sprintf("Shapiro-Wilk pre-test explorer, %s mode: %d train + %d test configurations, %d replications each.\nResults go to: %s\n",
              toupper(MODE), N_TRAIN, N_TEST, N_REPLICATIONS, OUTPUT_DIR))
  dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(CHECKPOINT_DIR))
    stop("Cannot create the output folder ", OUTPUT_DIR, ". Check that the drive exists or set SIM_OUTPUT_ROOT.", call. = FALSE)
  train <- sample_configs(N_TRAIN, BASE_SEED + 11L, 0L)
  test <- sample_configs(N_TEST, BASE_SEED + 22L, 500000000L)
  config <- list(version = CODE_VERSION, mode = MODE, train = N_TRAIN, test = N_TEST, replications = N_REPLICATIONS,
                 seed = BASE_SEED, alpha = ALPHA_LEVEL, pretest_alpha = PRETEST_ALPHA, train_configs = train, test_configs = test)
  config_file <- file.path(OUTPUT_DIR, "configuration.rds")
  if (file.exists(config_file) && RESUME) {
    old <- readRDS(config_file)
    if (!identical(old, config)) {
      archived <- paste0(OUTPUT_DIR, "_archived_", format(Sys.time(), "%Y%m%d%H%M%S"))
      cat(sprintf("The results in %s were made with other settings or code version %s.\nArchiving that folder as %s and starting a fresh run.\n\n",
                  OUTPUT_DIR, as.character(old$version), archived))
      if (!file.rename(OUTPUT_DIR, archived)) stop("Could not rename ", OUTPUT_DIR, ". Move it by hand and run again.", call. = FALSE)
      dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
      saveRDS(config, config_file)
    }
  } else saveRDS(config, config_file)
  writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info.txt"))
  out <- list()
  for (set_name in c("train", "test")) {
    configs <- if (set_name == "train") train else test
    offset <- if (set_name == "train") 0L else 500000000L
    out[[set_name]] <- run_set(set_name, configs, offset)
    write.csv(out[[set_name]], file.path(OUTPUT_DIR, sprintf("explore_pretest_%s.csv", set_name)), row.names = FALSE)
  }
  cat("\nFinished. Results:", normalizePath(OUTPUT_DIR), "\n")
  invisible(out)
}
explorer_results <- run_explorer()
