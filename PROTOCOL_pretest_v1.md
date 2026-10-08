# Protocol: frozen hypotheses for the confirmatory run

**Study.** When Missing Data Change the Normality Verdict: Shapiro-Wilk Pre-Testing and Test Choice in Two-Group Comparisons (Eungi Kim)

**Frozen on 1 October 2026 (KST), before the confirmatory run began.** No confirmatory result existed when this document was written. The hypotheses come from the exploratory run (toolkit 0.6, 300 training + 100 held-out designs, 100 replications, seed 2041).

## Confirmatory run

- Toolkit release 0.7, `explore_pretest.R` CODE_VERSION 0.5.0, mode `full`.
- 1,500 training + 500 held-out designs, 1,000 replications per design, seed 732041 (exploratory seed + 730000). The design distribution is the same as in the exploratory run (Table 1 of the manuscript).
- A design is harmful when its Type I error exceeds 0.075 (Bradley's liberal bound). The line at 0.10 is reported as a sensitivity check only and does not enter any verdict.
- All criteria are read from the tables written by `summarize_pretest.R` and are checked mechanically by `check_pretest_hypotheses.R`. Checksums of the frozen scripts:

| File | md5 |
|---|---|
| explore_pretest.R | c4e821ca5e352552cf55b04d61d87d0e |
| summarize_pretest.R | 26848d56c4a3df12fdb0cd8b102c2286 |
| check_pretest_hypotheses.R | 65b2a2b7dd7f8a9bf5f176de0934f006 |

## Definitions

- *Change* is the per-group Shapiro-Wilk rejection rate on the observed data minus the rate on the complete data of the same replications, averaged over designs without outliers in a shape × mechanism cell (`pretest_sw_rates.csv`, column `group_change`).
- *Upper-tail loss* is outcome-driven missingness with larger values more likely to be lost (MNAR_large_lost); *lower-tail loss* is the reverse (MNAR_small_lost).
- *Identical distributions* means equal population means and an SD ratio of 1.
- *Two-stage procedures* are the four combinations of the per-group or residual pre-test with Student's or Welch's t test, with the Mann-Whitney test when the pre-test rejects.
- A hypothesis is **supported** only if every one of its criteria is met. Point estimates are compared with the thresholds as they are; no tolerance is added for Monte Carlo error.

## Hypotheses

**H1 (masking).** Loss of the upper tail hides skewness from the pre-test. For mild skew and for strong skew separately: (a) the change under upper-tail loss is at most −0.10, and (b) it is at least 0.05 more negative than the change under MCAR. *(4 criteria. Exploratory values: −0.32 and −0.20; gaps −0.20 and −0.12.)*

**H2 (random loss and normal populations).** (a) Under MCAR and under group-dependent missingness, the change is at most −0.03 for each non-normal shape (mild skew, strong skew, heavy-tailed). (b) For normal populations, the absolute change is at most 0.04 under each of the four mechanisms. *(10 criteria. Exploratory values: −0.05 to −0.12; normal at most 0.027 in absolute value.)*

**H3 (no benefit of pre-testing).** (a) Welch's t test alone has a lower harmful share than each of the four two-stage procedures, in designs with identical distributions and in all designs with equal means. (b) In designs with identical distributions, Welch's t test alone has a lower mean Type I error than each of the four two-stage procedures. *(12 criteria. Exploratory: 6.2% against 10.2–13.3%; 16.3% against 25.1–28.5%.)*

**H4 (hidden skew leads to the t test when it is least reliable).** In skewed populations (mild or strong) with identical distributions (`pretest_masking.csv`): (a) upper-tail loss raises the share of replications in which the per-group pre-test passes by at least 0.15; (b) this rise exceeds the rise under MCAR and under group-dependent missingness by at least 0.10 each. (c) In designs with equal means (`pretest_conditional.csv`), the Type I error of Welch's t test in replications where the per-group pre-test passed exceeds its error where the pre-test failed by at least 0.02, for mild skew and for strong skew separately. *(5 criteria. Exploratory: rise 0.26; gaps 0.18 and 0.20; conditional gaps 0.054 and 0.047.)*

**H5 (no warning).** The observed per-group rejection rate does not identify outcome-driven missingness. Held-out AUC for outcome-driven against random or group-dependent missingness (`pretest_diagnostic.csv`, `auc_per_group_test`): (a) all designs, between 0.35 and 0.65; (b) normal populations without outliers, below 0.70. *(2 criteria. Exploratory: 0.52 and 0.61.)*

## Reporting

- Every hypothesis is reported with every criterion and its value, whether met or not. A hypothesis that fails is reported as not supported, with the value that failed, and the paper's conclusions are adjusted to it.
- Descriptive results beyond H1–H5 (for example the flip rates, power, and the warpbreaks demonstration) are reported as descriptive, not as confirmatory tests.
- Any change to the run, the scripts, or the criteria after this date is recorded below with its date and reason.

## Deviations

No deviation from the run or the criteria. Two later changes to the reporting scripts, made on 2 October 2026 after the confirmatory results were checked, leave every criterion and every verdict unchanged:

- `check_pretest_hypotheses.R` now counts a missing value as a criterion not met instead of stopping, and runs automatically at the end of the full pipeline (toolkit 0.9; md5 3a107ddd525bfa0b4ac3fd581358e1b8). On the confirmatory results it returns the same 33 values and verdicts.
- `summarize_pretest.R` additionally writes the design counts, the share of failed replications, the median planned sample size and an "unequal spreads" group; the tables used by the criteria are unchanged (toolkit 0.9; md5 f8f2de7bdc5de5bd0a79ad745bbf1f23).
- `make_pretest_figures.R` sets the letter *t* in italics in the Figure 2 axis label and the Figure 3 legend, following APA style (3 October 2026; md5 7535890cc76680bd39b3a69903db8d57). The figures show the same values.
