# Shapiro-Wilk pre-testing under missing data

Code for the paper *How Do Missing Data Change the Normality Verdict? Shapiro-Wilk Pre-Testing and Test Choice in Two-Group Comparisons* (Eungi Kim, Keimyung University).

| File | Purpose |
|---|---|
| `explore_pretest.R` | The simulation (designs, data, missingness, pre-tests, tests) |
| `summarize_pretest.R` | Summary tables reported in the paper |
| `check_pretest_hypotheses.R` | Checks the prespecified hypotheses (Table 4) |
| `make_pretest_figures.R` | Figures 1–3 (PNG) |
| `demo_warpbreaks.R` | Real-data demonstration (Table 3) |
| `PROTOCOL_pretest_v1.md` | Prespecified hypotheses, frozen before the confirmatory run, with the change log |

Requires R 4.0 or newer (base R only).

## How to run

Put the files in one folder, then in R:

```r
setwd("<folder with these files>")
Sys.setenv(SIM_OUTPUT_ROOT = "C:/pretest", SIM_MODE = "full")   # "medium" = exploratory run
source("explore_pretest.R")            # about 5 hours for "full"; resumes if interrupted
source("summarize_pretest.R")
source("check_pretest_hypotheses.R")
source("make_pretest_figures.R")
source("demo_warpbreaks.R")
```

Output folders are created automatically.

## Seeds

Exploratory run (`medium`): 2041. Confirmatory run (`full`): 732041.

## License

MIT (see `LICENSE`)
