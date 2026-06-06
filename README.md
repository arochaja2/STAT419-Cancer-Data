# STAT 419 — Prostate Cancer Discriminant Analysis (Group 2)

Discriminant and classification analysis of prostate cancer prognostic
measurements, for the STAT 419 Multivariate Analysis course project.

**[View the rendered PDF report](https://github.com/arochaja2/STAT419-Cancer-Data/blob/main/report.pdf)**

## Overview

The dataset records 8 measurements on 97 men with advanced prostate cancer.
The grouping variable is **DIAGN** — severity of seminal vesicle invasion
(Low / Moderate / Severe, nearly balanced at ~32 per group). The seven
quantitative predictors are PSA, Volume, Weight, Age, BPH, Capsule, and Gleason.

The analysis proceeded in four stages:

1. **Exploratory analysis** — histograms and summary statistics for each variable.
2. **Transformation & outlier removal** — log transforms applied to right-skewed
   variables (PSA, Volume, Weight); BPH and Capsule dropped for persistent
   non-normality; one high-leverage outlier (obs #41, extreme Weight) removed.
3. **Discriminant analysis** — two linear discriminant functions fit on the
   remaining five variables. LD1 is significant (Wilks' Lambda test, *p* < 0.025);
   LD2 is not. PSA is by far the most important variable by standardized
   coefficient and partial F-test.
4. **Classification** — linear classification functions built on the top 4
   variables by LD1 importance. Resubstitution ACCR ≈ 79%, with all
   misclassifications between adjacent severity groups (no Low ↔ Severe errors).

## Project structure

```
STAT419-Cancer-Data/
├── report.qmd                    # single source file — all analysis and prose
├── report.pdf                    # rendered output (auto-built by CI)
├── _all_customized_functions.R   # course-provided multivariate helper functions
├── data/
│   ├── cancer_419.csv            # dataset (97 obs × 8 variables)
│   └── data_description.pdf      # variable descriptions
└── .github/workflows/
    └── render-report.yml         # CI: renders report.pdf on push to main
```

## How to render

Requires R and [Quarto](https://quarto.org). The only R package needed is
**MASS** (ships with base R).

**RStudio:** open `stat419-cancer.Rproj`, then click *Render* on `report.qmd`, or:

```r
quarto::quarto_render("report.qmd")
```

**Command line:**

```sh
quarto render report.qmd --to pdf
```

The PDF is also built automatically by GitHub Actions on every push to `main`
that touches `report.qmd`.

## Key findings

- **PSA** is the dominant predictor — it is the only variable with a
  significant partial F-test after adjusting for the others.
- **LD1** alone captures most of the group separation; LD2 is not significant.
- No multicollinearity was flagged among the five retained variables (all
  pairwise |r| < 0.80).
