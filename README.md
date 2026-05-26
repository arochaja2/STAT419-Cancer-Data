# STAT 419 Project — Cancer Data (Group 2)

Discriminant and classification analysis of prostate cancer prognostic
measurements, for the STAT 419 Multivariate Analysis course project.

The data set records 8 measurements on 97 men with advanced prostate cancer.
The grouping variable is **DIAGN** — severity of seminal vesicle invasion
(1 = low, 2 = moderate, 3 = severe). The remaining seven variables are
quantitative predictors (PSA, Volume, Weight, Age, BPH, Capsule, Gleason).

> This repo contains **only the coding portion and project structure**. The
> written report (cover page, prose for Sections A–E) is produced separately.

## Project structure

```
stat419-cancer/
├── README.md
├── stat419-cancer.Rproj        # open this in RStudio so paths resolve to root
├── .gitignore
├── data/
│   ├── cancer_419.csv          # the data set
│   └── data_description.pdf     # variable descriptions (provided)
├── R/
│   ├── 00_setup.R              # load data + packages, define paths & helpers
│   ├── 01_explore.R            # Section B: histograms + summary statistics
│   ├── 02_correlation.R        # Section C.1: correlation / multicollinearity
│   ├── 03_discriminant.R       # Section C.2: discriminant analysis
│   ├── 04_classification.R     # Section C.3: classification analysis
│   └── run_all.R               # runs the whole pipeline in order
├── output/
│   ├── figures/                # generated PNGs (git-ignored)
│   └── tables/                 # generated CSVs (git-ignored)
└── report/                     # place the written report here
```

## How to run

You need R (the only required package, **MASS**, ships with base R).

**Option A — RStudio:** open `stat419-cancer.Rproj`, then in the console:

```r
source("R/run_all.R")
```

**Option B — command line**, from the project root:

```sh
Rscript R/run_all.R
```

Each script can also be run on its own; if its dependencies aren't loaded yet
it will source the earlier scripts automatically. Figures land in
`output/figures/`, tables in `output/tables/`, and the console prints the
results you paste into the report.

## How scripts map to report sections

| Report section | Script | Produces |
|---|---|---|
| B — Graphs & summary stats | `01_explore.R` | histogram per variable, summary stats table |
| C.1 — Multicollinearity | `02_correlation.R` | correlation matrix, scatterplot matrix, flagged pairs |
| C.2 — Discriminant analysis | `03_discriminant.R` | standardized coefficients, importance ranking, Wilks' Lambda + per-variable tests, LD1×LD2 plot |
| C.3 — Classification | `04_classification.R` | linear classification functions, Obs #1 prediction, confusion matrix, APER & correct rate |

## Decisions you must make and justify in the report

These are intentionally left as explicit, editable choices rather than
hard-coded — the project rewards justified judgement, not a single "right"
answer:

- **`DROP_VARS` in `02_correlation.R`** — which (if any) highly correlated
  variable to remove. The script flags pairs with |r| ≥ 0.80 (the threshold is
  also adjustable); you choose which member to keep and write up why. It is
  fine for this to be empty if nothing is sufficiently correlated.
- **`HIGH_CORR_THRESHOLD`** — the cutoff for "high" correlation.
- **`TOP_K` in `04_classification.R`** — fixed at 4 per the project spec
  ("top four most important variables").

## Notes

- The grouping variable is read as a labelled factor (Low / Moderate / Severe).
- Standardized discriminant coefficients are computed using the pooled
  within-group SD (the standard approach for comparing variable importance).
- The LD1 × LD2 plot uses a fixed-corner legend rather than the interactive
  `locator()` placement mentioned in the spec, to avoid the lock-up the
  instructions warn about. Adjust the legend position if it overlaps points.
- The hand-built linear classification functions in `04_classification.R` are
  cross-checked against `MASS::lda()`; the apparent accuracies should match.
```
