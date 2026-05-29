# =============================================================================
# 01_explore.R  --  SECTION B: Graphs and Summary Statistics
#
# Produces, for the non-grouping variables:
#   * a histogram for every quantitative variable      -> output/figures/
#   * a summary statistics table (mean, median, SD)     -> output/tables/
#   * a chi-square Q-Q plot to assess multivariate normality, using the
#     course function check.mvnorm.plot()               -> output/figures/
#
# Run with:  source("R/01_explore.R")   (after sourcing 00_setup.R, or it will
#            source it for you).
#
# NOTE: Section B has no dedicated course function for histograms / summary
#       statistics (those are plain base R: hist(), mean(), median(), sd()),
#       so those parts are unchanged. The one relevant custom function here is
#       check.mvnorm.plot(), used as a normality diagnostic.
# =============================================================================

if (!exists("cancer")) source(file.path("R", "00_setup.R"))

# --- B.1  Histograms -----------------------------------------------------------
# One PNG per quantitative variable. Comment on shape/skew in the report.
for (v in QUANT_VARS) {
  save_plot(
    filename = paste0("hist_", v, ".png"),
    expr = hist(cancer[[v]],
                main = paste("Histogram of", v),
                xlab = v,
                col  = "grey80",
                border = "white")
  )
}

# Optional: all histograms on a single panel for a quick overview in the report.
save_plot(
  filename = "hist_panel_all.png",
  width = 900, height = 700,
  expr = {
    op <- par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
    for (v in QUANT_VARS) {
      hist(cancer[[v]], main = v, xlab = "", col = "grey80", border = "white")
    }
    par(op)
  }
)

# --- B.2  Summary statistics ---------------------------------------------------
# Required: mean, median, SD for every quantitative variable.
summary_stats <- data.frame(
  Variable = QUANT_VARS,
  Mean     = sapply(cancer[QUANT_VARS], mean,   na.rm = TRUE),
  Median   = sapply(cancer[QUANT_VARS], median, na.rm = TRUE),
  SD       = sapply(cancer[QUANT_VARS], sd,     na.rm = TRUE),
  row.names = NULL
)

# Round for presentation, then save and print.
summary_stats[ , -1] <- round(summary_stats[ , -1], 3)
save_table("summary_statistics.csv", summary_stats)

cat("\n--- Section B: Summary Statistics ---\n")
print(summary_stats, row.names = FALSE)

# Base-R full summary (min / quartiles / max) for additional commentary.
cat("\n--- Full summary() for reference ---\n")
print(summary(cancer[QUANT_VARS]))

# --- B.3  Multivariate normality diagnostic (course function) ------------------
# check.mvnorm.plot() draws a chi-square probability (Q-Q) plot of the squared
# Mahalanobis distances. Points near the 45-degree line support multivariate
# normality, an assumption behind LDA / Hotelling-type procedures. Comment on
# departures from the line in the report.
save_plot(
  filename = "mvnorm_chisq_qq.png",
  width = 700, height = 600, res = 120,
  expr = check.mvnorm.plot(cancer[QUANT_VARS])
)
