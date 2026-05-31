# =============================================================================
# 02_correlation.R  --  SECTION C.1: Correlated Quantitative Variables
#
# Tasks (per project spec):
#   * correlation matrix for every pair of quantitative variables  [cor()]
#   * scatterplot matrix for all quantitative variables            [pairs()]
#   * flag highly correlated pairs; decide which member to DROP
#
# The decision of WHICH variable to drop is a judgement call you must justify
# in the report. A common rule of thumb: |r| >= 0.80 is "high". This script
# flags candidates; we choose and record the choice in DROP_VARS below.
# =============================================================================

if (!exists("cancer")) source(file.path("R", "00_setup.R"))

# Threshold for "high" correlation -- adjust + justify in the report.
HIGH_CORR_THRESHOLD <- 0.80

# --- C.1.1  Correlation matrix -------------------------------------------------
corr_matrix <- cor(cancer[QUANT_VARS], use = "complete.obs")
corr_matrix_round <- round(corr_matrix, 3)

save_table("correlation_matrix.csv", corr_matrix_round)
cat("\n--- Section C.1: Correlation Matrix ---\n")
print(corr_matrix_round)

# --- C.1.2  Scatterplot matrix -------------------------------------------------
save_plot(
  filename = "scatterplot_matrix.png",
  width = 1000, height = 1000, res = 120,
  expr = pairs(cancer[QUANT_VARS],
               main = "Scatterplot Matrix of Quantitative Variables",
               pch = 19, cex = 0.5,
               col = c("steelblue", "darkorange", "firebrick")[cancer$DIAGN])
)

# --- C.1.3  Flag highly correlated pairs --------------------------------------
# List every off-diagonal pair whose |r| meets/exceeds the threshold.
flag_high_corr <- function(cmat, thresh) {
  pairs_idx <- which(abs(cmat) >= thresh & upper.tri(cmat), arr.ind = TRUE)
  if (nrow(pairs_idx) == 0) {
    return(data.frame(Var1 = character(), Var2 = character(), r = numeric()))
  }
  data.frame(
    Var1 = rownames(cmat)[pairs_idx[, 1]],
    Var2 = colnames(cmat)[pairs_idx[, 2]],
    r    = round(cmat[pairs_idx], 3),
    row.names = NULL
  )
}

high_pairs <- flag_high_corr(corr_matrix, HIGH_CORR_THRESHOLD)
cat("\n--- Pairs with |r| >=", HIGH_CORR_THRESHOLD, "---\n")
if (nrow(high_pairs) == 0) {
  cat("None. No variables flagged for removal at this threshold.\n")
} else {
  print(high_pairs, row.names = FALSE)
  save_table("high_correlation_pairs.csv", high_pairs)
}

# --- C.1.4  Decide which variables to drop ------------------------------------
# >>> EDIT THIS based on high_pairs above and justification. <<<
# If nothing is sufficiently correlated, leave it as character(0).
# Example justification to record in the report: when two predictors are highly
# correlated, keep the one that is (a) cheaper/easier to measure clinically, or
# (b) more interpretable, or (c) more strongly associated with DIAGN.
DROP_VARS <- character(0)        # e.g. c("Capsule")

# The predictor set carried forward into the discriminant analysis (C.2).
ANALYSIS_VARS <- setdiff(QUANT_VARS, DROP_VARS)

cat("\nVariables dropped for C.2:",
    if (length(DROP_VARS)) paste(DROP_VARS, collapse = ", ") else "(none)", "\n")
cat("Variables retained for C.2:", paste(ANALYSIS_VARS, collapse = ", "), "\n")
