# =============================================================================
# 03_discriminant.R  --  SECTION C.2: Discriminant Analysis
#
# Depends on ANALYSIS_VARS / DROP_VARS defined in 02_correlation.R.
#
# Deliverables (per project spec):
#   1. Complete form of all discriminant functions using STANDARDIZED coeffs;
#      rank variable importance from those coefficients.
#   2. Significance tests for the discriminant functions
#      (Wilks' Lambda / MANOVA), with H0, H1, statistic, p-value, conclusion.
#   3. Significance test for each non-grouping variable, adjusted for the
#      others (Type-style MANOVA test per variable), with the same reporting.
#   4. Plot of the first two linear discriminant functions (LD1 vs LD2).
# =============================================================================

if (!exists("ANALYSIS_VARS")) source(file.path("R", "02_correlation.R"))

set.seed(419)   # reproducibility (only matters if you later add CV)

# Build a clean modelling frame: grouping factor + retained predictors.
model_df <- cancer[, c(GROUP_VAR, ANALYSIS_VARS)]

# -----------------------------------------------------------------------------
# C.2.0  Fit the LDA model
# -----------------------------------------------------------------------------
lda_fit <- lda(DIAGN ~ ., data = model_df)
cat("\n=== LDA fit ===\n")
print(lda_fit)

# Number of discriminant functions = min(k-1, p)
n_ld <- ncol(lda_fit$scaling)

# -----------------------------------------------------------------------------
# C.2.1  Standardized discriminant function coefficients + importance ranking
# -----------------------------------------------------------------------------
# lda() returns RAW (unstandardized) coefficients in $scaling. To compare
# variable importance fairly we standardize them by the WITHIN-GROUP standard
# deviation of each predictor (the standard approach in discriminant analysis).

within_sd <- function(x, grp) {
  # pooled within-group SD for one variable
  k  <- nlevels(grp)
  n  <- length(x)
  ss <- sum(tapply(seq_along(x), grp, function(i) sum((x[i] - mean(x[i]))^2)))
  sqrt(ss / (n - k))
}

w_sd <- sapply(ANALYSIS_VARS, function(v) within_sd(model_df[[v]], model_df$DIAGN))

# Standardized coefficients = raw coefficient * within-group SD of predictor
std_coef <- lda_fit$scaling * w_sd
colnames(std_coef) <- paste0("LD", seq_len(n_ld))

cat("\n--- C.2.1 Standardized discriminant coefficients ---\n")
print(round(std_coef, 3))
save_table("standardized_coefficients.csv", round(std_coef, 4))

# Variable importance: rank by |standardized coefficient| on LD1 (the function
# that explains the most between-group separation). You may also discuss LD2.
importance <- data.frame(
  Variable   = rownames(std_coef),
  abs_LD1    = abs(std_coef[, "LD1"]),
  row.names  = NULL
)
importance <- importance[order(-importance$abs_LD1), ]
importance$Rank <- seq_len(nrow(importance))

cat("\n--- C.2.1 Variable importance ranking (by |LD1| std. coef) ---\n")
print(importance, row.names = FALSE)
save_table("variable_importance.csv", importance)

# Proportion of between-group variance explained by each discriminant function.
prop_trace <- lda_fit$svd^2 / sum(lda_fit$svd^2)
cat("\nProportion of trace (separation) by function:\n")
print(round(setNames(prop_trace, paste0("LD", seq_along(prop_trace))), 4))

# -----------------------------------------------------------------------------
# C.2.2  Significance tests for the discriminant functions (Wilks' Lambda)
# -----------------------------------------------------------------------------
# Test whether the group mean vectors differ at all -> equivalent to testing
# the (joint) discriminant functions via a one-way MANOVA.
#
#   H0: mu_Low = mu_Moderate = mu_Severe   (group mean vectors are equal;
#                                            discriminant functions explain
#                                            no real separation)
#   H1: at least one group mean vector differs
#
# Reported with Wilks' Lambda (the classic discriminant-analysis statistic).
manova_fit <- manova(as.matrix(model_df[ANALYSIS_VARS]) ~ model_df$DIAGN)

cat("\n--- C.2.2 Test of discriminant functions: Wilks' Lambda ---\n")
wilks <- summary(manova_fit, test = "Wilks")
print(wilks)

# Also report Pillai (more robust) for completeness/discussion.
cat("\n(Reference) Pillai's Trace:\n")
print(summary(manova_fit, test = "Pillai"))

save_table("manova_wilks.csv", as.data.frame(wilks$stats))

# NOTE on "per-function" sequential tests:
# A full dimensionality test (is LD2 significant given LD1?) can be done with a
# sequential Bartlett chi-square on the eigenvalues. If your course covered it,
# the eigenvalues are (lda_fit$svd)^2; see report notes. Otherwise the MANOVA
# above is the expected deliverable.

# -----------------------------------------------------------------------------
# C.2.3  Significance of each variable, adjusted for the others
# -----------------------------------------------------------------------------
# For each predictor, test its contribution GIVEN the others are in the model.
# Implemented as: full MANOVA vs MANOVA without that predictor isn't the right
# framing for a grouping test; instead we test each predictor's group effect
# while adjusting for the rest using a univariate ANCOVA-style F (the standard
# "additional information" / partial-F test taught in discriminant analysis).
#
#   For predictor X_j:
#   H0: X_j adds no group-separating information beyond the other predictors
#   H1: X_j does add group-separating information
#
# We obtain this from the multivariate model's per-variable test: fit
# lm(X_j ~ other predictors + DIAGN) and test the DIAGN term (partial F).
partial_tests <- lapply(ANALYSIS_VARS, function(v) {
  others <- setdiff(ANALYSIS_VARS, v)
  rhs    <- paste(c(others, GROUP_VAR), collapse = " + ")
  fml    <- as.formula(paste(v, "~", rhs))
  fit    <- lm(fml, data = model_df)
  av     <- anova(fit)
  grp_row <- av[GROUP_VAR, , drop = FALSE]
  data.frame(
    Variable = v,
    F_value  = round(grp_row[["F value"]], 3),
    df1      = grp_row[["Df"]],
    df2      = av["Residuals", "Df"],
    p_value  = signif(grp_row[["Pr(>F)"]], 4),
    row.names = NULL
  )
})
partial_tests <- do.call(rbind, partial_tests)
partial_tests$Significant_05 <- ifelse(partial_tests$p_value < 0.05, "Yes", "No")

cat("\n--- C.2.3 Per-variable significance (adjusted for other variables) ---\n")
print(partial_tests, row.names = FALSE)
save_table("per_variable_significance.csv", partial_tests)

# -----------------------------------------------------------------------------
# C.2.4  Plot of the first two discriminant functions (LD1 vs LD2)
# -----------------------------------------------------------------------------
lda_scores <- as.data.frame(predict(lda_fit)$x)   # columns LD1, LD2, ...
lda_scores$DIAGN <- model_df$DIAGN

grp_cols <- c(Low = "steelblue", Moderate = "darkorange", Severe = "firebrick")

if (n_ld >= 2) {
  save_plot(
    filename = "lda_LD1_LD2.png",
    width = 800, height = 650, res = 120,
    expr = {
      plot(lda_scores$LD1, lda_scores$LD2,
           col = grp_cols[lda_scores$DIAGN], pch = 19,
           xlab = "LD1", ylab = "LD2",
           main = "First Two Linear Discriminant Functions")
      # NOTE: project text mentions an interactive legend(locator(1)). To avoid
      # locking up, we place the legend at a fixed corner instead.
      legend("topright", legend = levels(lda_scores$DIAGN),
             col = grp_cols, pch = 19, bty = "n")
    }
  )
} else {
  # Only one discriminant function exists -> plot LD1 by group instead.
  save_plot(
    filename = "lda_LD1_by_group.png",
    expr = stripchart(LD1 ~ DIAGN, data = lda_scores,
                      vertical = TRUE, method = "jitter", pch = 19,
                      col = grp_cols,
                      main = "LD1 Scores by Group", ylab = "LD1")
  )
}

cat("\nSection C.2 complete. Objects available: lda_fit, std_coef, importance,",
    "manova_fit, partial_tests, lda_scores.\n")
