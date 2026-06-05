# =============================================================================
# 03_discriminant.R  --  SECTION C.2: Discriminant Analysis
#
# Depends on ANALYSIS_VARS / DROP_VARS defined in 02_correlation.R.
#
# Uses the provided COURSE FUNCTIONS wherever possible:
#   discrim()    -- raw ($a) and standardized ($a.stand) discriminant coeffs
#   discr.sig()  -- sequential significance tests of the discriminant functions
#   partial.F()  -- partial-F significance of each variable (adjusted for others)
#   discr.plot() -- plot of the first two linear discriminant functions
#
# Deliverables (per project spec):
#   1. Complete form of all discriminant functions using STANDARDIZED coeffs;
#      rank variable importance from those coefficients.       -> discrim()
#   2. Significance tests for the discriminant functions, with H0, H1,
#      statistic, p-value, conclusion.                         -> discr.sig()
#   3. Significance test for each non-grouping variable, adjusted for the
#      others, with the same reporting.                        -> partial.F()
#   4. Plot of the first two linear discriminant functions.    -> discr.plot()
# =============================================================================

if (!exists("ANALYSIS_VARS")) source(file.path("R", "02_correlation.R"))

set.seed(419)   # reproducibility

# Build the predictor matrix (Y) and grouping vector (group) that the course
# functions expect as their two arguments.
Y     <- cancer[, ANALYSIS_VARS]
group <- cancer[[GROUP_VAR]]

# -----------------------------------------------------------------------------
# C.2.1  Discriminant function coefficients + importance ranking
# -----------------------------------------------------------------------------
# discrim() returns:
#   $a        raw (unstandardized) discriminant coefficients
#   $a.stand  STANDARDIZED coefficients (scaled by within-group SD) -- use these
#             both to WRITE the discriminant functions and to RANK importance.
da <- discrim(Y, group)

n_ld <- ncol(da$a.stand)

std_coef <- da$a.stand
rownames(std_coef) <- ANALYSIS_VARS
colnames(std_coef) <- paste0("LD", seq_len(n_ld))

raw_coef <- da$a
rownames(raw_coef) <- ANALYSIS_VARS
colnames(raw_coef) <- paste0("LD", seq_len(n_ld))

cat("\n--- C.2.1 Raw discriminant coefficients ($a) ---\n")
print(round(raw_coef, 3))

cat("\n--- C.2.1 Standardized discriminant coefficients ($a.stand) ---\n")
print(round(std_coef, 3))
save_table("standardized_coefficients.csv", round(std_coef, 4))

# Variable importance: rank by |standardized coefficient| on LD1 (the function
# that explains the most between-group separation). Discuss LD2 if relevant.
importance <- data.frame(
  Variable  = rownames(std_coef),
  abs_LD1   = abs(std_coef[, "LD1"]),
  row.names = NULL
)
importance <- importance[order(-importance$abs_LD1), ]
importance$Rank <- seq_len(nrow(importance))

cat("\n--- C.2.1 Variable importance ranking (by |LD1| std. coef) ---\n")
print(importance, row.names = FALSE)
save_table("variable_importance.csv", importance)

# -----------------------------------------------------------------------------
# C.2.2  Significance tests for the discriminant functions
# -----------------------------------------------------------------------------
# discr.sig() performs the sequential likelihood-ratio (Wilks' Lambda -> chi-
# square V) tests of the discriminant functions LD1, LD2, ...:
#
#   Test for LD1 (all functions): is there ANY group separation?
#     H0: the discriminant functions explain no separation
#         (equivalently mu_Low = mu_Moderate = mu_Severe)
#     H1: at least the first discriminant function is significant
#   Test for LD2 (after removing LD1): does LD2 add significant separation
#     beyond LD1?  H0: LD2 (and beyond) explain no further separation, etc.
#
# Output columns: Lambda, V (chi-square statistic), p.values.
# Use the p-values to decide which discriminant functions are significant.
discr_sig <- discr.sig(Y, group)

cat("\n--- C.2.2 Significance tests of the discriminant functions ---\n")
print(round(discr_sig, 4))
save_table("discriminant_significance.csv", round(discr_sig, 6))

# -----------------------------------------------------------------------------
# C.2.3  Significance of each variable, adjusted for the others
# -----------------------------------------------------------------------------
# partial.F() computes, for each predictor, the partial-F (and Wilks' Lambda)
# testing its contribution to group separation GIVEN the other predictors:
#
#   For predictor X_j:
#   H0: X_j adds no group-separating information beyond the other predictors
#   H1: X_j does add group-separating information
#
# Output columns: Lambda, F.stat, p.value -- rows are ordered by F (descending),
# i.e. most-to-least important after adjustment. Compare this ordering to the
# standardized-coefficient ranking in C.2.1 in your write-up.
partial_tests <- partial.F(Y, group)

cat("\n--- C.2.3 Per-variable significance (adjusted for other variables) ---\n")
print(round(partial_tests, 4))
save_table("per_variable_significance.csv", round(partial_tests, 6))

# -----------------------------------------------------------------------------
# C.2.4  Plot of the first two discriminant functions (LD1 vs LD2)
# -----------------------------------------------------------------------------
# discr.plot() plots LD1 vs LD2 with a different symbol per group. It calls
# legend(locator(1)), which waits for you to CLICK an empty spot on the graph
# to place the legend. Run this part INTERACTIVELY (not via run_all.R) so you
# can click. After clicking, export the plot from the RStudio Plots pane, or
# wrap it in png()/dev.off() and click within the device.
#
# If locator() locks up the system, edit discr.plot() in
# _all_customized_functions.R and delete the final legend(locator(1), ...) line,
# then annotate the legend manually.

save_plot("discrim.plot.png",
  expr = discr.plot(Y, group, title='Discriminant Plot for Diagnosis Groups', leg = levels(group)))

cat("\nSection C.2 complete. Objects available: da, std_coef, importance,",
    "discr_sig, partial_tests.\n")
