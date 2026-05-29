# =============================================================================
# 04_classification.R  --  SECTION C.3: Classification Analysis
#
# Uses the TOP FOUR most important variables from the C.2 importance ranking.
#
# This version uses the instructor-provided COURSE FUNCTIONS wherever possible:
#   lin.class()  -- linear classification function coefficients + constants
#                   for each group level (Fisher's LCFs under equal priors)
#   rates()      -- confusion matrix, correct classification rate, error rate
#
# Deliverables (per project spec):
#   1. Linear classification function for each group level.    -> lin.class()
#   2. Apply the functions to Observation #1; predict its group; compare to
#      the actual DIAGN; state whether the prediction was correct.
#   3. Confusion matrix, Apparent Error Rate, Apparent Correct
#      Classification Rate; comment on performance.            -> rates()
# =============================================================================

if (!exists("importance")) source(file.path("R", "03_discriminant.R"))

# --- C.3.0  Select the top-4 variables ----------------------------------------
TOP_K <- 4
top_vars <- head(importance$Variable, TOP_K)
cat("\nTop", TOP_K, "variables used for classification:",
    paste(top_vars, collapse = ", "), "\n")

# Inputs for the course functions: predictor matrix Y and grouping vector group.
Y     <- cancer[, top_vars]
group <- cancer[[GROUP_VAR]]
groups <- levels(group)

# -----------------------------------------------------------------------------
# C.3.1  Linear classification functions (Fisher's) for each group
# -----------------------------------------------------------------------------
# lin.class() returns:
#   $coefs : matrix of coefficients, one ROW per group, one COLUMN per variable
#            (b_g = Sp^{-1} %*% mbar_g)
#   $c.0   : vector of group constants  (c_g = -0.5 * mbar_g' Sp^{-1} mbar_g)
# An observation is assigned to the group with the LARGEST score
#   L_g(x) = c.0[g] + sum_j coefs[g, j] * x_j
lc <- lin.class(Y, group)

# Coefficients come out with groups as rows; transpose to a readable
# "columns = groups" table with a Constant row on top (matches how the
# discriminant functions were displayed in C.2).
coefs_byvar <- t(lc$coefs)              # rows = variables, cols = groups
rownames(coefs_byvar) <- top_vars
colnames(coefs_byvar) <- groups

lcf_table <- rbind(Constant = lc$c.0, coefs_byvar)
colnames(lcf_table) <- groups

cat("\n--- C.3.1 Linear Classification Functions (columns = groups) ---\n")
print(round(lcf_table, 4))
save_table("classification_functions.csv", round(lcf_table, 4))

# Helper to score one observation against every group's LCF (built from the
# lin.class() output so we stay consistent with the course function).
lcf_scores <- function(x_vec) {
  sapply(groups, function(g) lc$c.0[match(g, groups)] +
           sum(coefs_byvar[, g] * x_vec))
}

# -----------------------------------------------------------------------------
# C.3.2  Classify Observation #1 and check it
# -----------------------------------------------------------------------------
X           <- as.matrix(Y)
obs1_x      <- X[1, ]
obs1_scores <- lcf_scores(obs1_x)
obs1_pred   <- groups[which.max(obs1_scores)]
obs1_actual <- as.character(group[1])

cat("\n--- C.3.2 Classifying Observation #1 ---\n")
cat("Predictor values:\n"); print(round(obs1_x, 4))
cat("\nClassification scores L_g(x):\n"); print(round(obs1_scores, 3))
cat("\nPredicted group:", obs1_pred,
    "| Actual group:", obs1_actual,
    "|", if (obs1_pred == obs1_actual) "CORRECT" else "INCORRECT", "\n")

# -----------------------------------------------------------------------------
# C.3.3  Confusion matrix + apparent error / correct rates
# -----------------------------------------------------------------------------
# rates() (method = "l" for LDA) returns the apparent (resubstitution)
# confusion matrix, the correct classification rate, and the error rate.
cls <- rates(Y, group, method = "l")

confusion <- cls[["Confusion Matrix"]]
ACCR      <- cls[["Correct Class Rate"]]   # Apparent Correct Classification Rate
APER      <- cls[["Error Rate"]]           # Apparent Error Rate

cat("\n--- C.3.3 Confusion Matrix (apparent / resubstitution) ---\n")
print(confusion)
save_table("confusion_matrix.csv", as.data.frame.matrix(confusion))

n         <- nrow(Y)
n_correct <- round(ACCR * n)
cat(sprintf("\nApparent Correct Classification Rate: %.4f (%d/%d)\n",
            ACCR, n_correct, n))
cat(sprintf("Apparent Error Rate (APER):           %.4f\n", APER))
cat("Method used:", cls[["Method"]], "\n")

save_table("classification_rates.csv",
           data.frame(Metric = c("ACCR", "APER", "N", "Correct"),
                      Value  = c(round(ACCR, 4), round(APER, 4), n, n_correct)))

cat("\nSection C.3 complete.\n")
