# =============================================================================
# 04_classification.R  --  SECTION C.3: Classification Analysis
#
# Uses the TOP FOUR most important variables from the C.2 importance ranking.
#
# Deliverables (per project spec):
#   1. Linear classification function for each group level (Fisher's LCFs).
#   2. Apply the functions to Observation #1; predict its group; compare to
#      the actual DIAGN; state whether the prediction was correct.
#   3. Confusion matrix, Apparent Error Rate (APER), Apparent Correct
#      Classification Rate; comment on performance.
# =============================================================================

if (!exists("importance")) source(file.path("R", "03_discriminant.R"))

# --- C.3.0  Select the top-4 variables ----------------------------------------
TOP_K <- 4
top_vars <- head(importance$Variable, TOP_K)
cat("\nTop", TOP_K, "variables used for classification:",
    paste(top_vars, collapse = ", "), "\n")

class_df <- cancer[, c(GROUP_VAR, top_vars)]

# -----------------------------------------------------------------------------
# C.3.1  Linear classification functions (Fisher's) for each group
# -----------------------------------------------------------------------------
# Fisher's linear classification function for group g:
#     L_g(x) = c_g + sum_j b_gj * x_j
# where, with pooled within-group covariance Sp and group mean vector mbar_g:
#     b_g  =  Sp^{-1} %*% mbar_g
#     c_g  = -0.5 * mbar_g' Sp^{-1} mbar_g  (+ log prior, here equal priors)
# An observation is assigned to the group with the LARGEST L_g(x).

groups   <- levels(class_df$DIAGN)
X        <- as.matrix(class_df[top_vars])
grp      <- class_df$DIAGN
k        <- length(groups)
n        <- nrow(X)
p        <- length(top_vars)

# Group means (rows = groups).
group_means <- t(sapply(groups, function(g) colMeans(X[grp == g, , drop = FALSE])))

# Pooled within-group covariance Sp.
Sp <- Reduce(`+`, lapply(groups, function(g) {
  Xi <- X[grp == g, , drop = FALSE]
  (nrow(Xi) - 1) * cov(Xi)
})) / (n - k)

Sp_inv <- solve(Sp)

# Coefficients b_g (one column per group) and constants c_g.
lcf_b <- sapply(groups, function(g) Sp_inv %*% group_means[g, ])
rownames(lcf_b) <- top_vars
lcf_c <- sapply(groups, function(g) {
  m <- group_means[g, ]
  -0.5 * as.numeric(t(m) %*% Sp_inv %*% m)
})

# Assemble a readable coefficient table: constant + one row per variable.
lcf_table <- rbind(Constant = lcf_c, lcf_b)
colnames(lcf_table) <- groups

cat("\n--- C.3.1 Linear Classification Functions (columns = groups) ---\n")
print(round(lcf_table, 4))
save_table("classification_functions.csv", round(lcf_table, 4))

# Helper: score an observation against every group's LCF; return the scores.
lcf_scores <- function(x_vec) {
  sapply(groups, function(g) lcf_c[g] + sum(lcf_b[, g] * x_vec))
}

# -----------------------------------------------------------------------------
# C.3.2  Classify Observation #1 and check it
# -----------------------------------------------------------------------------
obs1_x      <- X[1, ]
obs1_scores <- lcf_scores(obs1_x)
obs1_pred   <- groups[which.max(obs1_scores)]
obs1_actual <- as.character(grp[1])

cat("\n--- C.3.2 Classifying Observation #1 ---\n")
cat("Predictor values:\n"); print(round(obs1_x, 4))
cat("\nClassification scores L_g(x):\n"); print(round(obs1_scores, 3))
cat("\nPredicted group:", obs1_pred,
    "| Actual group:", obs1_actual,
    "|", if (obs1_pred == obs1_actual) "CORRECT" else "INCORRECT", "\n")

# -----------------------------------------------------------------------------
# C.3.3  Confusion matrix + apparent error / correct rates
# -----------------------------------------------------------------------------
# Classify every observation with the LCFs (apparent = resubstitution).
all_scores <- t(apply(X, 1, lcf_scores))
pred_all   <- factor(groups[max.col(all_scores)], levels = groups)

confusion <- table(Actual = grp, Predicted = pred_all)
cat("\n--- C.3.3 Confusion Matrix (apparent / resubstitution) ---\n")
print(confusion)
save_table("confusion_matrix.csv", as.data.frame.matrix(confusion))

n_correct <- sum(diag(confusion))
APER      <- 1 - n_correct / n          # Apparent Error Rate
ACCR      <- n_correct / n              # Apparent Correct Classification Rate

cat(sprintf("\nApparent Correct Classification Rate: %.4f (%d/%d)\n",
            ACCR, n_correct, n))
cat(sprintf("Apparent Error Rate (APER):           %.4f\n", APER))

# Cross-check against MASS::lda with the SAME 4 variables (should match closely;
# lda uses the same Fisher rule under equal priors). Useful as a sanity test.
lda_top  <- lda(DIAGN ~ ., data = class_df)
lda_pred <- predict(lda_top)$class
cat("\n(Sanity check) lda() apparent accuracy:",
    round(mean(lda_pred == grp), 4), "\n")

save_table("classification_rates.csv",
           data.frame(Metric = c("ACCR", "APER", "N", "Correct"),
                      Value  = c(round(ACCR, 4), round(APER, 4), n, n_correct)))

cat("\nSection C.3 complete.\n")
