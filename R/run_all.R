# =============================================================================
# run_all.R  --  Master script: runs the entire analysis pipeline in order.
#
# From the project ROOT directory (the folder containing data/, R/, output/):
#
#     source("R/run_all.R")
#
# Each step writes its figures to output/figures/ and tables to output/tables/.
# Console output (the bits you paste into the report) is printed as it runs.
# =============================================================================

cat("==========================================================\n")
cat(" STAT 419 Project -- Cancer data -- full analysis pipeline \n")
cat("==========================================================\n\n")

steps <- c(
  "R/00_setup.R",          # load data, packages, helpers
  "R/01_explore.R",        # Section B: histograms + summary stats
  "R/02_correlation.R",    # Section C.1: correlation + scatterplot matrix
  "R/03_discriminant.R",   # Section C.2: LDA, coeffs, tests, LD plot
  "R/04_classification.R"  # Section C.3: LCFs, obs #1, confusion matrix
)

for (s in steps) {
  cat("\n\n>>>>>>>>>> Running", s, "<<<<<<<<<<\n")
  source(s, echo = FALSE)
}

cat("\n\n==========================================================\n")
cat(" Done. See output/figures/ and output/tables/ for results.\n")
cat("==========================================================\n")
