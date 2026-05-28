# =============================================================================
# run_all.R  --  Master script: runs the entire analysis pipeline in order.
#
# From the project ROOT directory (the folder containing data/, R/, output/):
#
#     source("R/run_all.R")
#
# Each step writes its figures to output/figures/ and tables to output/tables/.
# Console output (the bits you paste into the report) is printed as it runs.
#
# The pipeline relies on the instructor-provided course functions in
# R/_all_customized_functions.R (sourced by 00_setup.R), in particular:
#   check.mvnorm.plot, discrim, discr.sig, partial.F, discr.plot,
#   lin.class, rates.
#
# IMPORTANT -- the LD1 vs LD2 plot in C.2.4 uses discr.plot(), which calls
# legend(locator(1)) and waits for a MOUSE CLICK. It is therefore SKIPPED when
# sourcing this master script. After run_all.R finishes, run that one plot
# interactively:
#     source("R/03_discriminant.R")
#     discr.plot(Y, group, leg = levels(group))
#   then click an empty area of the plot to drop the legend.
# =============================================================================

cat("==========================================================\n")
cat(" STAT 419 Project -- Cancer data -- full analysis pipeline \n")
cat("==========================================================\n\n")

steps <- c(
  "R/00_setup.R",          # load data, packages, course functions, helpers
  "R/01_explore.R",        # Section B: histograms, summary stats, MVN Q-Q plot
  "R/02_correlation.R",    # Section C.1: correlation + scatterplot matrix
  "R/03_discriminant.R",   # Section C.2: discrim/discr.sig/partial.F (+ LD plot)
  "R/04_classification.R"  # Section C.3: lin.class, obs #1, rates/confusion
)

for (s in steps) {
  cat("\n\n>>>>>>>>>> Running", s, "<<<<<<<<<<\n")
  source(s, echo = FALSE)
}

cat("\n\n==========================================================\n")
cat(" Done. See output/figures/ and output/tables/ for results.\n")
cat(" Reminder: run discr.plot() interactively for the LD1 vs LD2 plot.\n")
cat("==========================================================\n")
