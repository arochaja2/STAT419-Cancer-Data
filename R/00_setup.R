# =============================================================================
# 00_setup.R
# STAT 419 Multivariate Analysis Project -- Cancer (Group 2)
#
# Purpose: Central setup script. Sourced by every other script so that paths,
#          packages, the dataset, the SHARED COURSE FUNCTIONS, and local helpers
#          are defined in ONE place.
#
# Usage:   source("R/00_setup.R")   # run from the project root directory
# =============================================================================

# --- 0. Project root -----------------------------------------------------------
# All scripts assume the working directory is the PROJECT ROOT (the folder that
# contains the data/, R/, and output/ subfolders). If you open the .Rproj file
# (or set the working directory manually) this will be handled for you.
#
#   setwd("C:/path/to/STAT419-Cancer-Data")   # <- edit if running outside RStudio proj

# --- 1. Packages ---------------------------------------------------------------
# MASS    : lda()/qda(); also used internally by the course functions
#           lin.class() and rates().
# Base R covers cor(), pairs(), hist(), manova(), summary.manova().
required_pkgs <- c("MASS")

install_if_missing <- function(pkgs) {
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(required_pkgs)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# --- 2. Paths ------------------------------------------------------------------
PATHS <- list(
  data      = file.path("data", "cancer_419.csv"),
  functions = file.path("R", "_all_customized_functions.R"),
  figures   = file.path("output", "figures"),
  tables    = file.path("output", "tables")
)

dir.create(PATHS$figures, showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$tables,  showWarnings = FALSE, recursive = TRUE)

# --- 3. Shared COURSE functions ------------------------------------------------
# These are the instructor-provided custom functions (formerly
# _all_customized_functions.txt). Per the project spec we use these "relevant
# functions in R" wherever possible. Save that file as an .R file alongside the
# scripts. Key functions used in this project:
#   check.mvnorm.plot()  -- chi-square Q-Q plot to assess multivariate normality
#   discrim()            -- raw + standardized discriminant coefficients
#   discr.sig()          -- significance tests of the discriminant functions
#   partial.F()          -- partial-F significance of each variable (adjusted)
#   discr.plot()         -- LD1 vs LD2 plot (interactive legend)
#   lin.class()          -- Fisher's linear classification functions
#   rates()              -- confusion matrix + correct/error classification rates
#   Box.M()              -- Box's M test for equality of covariance matrices
source(PATHS$functions)

# --- 4. Load data --------------------------------------------------------------
# DIAGN is the grouping variable (1 = low, 2 = moderate, 3 = severe seminal
# vesicle invasion). Everything else is a quantitative predictor.
# fileEncoding = "UTF-8-BOM" strips the byte-order-mark some editors prepend,
# which otherwise turns the first column name into "X...DIAGN".
cancer <- read.csv(PATHS$data, sep = ",", header = TRUE,
                   fileEncoding = "UTF-8-BOM")

# Make the grouping variable a labelled factor for modelling / plotting.
cancer$DIAGN <- factor(cancer$DIAGN,
                       levels = c(1, 2, 3),
                       labels = c("Low", "Moderate", "Severe"))

# Convenience vectors -----------------------------------------------------------
GROUP_VAR   <- "DIAGN"
QUANT_VARS  <- setdiff(names(cancer), GROUP_VAR)   # all non-grouping variables

# --- 5. Local helpers (file I/O only -- not statistical) -----------------------

# Save a base-R plot to output/figures as a PNG. `expr` is the plotting code.
save_plot <- function(filename, expr, width = 700, height = 500, res = 110) {
  path <- file.path(PATHS$figures, filename)
  png(path, width = width, height = height, res = res)
  on.exit(dev.off())
  force(expr)
  message("Saved figure: ", path)
  invisible(path)
}

# Write a data.frame / matrix to output/tables as CSV.
save_table <- function(filename, obj) {
  path <- file.path(PATHS$tables, filename)
  write.csv(obj, path, row.names = TRUE)
  message("Saved table: ", path)
  invisible(path)
}

# --- 6. Quick sanity check (only when run interactively) -----------------------
if (interactive()) {
  cat("Loaded", nrow(cancer), "observations and",
      ncol(cancer), "variables.\n")
  cat("Group sizes:\n"); print(table(cancer$DIAGN))
  cat("Quantitative variables:", paste(QUANT_VARS, collapse = ", "), "\n")
}
