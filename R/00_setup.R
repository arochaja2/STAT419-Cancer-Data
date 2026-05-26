# =============================================================================
# 00_setup.R
# STAT 419 Multivariate Analysis Project -- Cancer (Group 2)
#
# Purpose: Central setup script. Sourced by every other script so that paths,
#          packages, the dataset, and shared helpers are defined in ONE place.
#
# Usage:   source("R/00_setup.R")   # run from the project root directory
# =============================================================================

# --- 0. Project root -----------------------------------------------------------
# All scripts assume the working directory is the PROJECT ROOT (the folder that
# contains the data/, R/, and output/ subfolders). If you open the .Rproj file
# (or set the working directory manually) this will be handled for you.
#
#   setwd("C:/path/to/stat419-cancer")   # <- edit if running outside RStudio proj

# --- 1. Packages ---------------------------------------------------------------
# MASS    : lda() for linear discriminant analysis
# car     : (optional) Anova / multivariate test helpers
# Base R covers cor(), pairs(), hist(), manova(), summary.manova().
required_pkgs <- c("MASS")
optional_pkgs <- c("car", "GGally")   # nice-to-have, not strictly required

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
  data    = file.path("data", "cancer_419.csv"),
  figures = file.path("output", "figures"),
  tables  = file.path("output", "tables")
)

dir.create(PATHS$figures, showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$tables,  showWarnings = FALSE, recursive = TRUE)

# --- 3. Load data --------------------------------------------------------------
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

# --- 4. Shared helpers ---------------------------------------------------------

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

# --- 5. Quick sanity check (only when run interactively) -----------------------
if (interactive()) {
  cat("Loaded", nrow(cancer), "observations and",
      ncol(cancer), "variables.\n")
  cat("Group sizes:\n"); print(table(cancer$DIAGN))
  cat("Quantitative variables:", paste(QUANT_VARS, collapse = ", "), "\n")
}
