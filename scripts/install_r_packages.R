#!/usr/bin/env Rscript
# Install the public example's R dependencies into the current R library.
# R >= 4.4 is recommended. Versions are captured by each workflow's sessionInfo.
# This is an installer, not an environment lockfile.

args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--with-doublets", "--help"))) stop("Unknown installer option.")
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/install_r_packages.R [--with-doublets]\n")
  quit(status = 0)
}
options(repos = c(CRAN = "https://cloud.r-project.org"))
packages <- c("Matrix", "Seurat", "SoupX", "ggplot2")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing)
}
still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) {
  stop("Installation incomplete: ", paste(still_missing, collapse = ", "))
}
if ("--with-doublets" %in% args) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  optional <- c("SingleCellExperiment", "scDblFinder", "BiocParallel")
  missing_optional <- optional[!vapply(optional, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_optional)) BiocManager::install(missing_optional, ask = FALSE, update = FALSE)
  if (!all(vapply(optional, requireNamespace, logical(1), quietly = TRUE))) {
    stop("Optional Bioconductor package installation incomplete.")
  }
  packages <- c(packages, optional)
}
print(data.frame(
  package = packages,
  version = vapply(packages, function(x) as.character(packageVersion(x)), character(1))
))
