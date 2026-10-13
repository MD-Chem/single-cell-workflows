#!/usr/bin/env Rscript
# Generic ambient-RNA correction example for the repository's synthetic data.
# References: https://github.com/constantAmateur/SoupX
#             https://satijalab.org/seurat/articles/pbmc3k_tutorial
# No marker panel or cell-type labels are supplied by this script.

parse_args <- function(args) {
  defaults <- list("input-dir" = "generated", output = "outputs/soupx", rho = NULL)
  if ("--help" %in% args) {
    cat(paste0(
      "Usage: Rscript scripts/soupx_workflow.R [--input-dir generated] ",
      "[--output outputs/soupx] [--rho 0.10]\n",
      "Default: estimate contamination using SoupX::autoEstCont.\n",
      "--rho: explicitly assume a fixed contamination fraction for a demo.\n",
      "Input must contain raw/ and filtered/ with Matrix Market count files.\n"
    ))
    quit(status = 0)
  }
  if (length(args) %% 2L != 0L) stop("Each option requires a value; see --help.")
  for (pair in seq_len(length(args) %/% 2L)) {
    i <- 2L * pair - 1L
    key <- sub("^--", "", args[[i]])
    if (!startsWith(args[[i]], "--") || !key %in% names(defaults)) {
      stop("Unknown option: ", args[[i]])
    }
    defaults[[key]] <- args[[i + 1L]]
  }
  defaults
}

read_counts <- function(directory) {
  required <- file.path(directory, c("matrix.mtx", "features.tsv", "barcodes.tsv"))
  if (!all(file.exists(required))) stop("Missing matrix.mtx, features.tsv or barcodes.tsv in input.")
  counts <- as(Matrix::readMM(required[[1L]]), "dgCMatrix")
  features <- read.delim(required[[2L]], header = FALSE, stringsAsFactors = FALSE,
                         check.names = FALSE, quote = "", comment.char = "")
  barcodes <- read.delim(required[[3L]], header = FALSE, stringsAsFactors = FALSE,
                         check.names = FALSE, quote = "", comment.char = "")[[1L]]
  ids <- as.character(features[[1L]])
  barcodes <- as.character(barcodes)
  if (nrow(counts) != length(ids) || ncol(counts) != length(barcodes)) {
    stop("Matrix dimensions do not match feature and barcode tables.")
  }
  if (anyNA(ids) || anyNA(barcodes) || any(!nzchar(ids)) || any(!nzchar(barcodes)) ||
      anyDuplicated(ids) || anyDuplicated(barcodes)) {
    stop("Feature IDs and barcodes must be nonempty and unique.")
  }
  if (any(!is.finite(counts@x)) || any(counts@x < 0) ||
      any(abs(counts@x - round(counts@x)) > 1e-8)) {
    stop("Input counts must be finite, nonnegative integers.")
  }
  dimnames(counts) <- list(ids, barcodes)
  if (!all(grepl("^Gene[0-9]{4,}$", ids)) ||
      !all(grepl("^(SIMCELL|EMPTY)[0-9]{4,}$", barcodes))) {
    stop("This public example accepts artificial Gene/SIMCELL/EMPTY IDs only.")
  }
  if (ncol(features) != 3L || !identical(as.character(features[[2L]]), ids) ||
      !all(features[[3L]] == "Gene Expression")) {
    stop("Expected three-column synthetic feature table with identical artificial IDs/names.")
  }
  list(counts = counts, features = features)
}

read_demo_metadata <- function(directory, cell_ids) {
  provenance_path <- file.path(directory, "provenance.csv")
  metadata_path <- file.path(directory, "synthetic_metadata.csv")
  if (!file.exists(provenance_path) || !file.exists(metadata_path)) {
    stop("Run the repository's synthetic generator first; provenance and metadata are required.")
  }
  provenance <- read.csv(provenance_path, stringsAsFactors = FALSE)
  if (!all(c("key", "value") %in% names(provenance))) stop("Invalid provenance table.")
  origin <- provenance$value[provenance$key == "data_origin"]
  if (length(origin) != 1L || is.na(origin) || origin != "independent_simulation_v1") {
    stop("Expected the repository's independent synthetic-data provenance marker.")
  }
  metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)
  fields <- c("barcode", "simulation_population", "simulation_batch")
  if (!all(fields %in% names(metadata))) stop("Synthetic metadata is missing required columns.")
  metadata <- metadata[, fields, drop = FALSE]
  if (anyNA(metadata) || anyDuplicated(metadata$barcode) ||
      !setequal(metadata$barcode, cell_ids) ||
      !all(grepl("^population_[0-9]+$", metadata$simulation_population)) ||
      !all(grepl("^batch_[0-9]+$", metadata$simulation_batch))) {
    stop("Synthetic metadata does not match the count matrix or expected artificial labels.")
  }
  list(origin = origin, cells = metadata[match(cell_ids, metadata$barcode), , drop = FALSE])
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  for (package in c("Matrix", "Seurat", "SoupX")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Missing R package ", package, "; run scripts/install_r_packages.R first.")
    }
  }
  seed <- 42L
  set.seed(seed)
  rho <- NULL
  if (!is.null(args$rho)) {
    rho <- suppressWarnings(as.numeric(args$rho))
    if (length(rho) != 1L || !is.finite(rho) || rho < 0 || rho >= 1) {
      stop("--rho must be a number in [0, 1).")
    }
  }
  raw <- read_counts(file.path(args[["input-dir"]], "raw"))
  filtered <- read_counts(file.path(args[["input-dir"]], "filtered"))
  tod <- raw$counts
  toc <- filtered$counts
  demo <- read_demo_metadata(args[["input-dir"]], colnames(toc))
  if (!identical(rownames(tod), rownames(toc))) {
    stop("Raw and filtered feature IDs must have the same order.")
  }
  if (!all(colnames(toc) %in% colnames(tod))) {
    stop("Every filtered barcode must occur in the raw matrix.")
  }
  if (Matrix::nnzero(tod[, colnames(toc), drop = FALSE] - toc) != 0) {
    stop("Filtered counts must match the corresponding raw columns before correction.")
  }
  if (ncol(tod) <= ncol(toc)) stop("The raw matrix must include additional background droplets.")
  if (ncol(toc) < 20L || sum(Matrix::rowSums(toc > 0) >= 3) < 20L) {
    stop("This example needs at least 20 cells and 20 genes detected in at least 3 cells.")
  }
  if (any(Matrix::colSums(toc) == 0)) stop("Filtered cells must have positive count totals.")
  if (dir.exists(args$output) && length(list.files(args$output, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory is nonempty; choose a new --output directory.")
  }

  # Provisional, unsupervised clusters provide SoupX with expression structure.
  provisional <- Seurat::CreateSeuratObject(toc, project = "synthetic_demo")
  provisional <- Seurat::NormalizeData(provisional, verbose = FALSE)
  provisional <- Seurat::FindVariableFeatures(
    provisional, nfeatures = min(2000L, nrow(provisional)), verbose = FALSE
  )
  provisional <- Seurat::ScaleData(provisional, verbose = FALSE)
  npcs <- min(30L, ncol(provisional) - 1L,
              length(Seurat::VariableFeatures(provisional)) - 1L)
  if (npcs < 2L) stop("Too few variable features for provisional clustering.")
  provisional <- Seurat::RunPCA(provisional, npcs = npcs, seed.use = seed, verbose = FALSE)
  dims <- seq_len(ncol(Seurat::Embeddings(provisional, "pca")))
  if (length(dims) < 2L) stop("PCA produced fewer than two components.")
  provisional <- Seurat::FindNeighbors(
    provisional, dims = dims, k.param = min(20L, ncol(provisional) - 1L), verbose = FALSE
  )
  provisional <- Seurat::FindClusters(
    provisional, resolution = 0.5, random.seed = seed, verbose = FALSE
  )
  clusters <- setNames(as.character(Seurat::Idents(provisional)), colnames(provisional))
  channel <- SoupX::SoupChannel(tod = tod, toc = toc)
  channel <- SoupX::setClusters(channel, clusters)
  if (is.null(rho)) {
    # Auto-estimation can fail on homogeneous or small synthetic inputs.
    # Fail explicitly: never substitute a made-up estimate.
    channel <- tryCatch(
      SoupX::autoEstCont(channel, doPlot = FALSE),
      error = function(e) stop(
        "SoupX auto-estimation failed: ", conditionMessage(e),
        "\nFor the synthetic demo only, use --rho 0.10 to explicitly assume 10% contamination.",
        call. = FALSE
      )
    )
    mode <- "estimated_by_autoEstCont"
  } else {
    message("DEMO ASSUMPTION: fixed contamination fraction = ", rho,
            "; this is not an estimate from the data.")
    channel <- SoupX::setContaminationFraction(channel, rho)
    mode <- "user_supplied_demo_assumption"
  }
  corrected <- SoupX::adjustCounts(channel, roundToInt = TRUE)
  if (!identical(dimnames(corrected), dimnames(toc)) ||
      any(!is.finite(corrected@x)) || any(corrected@x < 0)) {
    stop("Corrected matrix failed alignment or nonnegative-value checks.")
  }

  dir.create(args$output, recursive = TRUE, showWarnings = FALSE)
  Matrix::writeMM(corrected, file.path(args$output, "matrix.mtx"))
  public_features <- data.frame(rownames(corrected), rownames(corrected), "Gene Expression")
  write.table(public_features, file.path(args$output, "features.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(colnames(corrected), file.path(args$output, "barcodes.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  metrics <- data.frame(
    barcode = colnames(toc),
    simulation_population = demo$cells$simulation_population,
    simulation_batch = demo$cells$simulation_batch,
    provisional_cluster = unname(clusters[colnames(toc)]),
    input_umis = as.numeric(Matrix::colSums(toc)),
    corrected_umis = as.numeric(Matrix::colSums(corrected)),
    contamination_fraction = channel$metaData[colnames(toc), "rho"],
    contamination_mode = mode
  )
  write.csv(metrics, file.path(args$output, "cell_metrics.csv"), row.names = FALSE)
  metadata <- data.frame(
    key = c("data_origin", "example", "contamination_mode", "seed", "cells", "features", "roundToInt"),
    value = c(demo$origin, "synthetic ambient-RNA workflow", mode, seed, ncol(toc), nrow(toc), TRUE)
  )
  write.table(metadata, file.path(args$output, "run_metadata.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(args$output, "sessionInfo.txt"))
  message("Wrote corrected count matrix and run metadata to ", args$output)
}

main()
