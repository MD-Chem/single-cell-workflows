#!/usr/bin/env Rscript
# Generic Seurat example using the repository's synthetic count matrix.
# https://satijalab.org/seurat/articles/sctransform_vignette
# Clusters and ranked artificial genes are demonstration outputs, not cell types.

parse_args <- function(args) {
  options <- list("input-dir" = "generated/filtered", output = "outputs/seurat",
                  doublets = "false")
  if ("--help" %in% args) {
    cat(paste0("Usage: Rscript scripts/seurat_workflow.R ",
               "[--input-dir generated/filtered] [--output outputs/seurat] ",
               "[--doublets true]\n",
               "Optional doublet scoring assumes a 5% demo doublet rate and retains all cells.\n"))
    quit(status = 0)
  }
  if (length(args) %% 2L != 0L) stop("Each option requires a value; see --help.")
  for (pair in seq_len(length(args) %/% 2L)) {
    i <- 2L * pair - 1L
    key <- sub("^--", "", args[[i]])
    if (!startsWith(args[[i]], "--") || !key %in% names(options)) {
      stop("Unknown option: ", args[[i]])
    }
    options[[key]] <- args[[i + 1L]]
  }
  options
}

read_counts <- function(directory) {
  required <- file.path(directory, c("matrix.mtx", "features.tsv", "barcodes.tsv"))
  if (!all(file.exists(required))) stop("Missing matrix.mtx, features.tsv or barcodes.tsv in input.")
  counts <- as(Matrix::readMM(required[[1L]]), "dgCMatrix")
  features <- read.delim(required[[2L]], header = FALSE, stringsAsFactors = FALSE,
                         check.names = FALSE, quote = "", comment.char = "")[[1L]]
  barcodes <- read.delim(required[[3L]], header = FALSE, stringsAsFactors = FALSE,
                         check.names = FALSE, quote = "", comment.char = "")[[1L]]
  features <- as.character(features)
  barcodes <- as.character(barcodes)
  if (nrow(counts) != length(features) || ncol(counts) != length(barcodes)) {
    stop("Matrix dimensions do not match feature and barcode tables.")
  }
  if (anyNA(features) || anyNA(barcodes) || any(!nzchar(features)) ||
      any(!nzchar(barcodes)) || anyDuplicated(features) || anyDuplicated(barcodes)) {
    stop("Feature IDs and barcodes must be nonempty and unique.")
  }
  if (any(!is.finite(counts@x)) || any(counts@x < 0) ||
      any(abs(counts@x - round(counts@x)) > 1e-8)) {
    stop("Input counts must be finite, nonnegative integers.")
  }
  if (ncol(counts) < 20L || sum(Matrix::rowSums(counts > 0) >= 5L) < 20L) {
    stop("This example needs at least 20 cells and 20 genes detected in at least 5 cells.")
  }
  if (any(Matrix::colSums(counts) == 0)) stop("Cells must have positive count totals.")
  dimnames(counts) <- list(features, barcodes)
  if (!all(grepl("^Gene[0-9]{4,}$", features)) ||
      !all(grepl("^SIMCELL[0-9]{4,}$", barcodes))) {
    stop("This public example accepts artificial Gene/SIMCELL IDs only.")
  }
  counts
}

read_demo_metadata <- function(directory, cell_ids) {
  # Accept either the generator's filtered/ directory or this repo's SoupX output.
  directory <- normalizePath(directory, mustWork = TRUE)
  origin_csv <- file.path(dirname(directory), "provenance.csv")
  origin_tsv <- file.path(directory, "run_metadata.tsv")
  if (file.exists(origin_csv)) {
    provenance <- read.csv(origin_csv, stringsAsFactors = FALSE)
    metadata_path <- file.path(dirname(directory), "synthetic_metadata.csv")
  } else if (file.exists(origin_tsv)) {
    provenance <- read.delim(origin_tsv, stringsAsFactors = FALSE)
    metadata_path <- file.path(directory, "cell_metrics.csv")
  } else {
    stop("Synthetic-data provenance is required; run the repository's generator first.")
  }
  if (!all(c("key", "value") %in% names(provenance))) stop("Invalid provenance table.")
  origin <- provenance$value[provenance$key == "data_origin"]
  if (length(origin) != 1L || is.na(origin) || origin != "independent_simulation_v1") {
    stop("Expected the repository's independent synthetic-data provenance marker.")
  }
  if (!file.exists(metadata_path)) stop("Synthetic cell metadata is missing.")
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
  metadata <- metadata[match(cell_ids, metadata$barcode), , drop = FALSE]
  rownames(metadata) <- metadata$barcode
  list(origin = origin, cells = metadata)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!args$doublets %in% c("true", "false")) stop("--doublets must be true or false.")
  for (package in c("Matrix", "Seurat", "ggplot2")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Missing R package ", package, "; run scripts/install_r_packages.R first.")
    }
  }
  counts <- read_counts(args[["input-dir"]])
  demo <- read_demo_metadata(args[["input-dir"]], colnames(counts))
  if (dir.exists(args$output) && length(list.files(args$output, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory is nonempty; choose a new --output directory.")
  }
  seed <- 42L
  set.seed(seed)
  object <- Seurat::CreateSeuratObject(counts, project = "synthetic_demo", meta.data = demo$cells)
  if (args$doublets == "true") {
    for (package in c("SingleCellExperiment", "scDblFinder", "BiocParallel")) {
      if (!requireNamespace(package, quietly = TRUE)) {
        stop("For optional doublet scoring, run scripts/install_r_packages.R --with-doublets.")
      }
    }
    # The generator has no doublet truth labels. This exercises the API only,
    # with a fixed illustrative prior; it does not benchmark detection accuracy.
    message("DEMO ASSUMPTION: doublet prior = 0.05. All cells are retained.")
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    sce <- scDblFinder::scDblFinder(
      sce, samples = demo$cells$simulation_batch, dbr = 0.05,
      BPPARAM = BiocParallel::SerialParam(RNGseed = seed)
    )
    object$demo_doublet_score <- sce$scDblFinder.score
    object$demo_doublet_class <- sce$scDblFinder.class
  }
  # RNA retains the supplied counts; SCT is a separate assay.
  object <- Seurat::SCTransform(
    object, assay = "RNA", new.assay.name = "SCT", vst.flavor = "v2",
    variable.features.n = min(2000L, nrow(object)), seed.use = seed, verbose = FALSE
  )
  npcs <- min(30L, ncol(object) - 1L, length(Seurat::VariableFeatures(object)) - 1L)
  if (npcs < 2L) stop("Too few variable features after SCTransform.")
  object <- Seurat::RunPCA(object, assay = "SCT", npcs = npcs,
                           seed.use = seed, verbose = FALSE)
  dims <- seq_len(ncol(Seurat::Embeddings(object, "pca")))
  if (length(dims) < 2L) stop("PCA produced fewer than two components.")
  object <- Seurat::FindNeighbors(
    object, dims = dims, k.param = min(20L, ncol(object) - 1L), verbose = FALSE
  )
  object <- Seurat::FindClusters(object, resolution = 0.5,
                                random.seed = seed, verbose = FALSE)
  object <- Seurat::RunUMAP(
    object, dims = dims, n.neighbors = min(30L, ncol(object) - 1L),
    seed.use = seed, verbose = FALSE
  )
  # Rank genes on log-normalized RNA, not on scaled residuals or a UMAP embedding.
  # This is exploratory cluster review; cells are not independent study replicates.
  object <- Seurat::NormalizeData(object, assay = "RNA", verbose = FALSE)
  if (length(unique(Seurat::Idents(object))) < 2L) {
    stop("Clustering produced one group; cluster-marker ranking requires at least two.")
  }
  markers <- Seurat::FindAllMarkers(
    object, assay = "RNA", slot = "data", only.pos = TRUE,
    test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE
  )
  if (!nrow(markers)) warning("No markers passed the demonstration thresholds.")

  dir.create(args$output, recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, file.path(args$output, "seurat_demo.rds"))
  write.csv(markers, file.path(args$output, "ranked_artificial_genes.csv"), row.names = FALSE)
  metadata <- object[[]]
  metadata$barcode <- rownames(metadata)
  write.csv(metadata, file.path(args$output, "cell_metrics.csv"), row.names = FALSE)
  embedding <- data.frame(barcode = colnames(object), Seurat::Embeddings(object, "umap"),
                          cluster = as.character(Seurat::Idents(object)))
  write.csv(embedding, file.path(args$output, "umap_coordinates.csv"), row.names = FALSE)
  umap <- Seurat::DimPlot(object, reduction = "umap", label = TRUE) +
    ggplot2::ggtitle("Synthetic example: unsupervised clusters")
  ggplot2::ggsave(file.path(args$output, "synthetic_umap.png"), umap,
                  width = 7, height = 5, dpi = 180)
  qc <- ggplot2::ggplot(metadata, ggplot2::aes(x = nCount_RNA, y = nFeature_RNA)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) + ggplot2::theme_minimal() +
    ggplot2::labs(title = "Synthetic example: count-matrix QC",
                   x = "Counts per cell", y = "Detected features per cell")
  ggplot2::ggsave(file.path(args$output, "synthetic_qc.png"), qc,
                  width = 7, height = 5, dpi = 180)
  run_metadata <- data.frame(
    key = c("data_origin", "example", "seed", "cells", "features", "pca_dimensions", "marker_assay",
            "doublet_scoring", "doublet_prior"),
    value = c(demo$origin, "synthetic Seurat workflow", seed, ncol(object), nrow(counts),
              length(dims), "RNA log-normalized data; exploratory cluster review",
              args$doublets, if (args$doublets == "true") "0.05 assumed for demo" else "not run")
  )
  write.table(run_metadata, file.path(args$output, "run_metadata.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(args$output, "sessionInfo.txt"))
  message("Wrote synthetic Seurat example outputs to ", args$output)
}

main()
