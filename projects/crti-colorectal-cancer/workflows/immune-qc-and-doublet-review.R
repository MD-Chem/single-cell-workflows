# Immune-cell QC, clustering, doublet review and annotation
# Source-derived edited analysis code; private data and outputs are not bundled.
# Sample names, group assignments, biological marker panels, fitted doublet
# parameters and recorded findings have been removed. This file has not been rerun.
# The retained algorithm/QC settings require review for the dataset being analysed.
# The DoubletFinder excerpt must be run on an independent capture, not pooled
# unrelated captures. This is a documented editorial change to input handling.

if (!identical(Sys.getenv("RUN_PRIVATE_ANALYSIS"), "yes")) {
  stop("Code-only portfolio excerpt. Configure private inputs before setting RUN_PRIVATE_ANALYSIS=yes.")
}

library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(DoubletFinder)
set.seed(1234)
input_root <- "data/private"
out_dir <- "outputs/private/immune-workflow"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

colors <- c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b',
            '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
            '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5',
            '#ad494a', '#8c6d31', '#FF34FF', '#763568', '#63FFAC', '#B79762',
            '#004D43', '#002dff', '#997D87', '#5A0007', '#809693', '#6A3A4C',
            '#1B4400', '#4FC601', '#3B5DFF', '#4A3B53', '#FF2F80', '#61615A',
            '#BA0900', '#6B7900', '#00C2A0', '#FFAA92', '#FF90C9', '#B903AA',
            '#D16100', '#DDEFFF', '#000035', '#7B4F4B', '#A1C299', '#300018',
            '#0AA6D8', '#013349', '#00846F', '#372101', '#FFB500', '#C2FFED',
            '#A079BF', '#CC0744', 'yellow', 'green', 'magenta', 'red')

# Load the original-count Seurat object and privately maintained cell metadata.
# Metadata columns: barcode, capture_id, sample_type, sample_merged.
raw_obj <- readRDS(file.path(input_root, "immune_counts.rds"))
cell_metadata <- read.csv(file.path(input_root, "immune_cell_metadata.csv"),
                           row.names = 1, stringsAsFactors = FALSE)
stopifnot(all(colnames(raw_obj) %in% rownames(cell_metadata)))
raw_obj <- AddMetaData(raw_obj, metadata = cell_metadata[colnames(raw_obj), , drop = FALSE])
raw_obj[["percent.mt"]] <- PercentageFeatureSet(raw_obj, pattern = "^MT-|^mt-")

# QC: violin and paired scatter plots before filtering.
pdf(file.path(out_dir, "qc_before_filtering.pdf"), width = 5, height = 4)
VlnPlot(raw_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0)
dev.off()
plot1 <- FeatureScatter(raw_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(raw_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(file.path(out_dir, "qc_scatter.png"), width = 1200, height = 800)
print(plot1 + plot2)
dev.off()

analysis_obj <- subset(raw_obj,
  subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
pdf(file.path(out_dir, "qc_after_filtering.pdf"), width = 5, height = 4)
VlnPlot(analysis_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
dev.off()

# Standard log-normalization sequence retained from the source.
analysis_obj <- NormalizeData(analysis_obj)
analysis_obj <- FindVariableFeatures(analysis_obj, selection.method = "vst", nfeatures = 2000)
analysis_obj <- ScaleData(analysis_obj)
analysis_obj <- RunPCA(analysis_obj, npcs = 100)
ElbowPlot(analysis_obj, ndims = 50)
DimHeatmap(analysis_obj, dims = 19:41, cells = 200, balanced = TRUE)

# Representative dimensionality-review branches from the original notebook.
# These are exploratory alternatives, not a claim that either is optimal.
object_pc20 <- analysis_obj
object_pc20 <- FindNeighbors(object_pc20, dims = 1:20, verbose = FALSE)
object_pc20 <- FindClusters(object_pc20, verbose = FALSE, resolution = 1)
object_pc20 <- RunUMAP(object_pc20, dims = 1:20, verbose = FALSE)
DimPlot(object_pc20, label = TRUE, reduction = "umap", cols = colors, pt.size = 0.2)

object_pc50 <- analysis_obj
object_pc50 <- FindNeighbors(object_pc50, dims = 1:50, verbose = FALSE)
object_pc50 <- FindClusters(object_pc50, verbose = FALSE, resolution = 1)
object_pc50 <- RunUMAP(object_pc50, dims = 1:50, verbose = FALSE)
DimPlot(object_pc50, label = TRUE, reduction = "umap", cols = colors,
        split.by = "sample_type", pt.size = 0.2)
# Select/review a working object locally; no study-specific selection is encoded.
analysis_obj <- readRDS(file.path(input_root, "reviewed_immune_embedding.rds"))

# Annotation: original named-vector / RenameIdents idiom, without label mapping.
Idents(analysis_obj) <- "seurat_clusters"
annotation <- read.csv(file.path(input_root, "immune_annotation.csv"),
                        colClasses = "character")
stopifnot(all(c("cluster", "broad_label") %in% names(annotation)))
stopifnot(!anyDuplicated(annotation$cluster),
          setequal(annotation$cluster, levels(analysis_obj)))
broad.cluster.ids <- annotation$broad_label[match(levels(analysis_obj), annotation$cluster)]
names(broad.cluster.ids) <- levels(analysis_obj)
analysis_obj <- RenameIdents(analysis_obj, broad.cluster.ids)
analysis_obj$broadtype <- Idents(analysis_obj)
DimPlot(analysis_obj, label = FALSE, reduction = "umap", split.by = "sample_type",
        cols = colors, pt.size = 0.2)

# DoubletFinder parameter sweep and homotypic adjustment.
# Run this section per capture. No original pK or expected doublet count is encoded.
doublet_obj <- readRDS(file.path(input_root, "single_capture_preprocessed.rds"))
stopifnot(length(unique(doublet_obj$capture_id)) == 1L)
sweep.res.list <- paramSweep(doublet_obj, PCs = 1:50, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
write.csv(bcmvn, file.path(out_dir, "doublet_parameter_sweep.csv"), row.names = FALSE)

# This local file is prepared after reviewing the sweep and capture loading.
# Fields: selected_pK, expected_doublet_rate. Values are not distributed.
doublet_settings <- readRDS(file.path(input_root, "doublet_settings.rds"))
stopifnot(is.finite(doublet_settings$selected_pK),
          doublet_settings$selected_pK > 0,
          doublet_settings$expected_doublet_rate > 0,
          doublet_settings$expected_doublet_rate < 1)
homotypic.prop <- modelHomotypic(doublet_obj@meta.data$seurat_clusters)
nExp_poi <- round(doublet_settings$expected_doublet_rate * nrow(doublet_obj@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

metadata_before <- names(doublet_obj@meta.data)
doublet_obj <- doubletFinder(doublet_obj, PCs = 1:50, pN = 0.25,
  pK = doublet_settings$selected_pK, nExp = nExp_poi,
  reuse.pANN = NULL, sct = FALSE)
# Find the new classification field instead of hard-coding its data-derived name.
classification_field <- grep("^DF.classifications", setdiff(names(doublet_obj@meta.data),
                                                            metadata_before), value = TRUE)
stopifnot(length(classification_field) == 1)
doublet_obj$doublet_class <- doublet_obj@meta.data[[classification_field]]
DimPlot(doublet_obj, label = TRUE, reduction = "umap", group.by = "doublet_class",
        cols = colors, pt.size = 0.2)
write.csv(data.frame(barcode = colnames(doublet_obj), doublet_class = doublet_obj$doublet_class),
          file.path(out_dir, "capture_doublet_calls.csv"), row.names = FALSE)
# nExp_poi.adj is available for a separately reviewed sensitivity run.

# Combine reviewed per-capture calls locally, then retain singlets in the main object.
doublet_calls <- read.csv(file.path(input_root, "reviewed_doublet_calls.csv"),
                           row.names = 1, stringsAsFactors = FALSE)
stopifnot(all(colnames(analysis_obj) %in% rownames(doublet_calls)))
analysis_obj$doublet_class <- doublet_calls[colnames(analysis_obj), "doublet_class"]
Idents(analysis_obj) <- "doublet_class"
singlets <- subset(analysis_obj, idents = "Singlet")
# Recompute normalization/PCA after filtering before final downstream use.
singlets <- NormalizeData(singlets)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(singlets)
singlets <- RunPCA(singlets, npcs = 100)
singlets <- FindNeighbors(singlets, dims = 1:50, verbose = FALSE)
singlets <- FindClusters(singlets, verbose = FALSE, resolution = 1)
singlets <- RunUMAP(singlets, dims = 1:50, verbose = FALSE)

# Descriptive composition plots; retained-cell ratios are not replicate-level inference.
DF <- singlets@meta.data
DF <- DF %>% group_by(sample_merged, broadtype) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(C = sum(Count)) %>% mutate(percent = Count / C * 100)
print(ggplot(DF, aes(x = sample_merged, y = percent, fill = broadtype)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = colors) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
Cellratio <- prop.table(table(singlets$broadtype, singlets$sample_merged), margin = 2)
Cellratio <- as.data.frame(Cellratio)
print(ggplot(Cellratio, aes(x = Var2, y = Freq)) +
  geom_bar(stat = "identity", fill = "forest green") +
  facet_wrap(~Var1, ncol = 6, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.8)))

# Both cluster-marker and group-marker review were present in the source.
# Cell-level tests across groups are exploratory and do not account for study replicates.
Idents(singlets) <- "seurat_clusters"
cluster.markers <- FindAllMarkers(singlets, only.pos = TRUE, min.pct = 0.25,
                                  logfc.threshold = 0.25)
Idents(singlets) <- "sample_merged"
group.markers <- FindAllMarkers(singlets, min.pct = 0.25, logfc.threshold = 0.25)
write.table(group.markers, file.path(out_dir, "group_markers.tsv"),
             quote = FALSE, sep = "	", row.names = FALSE)
group.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(singlets, features = top20$gene)

# Subset by reviewed clusters, then recover those barcodes from original counts.
# Local RDS is a named list of cluster-ID vectors, with no assignments bundled.
subset_plan <- readRDS(file.path(input_root, "immune_subset_plan.rds"))
Idents(singlets) <- "seurat_clusters"
for (subset_name in names(subset_plan)) {
  selected <- subset(singlets, idents = subset_plan[[subset_name]])
  selected_barcodes <- rownames(selected@meta.data)
  selected_raw <- subset(raw_obj, cells = selected_barcodes)
  saveRDS(selected_raw, file.path(out_dir, paste0(subset_name, "_counts.rds")))
}

# Marker panels are private input; no original gene lists or conclusions remain.
review_features <- readLines(file.path(input_root, "immune_review_features.txt"))
Idents(singlets) <- "broadtype"
DotPlot(singlets, features = review_features, dot.min = 0.1, col.min = 0.1,
        split.by = "sample_merged", cluster.idents = TRUE) +
  theme(axis.text.x = element_text(angle = 90))
FeaturePlot(singlets, features = review_features, cols = c("lightgrey", "red"),
            order = TRUE, split.by = "sample_merged")
VlnPlot(singlets, features = review_features, group.by = "sample_merged", pt.size = 0)

# Reference-based annotation: the original SingleR call is retained.
# Supply an appropriate, independently obtained reference; no original reference is bundled.
library(SingleR)
ref <- readRDS(file.path(input_root, "annotation_reference.rds"))
results <- SingleR(test = as.SingleCellExperiment(singlets, assay = "RNA"),
                  ref = ref, labels = ref$label.main)
singlets$singlr_labels <- results$labels
DimPlot(singlets, reduction = "umap", cols = colors, group.by = "singlr_labels", label = TRUE)
saveRDS(singlets, file.path(out_dir, "immune_reviewed_object.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
