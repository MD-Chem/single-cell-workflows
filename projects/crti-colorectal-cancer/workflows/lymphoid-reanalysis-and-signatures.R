# Lymphoid reanalysis, exploratory marker review and signature scoring
# Source-derived edited code. Private objects, labels, gene panels and results
# are not distributed; this excerpt has not been rerun. The original exploratory
# branches are represented without reproducing their biological conclusions.
# Retained parameters are analysis settings, not universally suitable defaults.
# Requires a Seurat v5 object with RNA counts, per-cell metadata and enough
# cells/features for the requested PCA dimensions. Outputs remain local/private.

if (!identical(Sys.getenv("RUN_PRIVATE_ANALYSIS"), "yes")) {
  stop("Code-only portfolio excerpt. Configure private inputs before setting RUN_PRIVATE_ANALYSIS=yes.")
}

library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(UCell)
library(ggpubr)
library(Matrix)
library(pheatmap)
set.seed(1234)
input_root <- "data/private"
out_dir <- "outputs/private/lymphoid-analysis"
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

lymphoid <- readRDS(file.path(input_root, "lymphoid_object.rds"))
# sample_merged/sample_type are local metadata fields; their assignments are not bundled.
stopifnot(all(c("sample_merged", "sample_type", "seurat_clusters") %in% names(lymphoid@meta.data)))
Idents(lymphoid) <- "seurat_clusters"

# Recompute normalization/scaling after subsetting.
# Editorial cleanup: use a fresh object from counts rather than mutating internal slots.
lymphoid <- CreateSeuratObject(
  counts = GetAssayData(lymphoid, assay = "RNA", layer = "counts"),
  meta.data = lymphoid@meta.data
)
lymphoid <- NormalizeData(lymphoid)
lymphoid <- FindVariableFeatures(lymphoid, selection.method = "vst", nfeatures = 2000)
lymphoid <- ScaleData(lymphoid)
lymphoid <- RunPCA(lymphoid, npcs = 100)
ElbowPlot(lymphoid, ndims = 50)
DimHeatmap(lymphoid, dims = 10:30, cells = 200, balanced = TRUE)

# Two representative PC-review branches from the source; neither is called optimal.
lymphoid_pc10 <- lymphoid
lymphoid_pc10 <- FindNeighbors(lymphoid_pc10, dims = 1:10, verbose = FALSE)
lymphoid_pc10 <- FindClusters(lymphoid_pc10, verbose = FALSE, resolution = 0.5)
lymphoid_pc10 <- RunUMAP(lymphoid_pc10, dims = 1:10, verbose = FALSE)
DimPlot(lymphoid_pc10, label = TRUE, reduction = "umap", split.by = "sample_merged",
        cols = colors, pt.size = 1)

lymphoid_pc15 <- lymphoid
lymphoid_pc15 <- FindNeighbors(lymphoid_pc15, dims = 1:15, verbose = FALSE)
lymphoid_pc15 <- FindClusters(lymphoid_pc15, verbose = FALSE, resolution = 0.5)
lymphoid_pc15 <- RunUMAP(lymphoid_pc15, dims = 1:15, verbose = FALSE)
DimPlot(lymphoid_pc15, label = TRUE, reduction = "umap", split.by = "sample_merged",
        cols = colors, pt.size = 1)

# Review/exclude clusters locally; no original cluster assignments are encoded.
lymphoid <- readRDS(file.path(input_root, "reviewed_lymphoid_embedding.rds"))
clusters_to_exclude <- readLines(file.path(input_root, "lymphoid_excluded_clusters.txt"))
Idents(lymphoid) <- "seurat_clusters"
if (length(clusters_to_exclude)) {
  excluded <- subset(lymphoid, idents = clusters_to_exclude)
  excluded_barcodes <- rownames(excluded@meta.data)
  write.table(excluded_barcodes, file.path(out_dir, "excluded_barcodes.txt"),
               quote = FALSE, col.names = FALSE, row.names = FALSE)
  lymphoid <- subset(lymphoid, idents = clusters_to_exclude, invert = TRUE)
}
# Revisit normalization, dimensions and clusters after any exclusion before interpretation.

# Marker review by group: this reproduces the source's exploratory cell-level test,
# not a replicate-aware treatment-effect analysis.
sc43 <- lymphoid
Idents(sc43) <- "sample_merged"
sc43.markers <- FindAllMarkers(sc43, min.pct = 0.25, logfc.threshold = 0.25)
write.table(sc43.markers, file.path(out_dir, "group_markers.tsv"),
             quote = FALSE, sep = "	", col.names = TRUE, row.names = FALSE)
sc43.markers <- read_tsv(file.path(out_dir, "group_markers.tsv"))
sc43.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
pdf(file.path(out_dir, "group_marker_heatmap.pdf"), width = 6, height = 5)
DoHeatmap(sc43, features = top20$gene)
dev.off()

# Gene panels and reviewed labels are supplied locally, never copied into this file.
review_features <- readLines(file.path(input_root, "lymphoid_review_features.txt"))
DotPlot(lymphoid, features = review_features, dot.min = 0.1, col.min = 0.1,
        group.by = "seurat_clusters", cluster.idents = TRUE) +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_line(color = "grey", linewidth = 0.01, linetype = 2))
VlnPlot(lymphoid, features = review_features, cols = colors,
        group.by = "sample_merged", pt.size = 0)
FeaturePlot(lymphoid, features = review_features, split.by = "sample_merged",
            cols = c("lightgrey", "red"), order = TRUE)

# Source named-vector / RenameIdents annotation sequence.
Idents(lymphoid) <- "seurat_clusters"
annotation <- read.csv(file.path(input_root, "lymphoid_annotation.csv"), colClasses = "character")
stopifnot(all(c("cluster", "broad_label") %in% names(annotation)))
stopifnot(!anyDuplicated(annotation$cluster), setequal(annotation$cluster, levels(lymphoid)))
broad.cluster.ids <- annotation$broad_label[match(levels(lymphoid), annotation$cluster)]
names(broad.cluster.ids) <- levels(lymphoid)
lymphoid <- RenameIdents(lymphoid, broad.cluster.ids)
lymphoid$broadtype <- Idents(lymphoid)
DimPlot(lymphoid, label = FALSE, reduction = "umap", split.by = "sample_merged",
        cols = colors, pt.size = 2)

# Descriptive composition with the original group-by / percentage idiom.
DF <- lymphoid@meta.data
DF <- DF %>% group_by(sample_merged, broadtype) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(C = sum(Count)) %>% mutate(percent = Count / C * 100)
print(ggplot(DF, aes(x = sample_merged, y = percent, fill = broadtype)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = colors) +
  coord_flip() + theme_classic())
Cellratio <- prop.table(table(lymphoid$broadtype, lymphoid$sample_merged), margin = 2)
Cellratio <- as.data.frame(Cellratio)
pdf(file.path(out_dir, "composition_by_annotation.pdf"), width = 15, height = 5)
print(ggplot(Cellratio, aes(x = Var2, y = Freq)) +
  geom_bar(stat = "identity", fill = "forest green") +
  facet_wrap(~Var1, ncol = 6, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.8)))
dev.off()

# Matrix/metadata exchange for downstream tools: this is not a velocity calculation.
lymphoid$barcode <- colnames(lymphoid)
lymphoid$UMAP_1 <- lymphoid@reductions$umap@cell.embeddings[, 1]
lymphoid$UMAP_2 <- lymphoid@reductions$umap@cell.embeddings[, 2]
write.csv(lymphoid@meta.data, file.path(out_dir, "metadata.csv"), row.names = FALSE)
counts_matrix <- GetAssayData(lymphoid, assay = "RNA", layer = "counts")
writeMM(counts_matrix, file.path(out_dir, "counts.mtx"))
write.csv(lymphoid@reductions[["pca"]]@cell.embeddings,
          file.path(out_dir, "pca.csv"), row.names = TRUE)
write.table(data.frame(gene = rownames(counts_matrix)), file.path(out_dir, "gene_names.tsv"),
             quote = FALSE, row.names = FALSE, col.names = FALSE)

# CytoTRACE (the original package/API used in this source).
# Conversion to a dense matrix can be memory intensive; inspect input size locally.
library(CytoTRACE)
cyto <- as.matrix(GetAssayData(lymphoid, assay = "RNA", layer = "counts"))
results <- CytoTRACE(cyto)
cyto_scores <- as.data.frame(results$CytoTRACE)
lymphoid <- AddMetaData(lymphoid, metadata = cyto_scores, col.name = "CytoTRACE_score")
pdf(file.path(out_dir, "cytotrace_embedding.pdf"), width = 5, height = 4)
FeaturePlot(lymphoid, features = "CytoTRACE_score", pt.size = 1) &
  scale_color_viridis_c(option = "C")
dev.off()
Idents(lymphoid) <- "sample_merged"
pdf(file.path(out_dir, "cytotrace_group_review.pdf"), width = 5, height = 4)
VlnPlot(lymphoid, features = "CytoTRACE_score", cols = colors, pt.size = 0) +
  stat_summary(fun = median, geom = "point", size = 6, colour = "black", shape = 95)
dev.off()
# Scores are exploratory outputs; no interpretation or group result is supplied.

# Gene-panel heatmap, retaining the source pheatmap and custom-gradient pattern.
# The source's direction-specific gene lists are deliberately not reproduced.
scaled_data <- GetAssayData(lymphoid, assay = "RNA", layer = "scale.data")[review_features, ]
annotation_col <- data.frame(Group = lymphoid$sample_merged, row.names = colnames(lymphoid))
pheatmap(scaled_data, color = colorRampPalette(c("navy", "white", "firebrick"))(100),
         breaks = seq(-2, 2, length.out = 101), cluster_rows = FALSE, cluster_cols = TRUE,
         annotation_col = annotation_col, show_rownames = TRUE, show_colnames = FALSE)
DoHeatmap(lymphoid, features = review_features, group.by = "sample_merged",
          slot = "scale.data", disp.min = -2, disp.max = 2) +
  scale_fill_gradientn(colors = c("navy", "white", "firebrick"),
                       limits = c(-2, 2), na.value = "grey")

# UCell rank-based signature scores. Private named lists replace the original
# direction-specific signatures, so their contents and findings are not disclosed.
signatures <- readRDS(file.path(input_root, "lymphoid_signature_sets.rds"))
stopifnot(is.list(signatures), length(signatures) > 0, !is.null(names(signatures)))
metadata_before <- names(lymphoid@meta.data)
lymphoid <- AddModuleScore_UCell(obj = lymphoid, features = signatures, name = "_UCell")
score_columns <- setdiff(names(lymphoid@meta.data), metadata_before)
FeaturePlot(lymphoid, features = score_columns, split.by = "sample_merged",
            cols = c("grey", "red"), pt.size = 0.5)

# Reviewed subpopulation branch: original fixed cluster IDs are private input.
Idents(lymphoid) <- "seurat_clusters"
subset_clusters <- readLines(file.path(input_root, "lymphoid_subset_clusters.txt"))
subset_obj <- subset(lymphoid, idents = subset_clusters)
DimPlot(subset_obj, label = TRUE, reduction = "umap", cols = colors, pt.size = 2)
Idents(subset_obj) <- "sample_merged"
subset.markers <- FindAllMarkers(subset_obj, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(subset.markers, file.path(out_dir, "subset_markers.csv"), row.names = FALSE)
subset.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(subset_obj, features = top10$gene)
VlnPlot(subset_obj, features = review_features, group.by = "sample_merged",
        sort = "increasing", log = TRUE, pt.size = 0)
saveRDS(lymphoid, file.path(out_dir, "lymphoid_reviewed.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
