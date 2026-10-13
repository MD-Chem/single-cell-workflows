# Research code index

The portfolio contains selected, edited excerpts from 13 analysis files: three HCC notebooks and 10 additional Python/R files. Repeated exploratory blocks and study-specific content have been omitted. A figure-only PDF supplied with the analysis material is not part of this code release.

## Anti-PD-L1 resistance in hepatocellular carcinoma

| File | Primary methods |
| --- | --- |
| [Ambient RNA and doublets](../workflows/ambient-rna-and-doublets.Rmd) | SoupX, Seurat, scDblFinder |
| [Normalization and clustering](../workflows/sctransform-and-clustering.Rmd) | LogNormalize, SCTransform, PCA, neighbor graphs, UMAP |
| [Annotation and composition](../workflows/broad-annotation.Rmd) | Marker ranking, annotation maps, heatmaps, composition summaries |

## CRTi colorectal-cancer project

| File | Primary methods |
| --- | --- |
| [Python preprocessing](../projects/crti-colorectal-cancer/notebooks/01-preprocessing.ipynb) | Scanpy, AnnData, 10x input, QC and barcode metadata |
| [Python filtering / Scrublet](../projects/crti-colorectal-cancer/notebooks/02-filtering-scrublet.ipynb) | QC review, per-capture doublet scoring and prediction alignment |
| [Python clustering](../projects/crti-colorectal-cancer/notebooks/03-clustering.ipynb) | Normalization, HVGs, PCA, UMAP, Leiden and marker-panel plotting |
| [Immune QC / doublet review](../projects/crti-colorectal-cancer/workflows/immune-qc-and-doublet-review.R) | Seurat, DoubletFinder and SingleR |
| [Lymphoid reanalysis](../projects/crti-colorectal-cancer/workflows/lymphoid-reanalysis-and-signatures.R) | UCell, CytoTRACE, markers, composition and interchange exports |
| [Myeloid processing](../projects/crti-colorectal-cancer/workflows/myeloid-processing.Rmd) | Seurat, SingleR, reviewed subsetting and CytoTRACE |
| [Supporting bulk analysis](../projects/crti-colorectal-cancer/workflows/bulk-rna-differential-expression.R) | tximeta, DESeq2, GLM-PCA, apeglm and fgsea; upstream credits retained |

## Human fatty-liver endothelial-cell project

| File | Primary methods |
| --- | --- |
| [Subset scVI / PCA notebook](../projects/fatty-liver-endothelial-cells/notebooks/endothelial-subset-scvi-and-pca.ipynb) | Sparse input, count layers, scVI model training and latent graphs; PCA comparison |
| [PCA analysis notebook](../projects/fatty-liver-endothelial-cells/notebooks/endothelial-pca-analysis.ipynb) | Scanpy PCA, Leiden, UMAP, summaries and ranking |
| [R annotation workflow](../projects/fatty-liver-endothelial-cells/workflows/01_scvi_annotation.Rmd) | Imported scVI embedding, Seurat, UCell, subset analysis and clusterProfiler |

## Separate executable examples

The programs in [`scripts/`](../scripts/) and tests in [`tests/`](../tests/) were added as synthetic portfolio demonstrations. They are not source files from the three studies. See [the companion guide](SYNTHETIC_EXAMPLE.md) and [validation record](VALIDATION.md).
