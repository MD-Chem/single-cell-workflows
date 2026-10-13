# Additional method excerpts

These selected excerpts demonstrate code from the supplied exploratory notebooks. They are presented separately because the sources have mixed or unclear project context; no specific project attribution or biological result is claimed here.

| Excerpt | Implemented method |
| --- | --- |
| [Harmony and subclustering](harmony-subclustering.ipynb) | PCA, Harmony adjustment, neighbors built from the adjusted representation, and downstream UMAP/Leiden analysis |
| [Seurat CCA exploration](seurat-cca-exploration.Rmd) | Seurat layer-based CCA integration, clustering of the integrated representation and comparison with an unintegrated branch |
| [Epithelial downstream analysis](epithelial-downstream-analysis.Rmd) | CytoTRACE2 and Monocle3 conversion, graph learning and pseudotime exploration |

The epithelial excerpt starts from an already integrated Seurat object. It does not implement the upstream integration, and its input alone is not evidence of a particular integration algorithm.

All three files retain substantive method calls. Private inputs, metadata values, selected genes, biological annotations and outputs are omitted. They have not been rerun against the withheld data. These are method examples from research code, distinct from the independently simulated companion in `scripts/`.

[Portfolio](../README.md) · [Code index](../docs/CODE_INDEX.md)
