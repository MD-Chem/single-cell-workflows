# Characterisation of endothelial cells in progressive fatty liver disease in humans

**Human single-cell/single-nucleus analysis · Python and R · Manuscript in preparation**

For this project, I worked with endothelial-cell expression data in Python and R. These excerpts show how I handled sparse matrices, trained a scVI model, explored PCA-based clustering, and carried embeddings and expression data into Seurat for annotation and downstream review.

## Python notebooks

| Notebook | What the code demonstrates |
| --- | --- |
| [Endothelial subset: scVI and PCA](notebooks/endothelial-subset-scvi-and-pca.ipynb) | Matrix Market/AnnData input, count and log-normalized layers, scVI setup and training, latent-space neighbors, Leiden/UMAP, a separate PCA route, marker-ranking recipes and expression summaries |
| [Endothelial PCA analysis](notebooks/endothelial-pca-analysis.ipynb) | Loading a precomputed AnnData object, normalization, highly variable genes, PCA, neighbor graphs, Leiden/UMAP, feature summaries and ranked-feature plotting |

The scVI excerpt contains the actual `SCVI.setup_anndata`, model construction, `train()` and `get_latent_representation()` calls. It uses a count layer and retains the source model's two hidden layers, 30 latent dimensions and negative-binomial likelihood. **The source setup supplies no batch key**, so this example demonstrates latent modelling rather than establishing batch correction. The second notebook implements PCA analysis and does not train scVI.

The scVI and PCA routes overwrite the active graph and cluster labels in sequence, as in the exploratory source; their outputs should be tracked separately when adapting the code. The notebooks describe which branch later steps use.

## R workflow

[scVI handoff, annotation and subset review](workflows/01_scvi_annotation.Rmd) starts from a Seurat object with an imported scVI reduction. It demonstrates:

- Building a neighbor graph and Leiden clusters from the imported embeddings.
- Normalization, UCell signature scoring and marker-panel plotting.
- Applying a locally reviewed annotation map and examining subsets.
- PCA sensitivity analysis, composition summaries and exploratory marker contrasts.
- GO/KEGG enrichment with clusterProfiler and Matrix Market export for Python.

This R excerpt consumes an existing embedding; model fitting is shown in the first Python notebook. The original biological labels, signature genes and enrichment results are withheld.

## Inputs and execution

Each file lists its required private inputs and dependencies. The Jupyter notebooks use a disabled execution flag, and R Markdown chunks have `eval=FALSE`. The optional Enrichr example has a separate disabled flag because it queries an external service. No gene list was submitted while preparing this portfolio.

These are edited research-code excerpts, not an executable public release of the study. No data, rendered figures, model weights, annotation assignments or result tables are included. Notebook structure and Python syntax were checked; the analyses were not rerun against the unpublished inputs.

[All three projects](../../README.md) · [Code index](../../docs/CODE_INDEX.md)
