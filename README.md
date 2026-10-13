# My single-cell analysis workflow

**Mohua Das, PhD · R and Python · Single-cell transcriptomics**

I use single-cell analysis to connect experimental work with changes in gene expression and cell states. This repository brings together selected excerpts from R notebooks I wrote for my research, covering ambient-RNA correction, doublet checks, normalization, clustering, and annotation.

The notebooks retain my analysis sequence and plotting patterns. I have generalized file paths and identifiers and removed study-specific annotations and embedded outputs for public sharing as the project manuscript is still is under preparation. Research data are not distributed with the code.

## Selected plots from my analysis

These are real analysis figures showing numbered clusters. The biological annotations and experimental labels are not shared. The figures are archived outputs; the public notebook excerpts have not been rerun to reproduce them.

| Clustering overview | Subset reclustering |
| --- | --- |
| ![UMAP overview with numbered clusters](assets/research/clustering-overview.png) | ![UMAP of a subset with numbered clusters](assets/research/subset-reclustering.png) |

For the normalization, clustering, and `DimPlot` code, see my [clustering notebook](workflows/sctransform-and-clustering.Rmd). The [annotation notebook](workflows/broad-annotation.Rmd) shows how I work from cluster markers to reviewed labels.

## Start with my R notebooks

| Notebook | What I worked on |
| --- | --- |
| [Ambient RNA and doublet checks](workflows/ambient-rna-and-doublets.Rmd) | Loading raw and filtered droplet matrices, running SoupX with provisional Seurat clusters, combining corrected samples, and checking doublets |
| [Normalization and clustering](workflows/sctransform-and-clustering.Rmd) | Comparing log normalization with SCTransform, reviewing PCs, building neighbor graphs, and visualizing clusters with UMAP |
| [Broad annotation and composition](workflows/broad-annotation.Rmd) | Ranking cluster markers, making heatmaps, assigning reviewed labels, plotting cell proportions, and preparing subsets for further analysis |

These are edited excerpts of my analysis code. They require local count matrices or Seurat objects and are intended to be read alongside their input notes. They do not include the underlying experiment or its biological conclusions.

## How I approach the analysis

**Start with the count matrices.** My SoupX workflow uses both raw droplets and the filtered matrix. I generate provisional clusters in Seurat, pass those clusters to SoupX, and export the corrected matrix for downstream work.

**Review preprocessing choices.** My notebooks include doublet checks and both log-normalization and SCTransform branches. I inspect PCA and UMAP views while working through dimensionality and clustering choices. The parameters shown in the excerpts need reassessment for a different input.

**Use markers to review clusters.** I rank features within clusters, inspect heatmaps, then apply a reviewed annotation map. I also use grouped summaries and proportion plots to examine the composition of the resulting object. The public notebook accepts an external annotation map because the original biological labels are not part of this release.

## Python / Scanpy companion

I have also included a separate, runnable Python companion for QC, normalization, PCA, UMAP, Leiden clustering, and feature ranking. This companion was added for the portfolio and uses independently simulated cells and genes. It is separate from the research notebooks above.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/make_synthetic_data.py --output generated
python scripts/scanpy_workflow.py --input generated/synthetic.h5ad --output outputs/scanpy
python -m pytest -q
```

Use Python 3.12. On Windows, activate the environment with `.venv\Scripts\Activate.ps1` in PowerShell. Output directories must be empty; choose a new `--output` path for another run.

The companion demonstrates AnnData, sparse matrices, pandas summaries, count-layer preservation, and reusable command-line scripts. Its eight Python checks passed. [Full companion instructions](docs/SYNTHETIC_EXAMPLE.md) include the separate synthetic SoupX and Seurat examples.

<details>
<summary>View the synthetic Scanpy example</summary>

![Synthetic Scanpy example, separate from the research notebooks](assets/synthetic_umap.png)

Every cell, gene, and label in this preview is artificial. It does not show a research result.

</details>

## Reading and using the code

- The R notebooks contain no embedded execution output. Their chunks are disabled by default, and the edited excerpts have not been rerun against the private inputs.
- [Data and code provenance](docs/DATA_PROVENANCE.md) explains which files come from my notebooks and which are new examples.
- [Validation notes](docs/VALIDATION.md) distinguish the executed Python checks from the R code review.
- The original research data, presentation, sample metadata, and biological findings are not included.

## About me

My background combines medicinal chemistry, chemical biology, and cancer genomics. I am interested in research that connects computational analysis with experimental biology.

[GitHub profile](https://github.com/MD-Chem) · [ORCID](https://orcid.org/0000-0002-8891-6333) · [LinkedIn](https://www.linkedin.com/in/dmohua/)

