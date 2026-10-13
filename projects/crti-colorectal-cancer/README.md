# Reversing immune tolerance in colorectal cancer using designed beta-catenin–TCF inhibitors (CRTi)

**Single-cell and supporting bulk transcriptomics · Python and R · Manuscript in preparation**

This project connects my interests in chemical biology and cancer immunology. I used Python and R to work through single-cell preprocessing, doublet review, clustering and downstream analysis. These excerpts make that computational work visible while the manuscript is in preparation.

## Python notebooks

| Order | Notebook | What to inspect |
| --- | --- | --- |
| 1 | [Count-matrix preprocessing](notebooks/01-preprocessing.ipynb) | 10x import, AnnData, barcode metadata, feature/cell filtering and mitochondrial QC |
| 2 | [Filtering and Scrublet](notebooks/02-filtering-scrublet.ipynb) | QC plots, per-capture doublet scoring and barcode-aligned predictions |
| 3 | [Normalization and clustering](notebooks/03-clustering.ipynb) | Count-layer preservation, normalization, variable genes, PCA, neighborhood sensitivity, UMAP and Leiden |

These three notebooks form a selected sequence. Repeated sample-specific blocks have been parameterized, and the public copy explicitly selects singlets before normalization. Run paths are relative to the repository root. The actual sample assignments and marker panels are private inputs.

## R workflows

The R files are complementary analysis excerpts; they are not a single pipeline to run in filename order. Each documents the object or configuration it expects.

| Workflow | Methods demonstrated |
| --- | --- |
| [Immune QC and doublet review](workflows/immune-qc-and-doublet-review.R) | Seurat preprocessing and clustering, DoubletFinder and SingleR reference annotation |
| [Lymphoid reanalysis and signatures](workflows/lymphoid-reanalysis-and-signatures.R) | Subset analysis, UCell signature scores, CytoTRACE, marker and composition plotting, matrix/metadata export |
| [Myeloid processing](workflows/myeloid-processing.Rmd) | Seurat re-normalization and clustering, reference annotation, marker review, reviewed subsetting and CytoTRACE |

## Supporting bulk RNA-seq analysis

[Bulk differential-expression workflow](workflows/bulk-rna-differential-expression.R) shows Salmon quantification import with tximeta, DESeq2, sample-distance diagnostics, GLM-PCA, effect-size shrinkage and fgsea. This is a supporting bulk-transcriptomics workflow, separate from the single-cell analyses. The source author's and tutorial's credits are preserved in the script; its inclusion does not imply sole authorship.

## Scope of this release

The code retains analysis calls and plotting patterns but includes no new figures, notebook outputs or biological findings. Private inputs supply marker panels, reference objects, annotation maps and experimental contrasts. No counts, sample metadata or result tables are included.

The excerpts received structural and source review, not execution against the unpublished data. The R workflows include both historical and newer package interfaces; inspect the notebook-specific requirements before reuse. A validated public execution environment is not provided for these research files.

[All three projects](../../README.md) · [Code index](../../docs/CODE_INDEX.md) · [Validation](../../docs/VALIDATION.md)
