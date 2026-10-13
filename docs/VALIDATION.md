# Validation record

This repository contains selected, edited research-code excerpts. The underlying private inputs are not distributed, and the edited analyses have not been rerun against them.

## Python notebooks

The five Jupyter notebooks under `projects/` were checked for notebook structure and Python syntax. Saved outputs, execution counts, attachments and original notebook metadata were removed. Paths, identifiers, biological panels and study-derived literal values were reviewed against the supplied sources.

These checks establish file structure and review scope. They do not establish numerical correctness, package compatibility or reproducibility of the studies. No new model was trained and no private dataset was analyzed while preparing the portfolio.

## R workflows

The R files under `workflows/` and `projects/` received source review and delimiter/quote checks. R was unavailable in the preparation environment, so these edited files were not parsed by R or executed with their analysis packages. R Markdown analysis chunks are disabled by default.

## Editorial changes

The public copies include documented refactoring and consistency fixes as well as redaction. For example, repeated per-sample Scrublet blocks were combined into a loop; explicit singlet selection and count-layer handling were added or clarified. These changes are not presented as the exact original execution history. Source-to-excerpt mappings remain private because they contain the original study details.

## Figures and reuse

The two selected HCC figures are archival outputs supplied by Mohua Das and shared with her approval. They were not regenerated or numerically revalidated here. No figures from the other two projects are included.

Before adapting an excerpt, review its expected matrix scale, object structure, required metadata and installed package interfaces. No public dependency lockfile reproduces the historical research environments.
