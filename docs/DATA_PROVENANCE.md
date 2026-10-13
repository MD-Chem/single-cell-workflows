# Data and code provenance

## Research notebook excerpts

The three `.Rmd` files in `workflows/` are edited excerpts from Mohua Das's own R analysis notebooks. They retain the source workflows and representative code patterns for SoupX correction, doublet checks, normalization, clustering, marker review, annotation, and visualization.

Editorial changes include neutral object names and local input paths, removal of embedded output and study-specific metadata, and small code corrections documented in the notebooks. Repeated exploratory sections have been shortened. The original marker panels, annotation assignments, sample-to-condition mapping, biological conclusions, and recorded dataset dimensions are not distributed.

These excerpts are not executable reproductions of the underlying study. They require private input matrices, Seurat objects, or user-supplied annotation tables. Chunks are disabled by default, and the edited code has not been rerun with those inputs. No raw or processed research data, presentation files, or real analysis figures are included in this release.

## Separate synthetic examples

The files in `scripts/`, the Python tests, and `assets/synthetic_umap.png` are new portfolio demonstrations. They are separate from the research notebook excerpts.

`scripts/make_synthetic_data.py` creates counts entirely from a NumPy random-number generator. It does not download data or read external expression matrices, presentations, marker panels, or research results. Its default dimensions and parameters are arbitrary demo settings.

Features use `Gene0001`-style names. Cells use `SIMCELL0001`-style names; background droplets use `EMPTY0001`-style names. Population and batch labels are invented. The first 15 features form an artificial QC proxy and have no biological identities.

The generator writes the origin tag `independent_simulation_v1`, its seed and dimensions. The synthetic readers validate artificial identifiers and permitted metadata. Those checks reduce accidental input mixing; a tag is not proof that an arbitrary matrix is synthetic. Renaming real data would not make it safe to publish.

## Public files and local outputs

Research inputs belong in the ignored `data/private/` directory. Routine outputs, generated matrices, notebook HTML, and serialized analysis objects stay local. The committed PNG is the independently simulated Scanpy preview.

Seeds and package versions are recorded by the synthetic scripts. The Python requirements pin direct dependencies from the executed example. The R installer is not a lockfile, and the research notebooks do not have a validated public execution environment.
