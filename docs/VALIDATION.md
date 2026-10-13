# Validation record

The executed checks below apply to the separate synthetic companion in `scripts/` and `tests/`. They do not validate the research excerpts in `workflows/`.

## Executed Python example

The default generator and Scanpy command were run with Python 3.12.14 and the direct dependency versions listed in `requirements.txt`.

| Check | Observed result |
| --- | --- |
| Generated input | 360 artificial cells, 600 artificial features, 180 background droplets |
| QC output | 351 cells retained; 600 features retained |
| Graph clustering | 3 demonstration clusters |
| Saved artifacts | AnnData, CSV tables, provenance, QC and UMAP PNGs |
| Automated checks | 8 tests passed |

These are results of an independent software simulation. They do not describe a research dataset.

The eight tests cover deterministic count generation and Matrix Market orientation; rejection of negative, nonfinite, and fractional raw counts; artificial-identifier checks and removal of unknown metadata; an end-to-end raw-count workflow; and a mocked SoupX-format export with fractional adjusted counts.

The interoperability test checks parsing, identifier alignment, provenance, invalid dimensions, and preservation of the adjusted-count layer through the Python analysis. It **does not execute SoupX** or establish correctness of ambient-RNA estimation. The preview in `assets/synthetic_umap.png` is from the raw synthetic Python run.

## R validation boundary

The synthetic R scripts received source review and delimiter/quote checks. The research R notebooks in `workflows/` are selected, edited excerpts and have not been rerun with their required private inputs. R was unavailable in the preparation environment, so the scripts were not parsed by R or executed with Seurat, SoupX, or scDblFinder. Their numerical output and package compatibility remain unverified.

When run, each R example records its parameters and `sessionInfo()`. The installer selects available packages; it does not pin a tested R environment. No R-generated result is committed or presented as verified.

## Re-run checks

From the repository root, with the Python dependencies installed:

```bash
python -m pytest -q tests/test_python.py
```

To reproduce the preview, run the generator and Scanpy commands in the README using fresh output directories. Small numerical or layout differences may occur across platforms and transitive dependency versions.

