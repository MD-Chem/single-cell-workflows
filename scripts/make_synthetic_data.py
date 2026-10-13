#!/usr/bin/env python3
"""Create artificial counts for a single-cell software demonstration.

No downloads, private files, real genes, or biological labels are used.
Population, batch, and QC-proxy effects are invented for this example.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse

ORIGIN = "independent_simulation_v1"


def simulate(cells: int = 360, genes: int = 600, empty_droplets: int = 180,
             seed: int = 7) -> tuple[ad.AnnData, sparse.csr_matrix, list[str]]:
    """Return cells × genes counts and raw droplets × genes counts."""
    if cells < 30 or genes < 150 or empty_droplets < 20:
        raise ValueError("Use at least 30 cells, 150 genes, and 20 empty droplets.")
    rng = np.random.default_rng(seed)
    populations = np.arange(cells) % 3
    rng.shuffle(populations)
    batches = rng.integers(0, 2, size=cells)
    baseline = rng.gamma(shape=1.5, scale=1.0, size=genes)
    baseline[:15] *= 0.2
    profiles = np.tile(baseline, (cells, 1))
    for population in range(3):
        profiles[populations == population, 20 + 30 * population:50 + 30 * population] *= 7
    profiles[batches == 1, 115:145] *= 1.25
    # The first 15 artificial features provide a QC plotting example only.
    # They are not mitochondrial genes and carry no biological interpretation.
    profiles[:max(3, cells // 40), :15] *= 350
    profiles /= profiles.sum(axis=1, keepdims=True)
    depth = rng.lognormal(np.log(1800), 0.3, size=cells)
    ambient_profile = baseline / baseline.sum()
    endogenous = rng.poisson(depth[:, None] * profiles)
    ambient = rng.poisson(depth[:, None] * 0.06 * ambient_profile)
    counts = sparse.csr_matrix(endogenous + ambient, dtype=np.int32)
    empties = sparse.csr_matrix(rng.poisson(
        rng.uniform(15, 65, size=empty_droplets)[:, None] * ambient_profile), dtype=np.int32)
    barcodes = [f"SIMCELL{i + 1:04d}" for i in range(cells)]
    gene_ids = [f"Gene{i + 1:04d}" for i in range(genes)]
    obs = pd.DataFrame({
        "simulation_population": pd.Categorical([f"population_{x + 1}" for x in populations]),
        "simulation_batch": pd.Categorical([f"batch_{x + 1}" for x in batches]),
    }, index=pd.Index(barcodes, name="barcode"))
    var = pd.DataFrame({"is_qc_proxy": np.arange(genes) < 15}, index=gene_ids)
    data = ad.AnnData(X=counts, obs=obs, var=var)
    data.uns["data_origin"] = ORIGIN
    data.uns["simulation_seed"] = seed
    data.uns["interpretation"] = "Artificial counts and labels; no biological interpretation."
    raw = sparse.vstack([counts, empties], format="csr")
    raw_barcodes = barcodes + [f"EMPTY{i + 1:04d}" for i in range(empty_droplets)]
    return data, raw, raw_barcodes


def write_10x(directory: Path, counts: sparse.spmatrix,
              genes: list[str], barcodes: list[str]) -> None:
    """Write genes × barcodes Matrix Market data using 10x-style filenames."""
    directory.mkdir(parents=True, exist_ok=False)
    io.mmwrite(directory / "matrix.mtx", counts.T.tocoo(), field="integer")
    pd.DataFrame({0: genes, 1: genes, 2: "Gene Expression"}).to_csv(
        directory / "features.tsv", sep="\t", header=False, index=False)
    pd.Series(barcodes).to_csv(directory / "barcodes.tsv", sep="\t", header=False, index=False)


def generate(output: Path, cells: int = 360, genes: int = 600,
             empty_droplets: int = 180, seed: int = 7) -> ad.AnnData:
    output = output.expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Output directory must be empty; choose a new directory.")
    data, raw, raw_barcodes = simulate(cells, genes, empty_droplets, seed)
    output.mkdir(parents=True, exist_ok=True)
    data.write_h5ad(output / "synthetic.h5ad", compression="gzip")
    write_10x(output / "raw", raw, data.var_names.tolist(), raw_barcodes)
    write_10x(output / "filtered", data.X, data.var_names.tolist(), data.obs_names.tolist())
    data.obs.to_csv(output / "synthetic_metadata.csv", index_label="barcode")
    pd.DataFrame({"key": ["data_origin", "seed", "cells", "genes", "empty_droplets"],
                  "value": [ORIGIN, seed, cells, genes, empty_droplets]}).to_csv(
                      output / "provenance.csv", index=False)
    return data


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("generated"))
    parser.add_argument("--cells", type=int, default=360)
    parser.add_argument("--genes", type=int, default=600)
    parser.add_argument("--empty-droplets", type=int, default=180)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()
    data = generate(args.output, args.cells, args.genes, args.empty_droplets, args.seed)
    print(f"Wrote synthetic data: {data.n_obs} cells × {data.n_vars} artificial genes.")


if __name__ == "__main__":
    main()
