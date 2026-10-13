"""Checks for data orientation, preservation, validation, and privacy boundaries."""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import io, sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from make_synthetic_data import generate, simulate
from scanpy_workflow import load_input, run, validate_counts


def test_generator_reproducible_and_10x_orientation(tmp_path):
    data = generate(tmp_path / "input", cells=60, genes=180, empty_droplets=20, seed=11)
    again, _, _ = simulate(cells=60, genes=180, empty_droplets=20, seed=11)
    assert (data.X != again.X).nnz == 0
    raw = io.mmread(tmp_path / "input/raw/matrix.mtx", spmatrix=True).tocsr()
    filtered = io.mmread(tmp_path / "input/filtered/matrix.mtx", spmatrix=True).tocsr()
    assert raw.shape == (180, 80)
    assert (filtered != data.X.T).nnz == 0
    assert (raw[:, :60] != filtered).nnz == 0
    with pytest.raises(FileExistsError):
        generate(tmp_path / "input")


@pytest.mark.parametrize("value", [-1, np.nan, np.inf, 0.2])
def test_invalid_raw_counts_rejected(value):
    matrix = sparse.csr_matrix([[1, value], [1, 2]])
    with pytest.raises(ValueError):
        validate_counts(matrix)


def test_non_demo_input_rejected_and_hidden_fields_dropped(tmp_path):
    data, _, _ = simulate(cells=60, genes=180, empty_droplets=20)
    data.obs["unapproved_annotation"] = "do not export"
    data.uns["unapproved_note"] = "do not export"
    data.layers["unapproved_layer"] = data.X.copy()
    filename = tmp_path / "synthetic.h5ad"
    data.write_h5ad(filename)
    clean, _, _ = load_input(filename)
    assert "unapproved_annotation" not in clean.obs
    assert "unapproved_note" not in clean.uns
    assert not clean.layers
    data.var_names = ["real_feature"] + data.var_names.tolist()[1:]
    data.write_h5ad(filename)
    with pytest.raises(ValueError, match="artificial"):
        load_input(filename)


def test_end_to_end_preserves_counts(tmp_path):
    data, _, _ = simulate(cells=60, genes=180, empty_droplets=20)
    result = run(data, tmp_path / "output", min_genes=20, n_hvg=100)
    expected = data[result.obs_names, result.var_names].X
    assert (result.layers["rawcounts"] != expected).nnz == 0
    assert np.isfinite(result.obsm["X_umap"]).all()
    saved = ad.read_h5ad(tmp_path / "output/synthetic_processed.h5ad")
    assert sparse.issparse(saved.X)
    assert saved.n_obs < data.n_obs  # deliberately high artificial QC-proxy cells
    assert (tmp_path / "output/synthetic_umap.png").stat().st_size > 1000
    assert "data_origin" in pd.read_csv(tmp_path / "output/provenance.csv")["key"].tolist()


def test_mock_adjusted_matrix_interoperability(tmp_path):
    """Exercise the R output file format, without claiming to execute SoupX."""
    generated = tmp_path / "input"
    original = generate(generated, cells=60, genes=180, empty_droplets=20)
    fixture = generated / "filtered"
    # A deliberate fractional transformation used only to test matrix import.
    mock_adjusted = original.X.astype(float) * 0.875
    io.mmwrite(fixture / "matrix.mtx", mock_adjusted.T)
    provenance = pd.DataFrame({"key": ["data_origin"], "value": ["independent_simulation_v1"]})
    provenance.to_csv(fixture / "run_metadata.tsv", sep="\t", index=False)
    metrics = original.obs.copy()
    metrics["unapproved_column"] = "not exported"
    metrics.to_csv(fixture / "cell_metrics.csv", index_label="barcode")
    loaded, layer, _ = load_input(input_dir=fixture)
    assert layer == "soupx_adjusted_counts"
    assert (loaded.X != mock_adjusted).nnz == 0
    assert "unapproved_column" not in loaded.obs
    assert loaded.obs["simulation_population"].tolist() == original.obs["simulation_population"].tolist()
    result = run(loaded, tmp_path / "adjusted-output", min_genes=20, n_hvg=100, count_layer=layer)
    expected = loaded[result.obs_names, result.var_names].X
    assert (result.layers[layer] != expected).nnz == 0
    assert "rawcounts" not in result.layers
    io.mmwrite(fixture / "matrix.mtx", mock_adjusted.T[:, :-1])
    with pytest.raises(ValueError, match="dimensions"):
        load_input(input_dir=fixture)
    io.mmwrite(fixture / "matrix.mtx", mock_adjusted.T)
    features = pd.read_csv(fixture / "features.tsv", sep="\t", header=None)
    features.iloc[0, 0] = "unapproved_identifier"
    features.to_csv(fixture / "features.tsv", sep="\t", header=False, index=False)
    with pytest.raises(ValueError, match="artificial"):
        load_input(input_dir=fixture)
    provenance["value"] = "unrecognized_origin"
    provenance.to_csv(fixture / "run_metadata.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="provenance"):
        load_input(input_dir=fixture)
    (fixture / "run_metadata.tsv").unlink()
    with pytest.raises(FileNotFoundError):
        load_input(input_dir=fixture)
