"""
Verify that generate_entropy_metrics (COO/bincount rewrite) produces results
identical to the upstream sparse-loop reference implementation.
"""
import sys
sys.path.insert(0, 'src')

import numpy as np
import pytest
from conftest import make_adata, entropy_reference
from ember_py.generate_entropy_metrics import generate_entropy_metrics


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_psi_matches_reference(seed):
    adata = make_adata(seed=seed)
    Psi_ref, _, _ = entropy_reference(adata, 'partition')
    Psi_new, _, _ = generate_entropy_metrics(adata, 'partition')

    valid = (Psi_ref != -1) & (np.asarray(Psi_new) != -1)
    assert valid.sum() > 0, "No valid genes found"
    assert np.allclose(Psi_ref[valid], np.asarray(Psi_new)[valid], rtol=1e-4, atol=1e-6), \
        f"Psi mismatch (seed={seed}): max_err={np.max(np.abs(Psi_ref[valid] - np.asarray(Psi_new)[valid])):.2e}"


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_zeta_matches_reference(seed):
    adata = make_adata(seed=seed)
    _, _, Zeta_ref = entropy_reference(adata, 'partition')
    Psi_ref, _, _ = entropy_reference(adata, 'partition')
    _, _, Zeta_new = generate_entropy_metrics(adata, 'partition')

    valid = Psi_ref != -1
    assert np.allclose(Zeta_ref[valid], np.asarray(Zeta_new)[valid], rtol=1e-4, atol=1e-6), \
        f"Zeta mismatch (seed={seed}): max_err={np.max(np.abs(Zeta_ref[valid] - np.asarray(Zeta_new)[valid])):.2e}"


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_psi_block_matches_reference(seed):
    adata = make_adata(seed=seed)
    Psi_ref, Pb_ref, _ = entropy_reference(adata, 'partition')
    _, Pb_new, _ = generate_entropy_metrics(adata, 'partition')

    valid = Psi_ref != -1
    common = [c for c in Pb_ref.columns if c in Pb_new.columns]
    assert len(common) > 0, "No common block columns found"
    assert np.allclose(Pb_ref[common].values[valid], Pb_new[common].values[valid], rtol=1e-4, atol=1e-6), \
        f"Psi_block mismatch (seed={seed})"


def test_invalid_genes_marked_negative_one():
    """Genes with zero total counts must receive Psi = -1 (sentinel), not NaN or a score."""
    adata = make_adata(seed=0)
    # Force two genes to all-zero so we know exactly which should be -1.
    adata.X[:, 0] = 0
    adata.X[:, 5] = 0
    Psi, _, _ = generate_entropy_metrics(adata, 'partition')
    assert Psi[0] == -1, "Zero-count gene 0 should have Psi = -1"
    assert Psi[5] == -1, "Zero-count gene 5 should have Psi = -1"
    assert all(p == -1 or 0 <= p <= 1 for p in Psi), \
        "All Psi values should be in [0, 1] or -1"


def test_invalid_gene_set_matches_reference():
    """The exact set of genes marked invalid should agree with the reference."""
    adata = make_adata(seed=0)
    Psi_ref, _, _ = entropy_reference(adata, 'partition')
    Psi_new, _, _ = generate_entropy_metrics(adata, 'partition')
    invalid_ref = set(np.where(np.asarray(Psi_ref) == -1)[0])
    invalid_new = set(np.where(np.asarray(Psi_new) == -1)[0])
    assert invalid_ref == invalid_new, \
        f"Invalid-gene sets differ: ref={len(invalid_ref)}, new={len(invalid_new)}"


def test_does_not_mutate_adata():
    """generate_entropy_metrics must not write any keys into adata.var."""
    adata = make_adata(seed=0)
    var_cols_before = set(adata.var.columns)
    generate_entropy_metrics(adata, 'partition')
    assert set(adata.var.columns) == var_cols_before, \
        f"adata.var was mutated; new keys: {set(adata.var.columns) - var_cols_before}"


def test_psi_block_sums_to_one():
    """Psi_block rows should sum to ~1 for genes with valid Psi."""
    adata = make_adata(seed=0)
    Psi, Pb, _ = generate_entropy_metrics(adata, 'partition')
    valid = np.asarray(Psi) != -1
    row_sums = Pb.values[valid].sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-5), \
        f"Psi_block rows don't sum to 1: min={row_sums.min():.4f}, max={row_sums.max():.4f}"


def test_zeta_in_range():
    adata = make_adata(seed=0)
    Psi, _, Zeta = generate_entropy_metrics(adata, 'partition')
    valid = np.asarray(Psi) != -1
    assert np.all(np.asarray(Zeta)[valid] >= 0), "Zeta should be >= 0"
    assert np.all(np.asarray(Zeta)[valid] <= 1 + 1e-6), "Zeta should be <= 1"


@pytest.mark.parametrize("n_blocks", [2, 5, 10])
def test_varying_block_counts(n_blocks):
    """COO/bincount path should handle any number of blocks without a loop."""
    adata = make_adata(n_blocks=n_blocks, seed=0)
    Psi_ref, Pb_ref, Zeta_ref = entropy_reference(adata, 'partition')
    Psi_new, Pb_new, Zeta_new = generate_entropy_metrics(adata, 'partition')

    valid = (Psi_ref != -1) & (np.asarray(Psi_new) != -1)
    common = [c for c in Pb_ref.columns if c in Pb_new.columns]
    assert np.allclose(Psi_ref[valid], np.asarray(Psi_new)[valid], rtol=1e-4, atol=1e-6)
    assert np.allclose(Zeta_ref[valid], np.asarray(Zeta_new)[valid], rtol=1e-4, atol=1e-6)
    assert np.allclose(Pb_ref[common].values[valid], Pb_new[common].values[valid], rtol=1e-4, atol=1e-6)
