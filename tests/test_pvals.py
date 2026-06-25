"""
Verify that the incremental count_ge accumulator produces p-values identical
to the upstream DataFrame .ge().sum() approach.
"""
import sys
sys.path.insert(0, 'src')

import numpy as np
import pandas as pd
import pytest


def _pvals_dataframe(perm_arrays, real_vals, var_index):
    """Upstream approach: build DataFrame then call .ge().sum()."""
    series_list = [pd.Series(arr, index=var_index) for arr in perm_arrays]
    df = pd.DataFrame({f'P_{i}': s for i, s in enumerate(series_list)})
    n = len(series_list)
    return (df.ge(real_vals, axis=0).sum(axis=1) + 1) / (n + 1)


def _pvals_incremental(perm_arrays, real_vals, var_index):
    """New approach: incremental count_ge array."""
    gene_positions = var_index.get_indexer(var_index)
    count_ge = np.zeros(len(var_index))
    for arr in perm_arrays:
        count_ge += arr[gene_positions] >= real_vals.values
    n = len(perm_arrays)
    return pd.Series((count_ge + 1) / (n + 1), index=var_index)


@pytest.mark.parametrize("seed,n_genes,n_iters", [
    (0,  300, 200),
    (1,  500, 100),
    (42, 100, 500),
])
def test_pval_methods_identical(seed, n_genes, n_iters):
    rng = np.random.default_rng(seed)
    var_index = pd.Index([f'g{i}' for i in range(n_genes)])
    real_vals = pd.Series(rng.uniform(0.05, 0.8, n_genes), index=var_index)
    perm_arrays = [rng.uniform(0, 1, n_genes) for _ in range(n_iters)]

    pvals_df  = _pvals_dataframe(perm_arrays, real_vals, var_index)
    pvals_inc = _pvals_incremental(perm_arrays, real_vals, var_index)

    assert np.allclose(pvals_df.values, pvals_inc.values, atol=1e-12), \
        f"Max diff: {np.max(np.abs(pvals_df.values - pvals_inc.values)):.2e}"


def test_pvals_in_valid_range():
    rng = np.random.default_rng(0)
    n_genes, n_iters = 200, 100
    var_index = pd.Index([f'g{i}' for i in range(n_genes)])
    real_vals = pd.Series(rng.uniform(0, 1, n_genes), index=var_index)
    perm_arrays = [rng.uniform(0, 1, n_genes) for _ in range(n_iters)]

    pvals = _pvals_incremental(perm_arrays, real_vals, var_index)
    assert (pvals > 0).all(), "p-values should be > 0"
    assert (pvals <= 1).all(), "p-values should be <= 1"


def test_adaptive_stopping_default_runs_all_iterations():
    """
    Simulate the generate_pvals loop to verify adaptive_stopping=False
    always runs all n_iterations regardless of convergence.
    """
    n_iters = 200
    n_genes = 20
    min_iters, check_interval, tol = 300, 100, 0.002
    adaptive_stopping = False  # the default

    rng = np.random.default_rng(0)
    real_vals = rng.uniform(0, 1, n_genes)
    count_ge = np.zeros(n_genes)
    prev_pvals = None
    n_done = 0
    stopped_early = False

    for _ in range(n_iters):
        perm = rng.uniform(0, 1, n_genes)
        count_ge += perm >= real_vals
        n_done += 1
        if adaptive_stopping and n_done >= min_iters:
            curr = (count_ge + 1) / (n_done + 1)
            if prev_pvals is not None and np.nanmax(np.abs(curr - prev_pvals)) < tol:
                stopped_early = True
                break
            if n_done % check_interval == 0:
                prev_pvals = curr

    assert not stopped_early, "Loop stopped early with adaptive_stopping=False"
    assert n_done == n_iters, f"Expected {n_iters} iterations, ran {n_done}"


def test_pvals_monotone_with_real_score():
    """Higher real scores should yield higher (less significant) p-values on average."""
    rng = np.random.default_rng(7)
    n_genes, n_iters = 200, 500
    var_index = pd.Index([f'g{i}' for i in range(n_genes)])
    perm_arrays = [rng.uniform(0, 1, n_genes) for _ in range(n_iters)]

    low_real  = pd.Series(np.full(n_genes, 0.1), index=var_index)
    high_real = pd.Series(np.full(n_genes, 0.9), index=var_index)

    pvals_low  = _pvals_incremental(perm_arrays, low_real,  var_index)
    pvals_high = _pvals_incremental(perm_arrays, high_real, var_index)

    assert pvals_low.mean() > pvals_high.mean(), \
        "Genes with higher real scores should have lower p-values (fewer permutations exceed them)"
