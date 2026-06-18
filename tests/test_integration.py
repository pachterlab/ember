"""
End-to-end integration tests for light_ember and generate_pvals.
Uses small synthetic data and minimal iterations to keep runtime short.
"""
import sys
sys.path.insert(0, 'src')

import os
import numpy as np
import pandas as pd
import pytest
from conftest import make_adata
from ember_py import light_ember, generate_pvals


@pytest.fixture
def h5ad_path(tmp_path):
    adata = make_adata(n_cells=400, n_genes=200, n_blocks=4, seed=1)
    path = tmp_path / 'test.h5ad'
    adata.write_h5ad(path)
    return str(path)


@pytest.fixture
def run_light_ember(h5ad_path, tmp_path):
    """Run light_ember with minimal settings and return the output directory."""
    save_dir = str(tmp_path / 'output')
    light_ember(
        h5ad_dir=h5ad_path,
        partition_label='partition',
        sample_id_col='sample_id',
        category_col='category',
        condition_col='condition',
        save_dir=save_dir,
        num_draws=10,
        n_pval_iterations=50,
        n_cpus=2,
        partition_pvals=True,
    )
    return save_dir


def test_entropy_csv_created(run_light_ember):
    path = os.path.join(run_light_ember, 'entropy_metrics_partition.csv')
    assert os.path.exists(path), "entropy_metrics_partition.csv not created"


def test_entropy_csv_has_expected_columns(run_light_ember):
    path = os.path.join(run_light_ember, 'entropy_metrics_partition.csv')
    df = pd.read_csv(path, index_col=0)
    for col in ['Psi_mean_partition', 'Psi_std_partition', 'Zeta_mean_partition', 'Zeta_std_partition']:
        assert col in df.columns, f"Missing column: {col}"


def test_pval_columns_merged(run_light_ember):
    path = os.path.join(run_light_ember, 'entropy_metrics_partition.csv')
    df = pd.read_csv(path, index_col=0)
    for col in ['Psi p-value', 'Zeta p-value', 'Psi q-value', 'Zeta q-value']:
        assert col in df.columns, f"Missing p-value column: {col}"


def test_psi_values_in_range(run_light_ember):
    path = os.path.join(run_light_ember, 'entropy_metrics_partition.csv')
    df = pd.read_csv(path, index_col=0)
    psi = df['Psi_mean_partition'].dropna()
    assert (psi >= 0).all() and (psi <= 1).all(), \
        f"Psi_mean out of [0,1]: min={psi.min():.3f}, max={psi.max():.3f}"


def test_zeta_values_in_range(run_light_ember):
    path = os.path.join(run_light_ember, 'entropy_metrics_partition.csv')
    df = pd.read_csv(path, index_col=0)
    zeta = df['Zeta_mean_partition'].dropna()
    assert (zeta >= 0).all() and (zeta <= 1 + 1e-6).all(), \
        f"Zeta_mean out of [0,1]: min={zeta.min():.3f}, max={zeta.max():.3f}"


def test_pvals_in_range(run_light_ember):
    path = os.path.join(run_light_ember, 'entropy_metrics_partition.csv')
    df = pd.read_csv(path, index_col=0)
    for col in ['Psi p-value', 'Zeta p-value']:
        vals = df[col].dropna()
        assert (vals > 0).all() and (vals <= 1).all(), \
            f"{col} out of (0,1]: min={vals.min():.3f}, max={vals.max():.3f}"


def test_psi_block_df_created(run_light_ember):
    path = os.path.join(run_light_ember, 'Psi_block_df', 'mean_Psi_block_df_partition.csv')
    assert os.path.exists(path), "mean_Psi_block_df_partition.csv not created"


def test_psi_block_columns_are_partition_values(run_light_ember, h5ad_path):
    import anndata as ad
    adata = ad.read_h5ad(h5ad_path)
    expected_blocks = set(adata.obs['partition'].unique())

    path = os.path.join(run_light_ember, 'Psi_block_df', 'mean_Psi_block_df_partition.csv')
    df = pd.read_csv(path, index_col=0)
    assert set(df.columns) == expected_blocks, \
        f"Psi_block columns {set(df.columns)} don't match partition values {expected_blocks}"


def test_psi_block_rows_sum_to_one(run_light_ember):
    path = os.path.join(run_light_ember, 'Psi_block_df', 'mean_Psi_block_df_partition.csv')
    df = pd.read_csv(path, index_col=0)
    row_sums = df.sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-4), \
        f"Psi_block rows don't sum to 1: min={row_sums.min():.4f}, max={row_sums.max():.4f}"


def test_rerun_does_not_raise(run_light_ember, h5ad_path, tmp_path):
    """Running light_ember twice to the same directory should not raise (re-run safety)."""
    light_ember(
        h5ad_dir=h5ad_path,
        partition_label='partition',
        sample_id_col='sample_id',
        category_col='category',
        condition_col='condition',
        save_dir=run_light_ember,
        num_draws=10,
        n_pval_iterations=50,
        n_cpus=2,
        partition_pvals=True,
    )
