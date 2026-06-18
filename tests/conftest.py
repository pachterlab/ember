import numpy as np
import pandas as pd
import anndata as ad
import pytest
from scipy.sparse import random as sparse_random, csr_matrix
from scipy import sparse


def make_adata(n_cells=800, n_genes=400, n_blocks=5, density=0.12, seed=0):
    rng = np.random.default_rng(seed)
    X = sparse_random(
        n_cells, n_genes, density=density, format='csr',
        data_rvs=lambda s: rng.exponential(2, size=s)
    ).astype(np.float32)
    labels = rng.choice([f'B{i}' for i in range(n_blocks)], size=n_cells)
    obs = pd.DataFrame(
        {
            'partition': labels,
            'sample_id': [f's{i % 10}' for i in range(n_cells)],
            'category':  rng.choice(['cat_A', 'cat_B'], size=n_cells).tolist(),
            'condition': rng.choice(['ctrl', 'treat'], size=n_cells).tolist(),
        },
        index=[str(i) for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f'g{i}' for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def small_adata():
    return make_adata(n_cells=400, n_genes=200, n_blocks=4, seed=0)


@pytest.fixture
def h5ad_path(tmp_path):
    """Write a small synthetic dataset to a temporary .h5ad file."""
    adata = make_adata(n_cells=400, n_genes=200, n_blocks=4, seed=1)
    path = tmp_path / 'test_data.h5ad'
    adata.write_h5ad(path)
    return str(path)


# ── Reference entropy implementation (replicates upstream sparse-loop logic) ──

def _safe_log2(mat):
    m = mat.copy()
    m.data = np.log2(m.data)
    m.data[~np.isfinite(m.data)] = 0
    return m


def _safe_div(num, den):
    d = np.copy(den)
    d[np.isclose(d, 0)] = np.inf
    r = num.multiply(1.0 / d)
    r.eliminate_zeros()
    return r


def entropy_reference(adata, pl):
    c = adata.X
    if not sparse.issparse(c):
        c = csr_matrix(c)
    c = c.astype(np.float64)
    tot = np.asarray(c.sum(axis=0)).ravel()

    pi = _safe_div(c, tot)
    lpi = _safe_log2(pi)
    ET = -pi.multiply(lpi).sum(axis=0).A1
    ET[np.isclose(tot, 0)] = -1

    blks = adata.obs[pl].unique()
    ng, nb = adata.shape[1], len(blks)
    op = adata.obs[pl].values
    EW = np.zeros(ng)
    PBN = np.zeros((ng, nb))

    for idx, blk in enumerate(blks):
        mask = op == blk
        bc = c[mask, :].astype(np.float64)
        bs = np.asarray(bc.sum(axis=0)).ravel()
        q = _safe_div(bc, bs)
        lq = _safe_log2(q)
        ent = -q.multiply(lq).sum(axis=0).A1
        pcj = np.divide(bs, tot, out=np.zeros_like(bs), where=~np.isclose(tot, 0))
        we = ent * pcj
        EW += we
        PBN[:, idx] = we

    with np.errstate(invalid='ignore', divide='ignore'):
        Psi = np.where(ET > 0, EW / ET, -1)

    EW_safe = np.where(EW > 0, EW, np.inf)
    Pb = PBN / EW_safe[:, None]
    Pb[~np.isfinite(Pb)] = 0
    Pb_df = pd.DataFrame(Pb, index=adata.var.index, columns=blks)

    safe_pb = np.where(Pb > 0, Pb, 1.0)
    log_pb = np.log2(safe_pb)
    log_pb[Pb <= 0] = 0
    Zeta = 1 - (-np.sum(Pb * log_pb, axis=1) / np.log2(nb))

    return Psi, Pb_df, Zeta
