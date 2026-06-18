"""
Verify shared-memory CSR round-trips preserve exact values and that the
registry correctly tracks and cleans up blocks.
"""
import sys
sys.path.insert(0, 'src')

import numpy as np
import pytest
from scipy.sparse import csr_matrix, random as sparse_random
from ember_py._shm_utils import create_shm_csr, read_shm_csr_rows, release_shm, \
    cleanup_shm, _reg_read, _reg_write, _REGISTRY
import os


@pytest.fixture
def csr_matrix_fixture():
    rng = np.random.default_rng(5)
    X = sparse_random(1000, 500, density=0.12, format='csr',
                      data_rvs=lambda s: rng.exponential(2, size=s)).astype(np.float32)
    return csr_matrix(X)


def test_full_matrix_round_trip(csr_matrix_fixture):
    csr = csr_matrix_fixture
    handles, desc = create_shm_csr(csr)
    try:
        all_rows = np.arange(csr.shape[0], dtype=np.int64)
        recovered = read_shm_csr_rows(desc, all_rows)
        assert np.allclose(csr.toarray(), recovered.toarray(), atol=1e-6)
    finally:
        release_shm(handles)


def test_partial_row_slice(csr_matrix_fixture):
    csr = csr_matrix_fixture
    handles, desc = create_shm_csr(csr)
    try:
        row_idx = np.array([0, 5, 99, 200, 499, 999], dtype=np.int64)
        partial = read_shm_csr_rows(desc, row_idx)
        assert np.allclose(csr[row_idx].toarray(), partial.toarray(), atol=1e-6)
    finally:
        release_shm(handles)


def test_shape_preserved(csr_matrix_fixture):
    csr = csr_matrix_fixture
    handles, desc = create_shm_csr(csr)
    try:
        all_rows = np.arange(csr.shape[0], dtype=np.int64)
        recovered = read_shm_csr_rows(desc, all_rows)
        assert recovered.shape == csr.shape
    finally:
        release_shm(handles)


def test_dtype_preserved(csr_matrix_fixture):
    csr = csr_matrix_fixture
    handles, desc = create_shm_csr(csr)
    try:
        all_rows = np.arange(csr.shape[0], dtype=np.int64)
        recovered = read_shm_csr_rows(desc, all_rows)
        assert recovered.dtype == np.float32
    finally:
        release_shm(handles)


def test_registry_populated_on_create(csr_matrix_fixture, tmp_path, monkeypatch):
    monkeypatch.setattr('ember_py._shm_utils._REGISTRY',
                        str(tmp_path / 'shm_registry.json'))
    import ember_py._shm_utils as shm_mod
    shm_mod._REGISTRY = str(tmp_path / 'shm_registry.json')

    handles, desc = create_shm_csr(csr_matrix_fixture)
    try:
        names = shm_mod._reg_read()
        assert len(names) == 3, "Registry should contain 3 block names (data, indices, indptr)"
        assert desc['d'][0] in names
        assert desc['i'][0] in names
        assert desc['p'][0] in names
    finally:
        release_shm(handles)


def test_registry_cleared_on_release(csr_matrix_fixture, tmp_path, monkeypatch):
    import ember_py._shm_utils as shm_mod
    shm_mod._REGISTRY = str(tmp_path / 'shm_registry.json')

    handles, desc = create_shm_csr(csr_matrix_fixture)
    release_shm(handles)
    names = shm_mod._reg_read()
    assert len(names) == 0, "Registry should be empty after release"


def test_cleanup_removes_stale_blocks(csr_matrix_fixture, tmp_path, monkeypatch):
    import ember_py._shm_utils as shm_mod
    shm_mod._REGISTRY = str(tmp_path / 'shm_registry.json')

    # Create blocks but don't release through release_shm (simulate a crash)
    handles, desc = create_shm_csr(csr_matrix_fixture)
    for h in handles:
        h.close()  # close without unlinking — simulates process exit without cleanup

    # Registry still has the names
    assert len(shm_mod._reg_read()) == 3

    # cleanup_shm should unlink them and clear the registry
    cleanup_shm()
    assert len(shm_mod._reg_read()) == 0


def test_empty_matrix_does_not_crash():
    """A matrix with no non-zeros (all-zero sparse) should not raise."""
    csr = csr_matrix((10, 20), dtype=np.float32)
    handles, desc = create_shm_csr(csr)
    try:
        rows = np.arange(10, dtype=np.int64)
        recovered = read_shm_csr_rows(desc, rows)
        assert recovered.shape == (10, 20)
        assert recovered.nnz == 0
    finally:
        release_shm(handles)
