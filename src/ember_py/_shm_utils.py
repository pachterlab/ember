"""Shared-memory helpers for zero-copy sparse matrix IPC across worker processes."""
import numpy as np
from multiprocessing import shared_memory
from scipy.sparse import csr_matrix


def create_shm_csr(csr):
    """
    Copy a CSR matrix into three shared-memory blocks (data / indices / indptr).

    Returns
    -------
    handles : list[SharedMemory]
        Owned by the calling process — pass to release_shm() when done.
    descriptor : dict
        Lightweight dict (names + shapes) that workers can pickle cheaply and
        use to re-attach without copying the full matrix through the IPC pipe.
    """
    data    = np.ascontiguousarray(csr.data,    dtype=np.float32)
    indices = np.ascontiguousarray(csr.indices, dtype=np.int32)
    indptr  = np.ascontiguousarray(csr.indptr,  dtype=np.int32)

    shm_d = shared_memory.SharedMemory(create=True, size=max(data.nbytes,    1))
    shm_i = shared_memory.SharedMemory(create=True, size=max(indices.nbytes, 1))
    shm_p = shared_memory.SharedMemory(create=True, size=max(indptr.nbytes,  1))

    if data.nbytes:
        np.ndarray(data.shape,    dtype=np.float32, buffer=shm_d.buf)[:] = data
    if indices.nbytes:
        np.ndarray(indices.shape, dtype=np.int32,   buffer=shm_i.buf)[:] = indices
    np.ndarray(indptr.shape, dtype=np.int32, buffer=shm_p.buf)[:] = indptr

    descriptor = {
        'd': (shm_d.name, data.shape),
        'i': (shm_i.name, indices.shape),
        'p': (shm_p.name, indptr.shape),
        'shape': csr.shape,
    }
    return [shm_d, shm_i, shm_p], descriptor


def read_shm_csr_rows(desc, row_indices):
    """
    Attach to shared memory, extract ``row_indices`` rows, return a new CSR copy.
    Shared-memory handles are closed before returning — single-use worker helper.
    """
    shm_d = shared_memory.SharedMemory(name=desc['d'][0])
    shm_i = shared_memory.SharedMemory(name=desc['i'][0])
    shm_p = shared_memory.SharedMemory(name=desc['p'][0])
    try:
        data    = np.ndarray(desc['d'][1], dtype=np.float32, buffer=shm_d.buf)
        indices = np.ndarray(desc['i'][1], dtype=np.int32,   buffer=shm_i.buf)
        indptr  = np.ndarray(desc['p'][1], dtype=np.int32,   buffer=shm_p.buf)
        full    = csr_matrix((data, indices, indptr), shape=desc['shape'], copy=False)
        return csr_matrix(full[row_indices], dtype=np.float32, copy=True)
    finally:
        shm_d.close()
        shm_i.close()
        shm_p.close()


def release_shm(handles):
    """Close and unlink all owned shared-memory blocks (main process only)."""
    for shm in handles:
        try:
            shm.close()
            shm.unlink()
        except Exception:
            pass
