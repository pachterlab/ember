"""Shared-memory helpers for zero-copy sparse matrix IPC across worker processes."""
import json
import os
import numpy as np
from multiprocessing import shared_memory
from scipy.sparse import csr_matrix

_REGISTRY = os.path.expanduser('~/.ember/shm_registry.json')


def _reg_read():
    try:
        with open(_REGISTRY) as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return []


def _reg_write(names):
    os.makedirs(os.path.dirname(_REGISTRY), exist_ok=True)
    with open(_REGISTRY, 'w') as f:
        json.dump(names, f)


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

    # SharedMemory(size=0) raises ValueError on Linux; 1-byte minimum is safe.
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

    # Persist block names before any worker attaches so cleanup_shm() can
    # recover them if the main process is killed before release_shm() runs.
    _reg_write(_reg_read() + [shm_d.name, shm_i.name, shm_p.name])

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
        # copy=True: the slice must own its data before the shm handles are closed below.
        return csr_matrix(full[row_indices], dtype=np.float32, copy=True)
    finally:
        shm_d.close()
        shm_i.close()
        shm_p.close()


def release_shm(handles):
    """Close and unlink all owned shared-memory blocks (main process only)."""
    released = set()
    for shm in handles:
        try:
            released.add(shm.name)
            shm.close()
            shm.unlink()
        except Exception:
            pass
    _reg_write([n for n in _reg_read() if n not in released])


def cleanup_shm():
    """
    Attach to and unlink any shared memory blocks left by crashed ember runs.
    Safe to call at any time; blocks already released are silently skipped.
    """
    names = _reg_read()
    if not names:
        print("No stale ember shared memory blocks found.")
        return

    n_cleaned = 0
    for name in names:
        try:
            shm = shared_memory.SharedMemory(name=name)
            shm.close()
            shm.unlink()
            n_cleaned += 1
        except FileNotFoundError:
            pass

    _reg_write([])
    print(f"Cleaned up {n_cleaned} stale shared memory block(s).")
