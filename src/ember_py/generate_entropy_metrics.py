import pandas as pd
import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix


def safe_log2_sparse(mat):
    """
    Helper function applies an element-wise log2 transform only to the non-zero
    *stored* entries in a SciPy CSR sparse matrix. Zeros in a sparse matrix are
    implicit (not stored in `.data`) and remain zero in the output.

    Any non-finite results from the log transform (e.g., `-inf` from explicit
    stored zeros, or `NaN` from negative values) are replaced with 0.
    Parameters
    """
    log_mat = mat.copy()
    log_mat.data = np.log2(log_mat.data)
    log_mat.data[~np.isfinite(log_mat.data)] = 0
    return log_mat


def safe_divide_sparse(numerator, denominator):
    """
    Helper frunction performs element-wise division of a matrix by a vector (or array) element-wise while
    avoiding division-by-zero issues. Supports both SciPy sparse matrices and
    dense NumPy arrays.

    Sparse case:
        - Zeros in the denominator are replaced with `inf` so the result is 0.
        - Zero entries resulting from division are removed from the sparse representation.

    Dense case:
        - Division is performed with suppressed warnings.
        - Positions with zero denominator are set to 0.
    """
    if sparse.issparse(numerator):
        denom_safe = np.copy(denominator)
        denom_safe[np.isclose(denom_safe, 0)] = np.inf
        result = numerator.multiply(1.0 / denom_safe)
        result.eliminate_zeros()
        return result
    else:
        with np.errstate(divide='ignore', invalid='ignore'):
            result = np.divide(numerator, denominator)
            result = np.where(np.isclose(denominator, 0), 0, result)
        return result


def generate_entropy_metrics(adata, partition_label):
    """
    Generate entropy metrics Psi, Psi_blocks (dataframe) and Zeta.

    Entropy metrics generated:
        - Psi : Fraction of infromation explained by partition of choice
        - Psi_block : Specificity of infromation to a block
        - Zeta : Speicifcty to a partition/ distance of Psi_blocks distribution from uniform

    Parameters
    ----------
    adata : AnnData
        Annotated data object with `.obs` containing metadata.
    partition_label : str
        Column in `.obs` to partition by when calculating entropy metrics.


    Returns
    -------
    Psi : np.Array
        A list of Psi scores (between 0 and 1) corresponding to the selected partition for all genes in `.var`.
    Psi_block_df : pd.Dataframe
        A dataframe of Psi_block scores (between 0 and 1) corresponding to the selected partition for all genes in `.var`.
        Scores are caluclated for all blocks, each column of the dataframe corresponds to one block.
    Zeta : np.Array
        A list of Zeta scores (between 0 and 1) corresponding to the selected partition for all genes in `.var`.

    """

    counts = adata.X
    if sparse.issparse(counts):
        if not np.issubdtype(counts.dtype, np.number):
            counts = counts.astype(np.float32)
    else:
        counts = csr_matrix(np.asarray(counts, dtype=np.float32))

    # Convert to COO once — all entropy work runs on the three raw arrays
    # (row, col, data) without creating any additional sparse matrices.
    coo = counts.tocoo()
    nz_row = np.asarray(coo.row)
    nz_col = np.asarray(coo.col)
    nz_data = coo.data.astype(np.float64)
    n_genes = counts.shape[1]

    total_counts_per_gene = np.asarray(counts.sum(axis=0)).ravel()

    # ── E_T: total entropy per gene ──────────────────────────────────────────
    # p_{ij} = counts_{ij} / total_j;  contribution = -p * log2(p)
    p = nz_data / np.where(total_counts_per_gene[nz_col] > 0,
                            total_counts_per_gene[nz_col], 1.0)
    log_p = np.log2(np.where(p > 0, p, 1.0))
    log_p[p <= 0] = 0
    E_T = np.bincount(nz_col, weights=-p * log_p, minlength=n_genes)
    E_T[np.isclose(total_counts_per_gene, 0)] = -1

    # ── Within-entropy and Psi_block — single vectorised pass ────────────────
    # Tag every non-zero with an integer block ID, then use a compound
    # (block * n_genes + gene) index so np.bincount handles all blocks at once.
    blocks = adata.obs[partition_label].unique()
    n_blocks = len(blocks)
    obs_partition = adata.obs[partition_label].values

    b2id = {b: i for i, b in enumerate(blocks)}
    cell_bid = np.array([b2id[b] for b in obs_partition])
    nz_bid = cell_bid[nz_row]
    bg_idx = nz_bid * n_genes + nz_col   # linearised (block, gene) index

    # Block-gene totals
    bg_sum = np.bincount(bg_idx, weights=nz_data, minlength=n_blocks * n_genes)

    # Within-block probabilities and their entropy contributions
    denom = bg_sum[bg_idx]
    pw = np.where(denom > 0, nz_data / denom, 0.0)
    log_pw = np.log2(np.where(pw > 0, pw, 1.0))
    log_pw[pw <= 0] = 0
    bg_ent = np.bincount(bg_idx, weights=-pw * log_pw, minlength=n_blocks * n_genes)

    bg_sum_2d = bg_sum.reshape(n_blocks, n_genes)   # (n_blocks, n_genes)
    bg_ent_2d = bg_ent.reshape(n_blocks, n_genes)

    # p_c_j = fraction of gene j's total counts belonging to block b
    p_c_j = np.divide(bg_sum_2d, total_counts_per_gene,
                      out=np.zeros((n_blocks, n_genes)),
                      where=~np.isclose(total_counts_per_gene, 0))

    weighted = bg_ent_2d * p_c_j          # (n_blocks, n_genes)
    E_W = weighted.sum(axis=0)             # (n_genes,)
    Psi_block_num = weighted.T             # (n_genes, n_blocks)

    # ── Psi score ─────────────────────────────────────────────────────────────
    with np.errstate(invalid='ignore', divide='ignore'):
        Psi = np.where(E_T > 0, E_W / E_T, -1)

    # ── Psi_block scores ──────────────────────────────────────────────────────
    # Avoid np.divide(..., where=) without out= — its output shape is
    # implementation-defined when the mask is all-False in some numpy builds.
    E_W_safe = np.where(E_W > 0, E_W, np.inf)   # div by inf → 0 for dead genes
    Psi_block = Psi_block_num / E_W_safe[:, None]
    Psi_block[~np.isfinite(Psi_block)] = 0

    Psi_block_df = pd.DataFrame(Psi_block, index=adata.var.index, columns=blocks)

    # ── Zeta score — pure numpy, no DataFrame operations ─────────────────────
    safe_pb = np.where(Psi_block > 0, Psi_block, 1.0)
    log_pb = np.log2(safe_pb)
    log_pb[Psi_block <= 0] = 0
    Zeta = 1 - (-np.sum(Psi_block * log_pb, axis=1) / np.log2(n_blocks))

    return Psi, Psi_block_df, Zeta
