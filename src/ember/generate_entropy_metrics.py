import os
import pandas as pd
import numpy as np
import random
from collections import defaultdict
from itertools import product
import scanpy as sc
from scipy import sparse
import warnings
from scipy.sparse import csr_matrix, lil_matrix


def safe_log2_sparse(mat):
    """
    
    Helper function applies an element-wise log2 transform only to the non-zero
    *stored* entries in a SciPy CSR sparse matrix. Zeros in a sparse matrix are
    implicit (not stored in `.data`) and remain zero in the output.

    Any non-finite results from the log transform (e.g., `-inf` from explicit
    stored zeros, or `NaN` from negative values) are replaced with 0.
    Parameters
    ----------
    mat : scipy.sparse.csr_matrix
        Input sparse matrix to be log2-transformed.

    Returns
    -------
    scipy.sparse.csr_matrix
        Sparse matrix of the same shape with log2-transformed data.
        Any NaN or infinite values are replaced with 0.
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

    Parameters
    ----------
    numerator : array-like or scipy.sparse.spmatrix
        Matrix or sparse matrix to be divided.
    denominator : array-like
        Vector or array to divide by. Must be broadcast-compatible with numerator.

    Returns
    -------
    array-like or scipy.sparse.spmatrix
        Result of the element-wise division with safe handling of zero denominators.
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
    # Input sparse matrix
    if not sparse.issparse(counts):
        counts = csr_matrix(counts)

    # Row-wise sum
    total_counts_per_gene = np.asarray(counts.sum(axis=0)).ravel()

    # Skip if total entropy has already been calculated 
    if 'E_T' not in adata.var.columns:
        # p_i = P(gene | total) across all cells
        p_i = safe_divide_sparse(counts, total_counts_per_gene)
        log_p_i = safe_log2_sparse(p_i)
        
        # Total entropy (E_T)
        entropy_per_gene = p_i.multiply(log_p_i).sum(axis=0).A1
        E_T = -entropy_per_gene
        E_T[np.isclose(total_counts_per_gene, 0)] = -1
        adata.var['E_T'] = E_T
    else:
        E_T = adata.var['E_T'].values

    blocks = adata.obs[partition_label].unique()
    n_genes = adata.shape[1]
    n_blocks = len(blocks)
    
    # Initialize E_W (within entropy) and Psi_block numerator (Psi_block_num)
    E_W = np.zeros(n_genes)
    Psi_block_num = np.zeros((n_genes, n_blocks))

    for idx, block in enumerate(blocks):
        mask = adata.obs[partition_label] == block
        block_counts = adata[mask, :].X
        # Calculate entripy within each block
        block_sum = np.asarray(block_counts.sum(axis=0)).ravel()
        q_j = safe_divide_sparse(block_counts, block_sum)
        log_q_j = safe_log2_sparse(q_j)
        entropy = -q_j.multiply(log_q_j).sum(axis=0).A1

        # Contribution of this block to Within Entropy (E_W)
        p_c_j = np.divide(
            block_sum,
            total_counts_per_gene,
            out=np.zeros_like(block_sum),
            where=~np.isclose(total_counts_per_gene, 0)
        )

        weighted_entropy = entropy * p_c_j
        E_W += weighted_entropy
        Psi_block_num[:, idx] = weighted_entropy

    # Psi score
    with np.errstate(invalid='ignore', divide='ignore'):
        Psi = np.where(E_T > 0, E_W / E_T, -1)

    # Psi_block scores
    with np.errstate(divide='ignore', invalid='ignore'):
        Psi_block = np.divide(Psi_block_num, E_W[:, None], where=E_W[:, None] != 0)
        Psi_block[~np.isfinite(Psi_block)] = 0

    Psi_block_df = pd.DataFrame(Psi_block, index=adata.var.index, columns=blocks)
    
    # Zeta score 
    entropy = -np.nansum(Psi_block_df * np.log2(Psi_block_df.where(Psi_block_df > 0)), axis=1)
    Zeta = 1 - (entropy / np.log2(len(Psi_block_df.columns)))
    

    return Psi, Psi_block_df, Zeta


def generate_entropy_metrics_optimized(adata, partition_label):
    """
    Memory efficient version!!
    
    Generate entropy metrics Psi, Psi_blocks (dataframe) and Zeta.
    
    Entropy metrics generated:
        - Psi : Fraction of infromation explained by partition of choice
        - Psi_block : Specificity of infromation to a block
        - Zeta : Speicifcty to a partition/ distance of Psi_blocks distribution from uniform
        
    This version minimizes memory usage by avoiding large intermediate sparse matrices
    and directly integrates the logic from the safe helper functions.
    
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
    if not sparse.issparse(counts):
        counts = csr_matrix(counts)

    total_counts_per_gene = np.asarray(counts.sum(axis=0)).ravel()

    if 'E_T' not in adata.var.columns:
        # Create a copy to modify its .data attribute directly
        p_i = counts.copy().astype(np.float64)
        
        # Apply the logic from safe_divide_sparse
        p_i.data /= np.take(total_counts_per_gene, p_i.indices)
        
        # Apply the logic from safe_log2_sparse
        with np.errstate(divide='ignore'):
            logged_data = np.log2(p_i.data)
        logged_data[~np.isfinite(logged_data)] = 0
        p_i.data *= logged_data
        
        E_T = -np.asarray(p_i.sum(axis=0)).ravel()
        E_T[np.isclose(total_counts_per_gene, -1)] = 0 
        adata.var['E_T'] = E_T
    else:
        E_T = adata.var['E_T'].values

    blocks = adata.obs[partition_label].unique()
    n_genes = adata.shape[1]
    
    E_W = np.zeros(n_genes)
    Psi_block_num = np.zeros((n_genes, len(blocks)))

    for idx, block in enumerate(blocks):
        mask = adata.obs[partition_label] == block
        block_counts = counts[mask, :]
        block_sum = np.asarray(block_counts.sum(axis=0)).ravel()

        if np.all(np.isclose(block_sum, 0)):
            continue

        q_j = block_counts.copy().astype(np.float64)
        q_j.data /= np.take(block_sum, q_j.indices)
        
        with np.errstate(divide='ignore'):
            logged_data = np.log2(q_j.data)
        logged_data[~np.isfinite(logged_data)] = 0
        q_j.data *= logged_data
        
        entropy = -np.asarray(q_j.sum(axis=0)).ravel()
        entropy[np.isclose(block_sum, 0)] = 0
        
        p_c_j = np.divide(block_sum, total_counts_per_gene, 
                          out=np.zeros_like(block_sum, dtype=float), 
                          where=total_counts_per_gene > 0)

        weighted_entropy = entropy * p_c_j
        E_W += weighted_entropy
        Psi_block_num[:, idx] = weighted_entropy

    with np.errstate(invalid='ignore', divide='ignore'):
        Psi = np.divide(E_W, E_T, out=np.zeros_like(E_T, dtype=float), where=E_T > 0)

    with np.errstate(invalid='ignore', divide='ignore'):
        E_W_expanded = E_W[:, np.newaxis]
        Psi_block = np.divide(Psi_block_num, E_W_expanded, 
                              out=np.zeros_like(Psi_block_num, dtype=float), 
                              where=E_W_expanded != 0)

    Psi_block_df = pd.DataFrame(Psi_block, index=adata.var.index, columns=blocks)
    
    with np.errstate(divide='ignore'):
        temp_psi_block = Psi_block_df.where(Psi_block_df > 0)
        log_psi_block = np.log2(temp_psi_block)
    
    block_entropy = -np.nansum(Psi_block_df * log_psi_block, axis=1)
    
    Zeta = 1 - (block_entropy / np.log2(len(blocks)))

    return Psi, Psi_block_df, Zeta
