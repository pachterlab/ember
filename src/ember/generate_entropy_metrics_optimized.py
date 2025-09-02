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

def generate_entropy_metrics_optimized(adata, partition_label):
    """
    A memory-optimized version of generate_entropy_metrics that avoids creating
    large intermediate sparse matrix copies.
    """
    counts = adata.X
    if not sparse.isspmatrix_csr(counts):
        counts = counts.tocsr()

    total_counts_per_gene = np.asarray(counts.sum(axis=0)).ravel()

    # --- Optimized E_T Calculation ---
    # Perform calculations in-place on a single copy to find total entropy
    if 'E_T' not in adata.var.columns:
        # 1. Create a single copy to hold p_i * log2(p_i)
        entropy_matrix = counts.copy().astype(np.float64, copy=False)

        # 2. Directly calculate probabilities on the non-zero data
        #    This is the key step: it avoids creating a new matrix for p_i
        entropy_matrix.data /= total_counts_per_gene[entropy_matrix.indices]
        
        # 3. Now calculate the entropy term, reusing the same matrix
        #    This avoids creating a new matrix for log_p_i
        entropy_matrix.data *= np.log2(entropy_matrix.data)
        
        # 4. Sum to get the final entropy per gene
        E_T = -np.asarray(entropy_matrix.sum(axis=0)).ravel()
        E_T[np.isclose(total_counts_per_gene, 0)] = -1
        adata.var['E_T'] = E_T
    else:
        E_T = adata.var['E_T'].values
    
    # --- Optimized Main Loop ---
    blocks = adata.obs[partition_label].unique()
    n_genes = adata.shape[1]
    n_blocks = len(blocks)
    
    E_W = np.zeros(n_genes, dtype=np.float64)
    Psi_block_num = np.zeros((n_genes, n_blocks), dtype=np.float64)

    for idx, block in enumerate(blocks):
        mask = (adata.obs[partition_label] == block).values
        block_counts = counts[mask, :]
        
        block_sum = np.asarray(block_counts.sum(axis=0)).ravel()
        
        # Skip blocks with no counts to avoid division by zero
        if np.all(np.isclose(block_sum, 0)):
            continue
            
        # Use the same in-place optimization for the within-block entropy
        block_entropy_matrix = block_counts.copy().astype(np.float64, copy=False)
        
        # Use a mask to handle genes with zero sum in this block
        valid_genes_mask = ~np.isclose(block_sum[block_entropy_matrix.indices], 0)
        
        # Filter data and indices to only valid genes
        original_data = block_entropy_matrix.data[valid_genes_mask]
        original_indices = block_entropy_matrix.indices[valid_genes_mask]

        # Calculate q_j and entropy term in one go
        q_j = original_data / block_sum[original_indices]
        entropy_data = -q_j * np.log2(q_j)
        
        # Reconstruct the sparse entropy matrix for summing
        entropy_matrix_block = csr_matrix((entropy_data, original_indices, block_entropy_matrix.indptr), shape=block_entropy_matrix.shape)
        entropy = np.asarray(entropy_matrix_block.sum(axis=0)).ravel()
        
        p_c_j = np.divide(
            block_sum,
            total_counts_per_gene,
            out=np.zeros_like(block_sum, dtype=np.float64),
            where=~np.isclose(total_counts_per_gene, 0)
        )
        
        weighted_entropy = entropy * p_c_j
        E_W += weighted_entropy
        Psi_block_num[:, idx] = weighted_entropy

    # --- Final calculations (unchanged) ---
    with np.errstate(invalid='ignore', divide='ignore'):
        Psi = np.where(E_T > 0, E_W / E_T, -1)

    with np.errstate(divide='ignore', invalid='ignore'):
        Psi_block = np.divide(Psi_block_num, E_W[:, None], where=E_W[:, None] != 0)
        Psi_block[~np.isfinite(Psi_block)] = 0

    Psi_block_df = pd.DataFrame(Psi_block, index=adata.var.index, columns=blocks)
    
    entropy = -np.nansum(Psi_block_df * np.log2(Psi_block_df.where(Psi_block_df > 0)), axis=1)
    Zeta = 1 - (entropy / np.log2(len(Psi_block_df.columns)))
    
    return Psi, Psi_block_df, Zeta