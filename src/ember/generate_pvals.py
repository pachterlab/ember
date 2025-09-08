import pandas as pd
import numpy as np
import os
import gc
import multiprocessing
import anndata as ad
from tqdm import tqdm 
import scanpy as sc
from statsmodels.stats.multitest import multipletests
from .generate_entropy_metrics import generate_entropy_metrics
from .sample_replicates import generate_balanced_draws

def worker_wrapper(args):
    """Unpacks arguments and calls the real worker function."""
    return run_permutation_task(*args)

# 3. THE PRODUCER: Prepares all data and arguments for each task
def prepare_permutation_tasks(h5ad_dir, sets, sample_id_col, partition_label, block_label, n_iterations):
    """
    Generator that handles all file I/O. It prepares and yields a tuple
    of all arguments needed for each permutation task.
    """
    adata = ad.read_h5ad(h5ad_dir, backed='r')
    
    for i in range(n_iterations):
        # Use the provided sets to define the data for this iteration
        # This assumes len(sets) is sufficient for n_iterations
        current_set_index = i % len(sets)
        subset_ids_index = adata.obs[adata.obs[sample_id_col].isin(sets[current_set_index])].index
        subset_in_memory = adata[subset_ids_index, :].to_memory().copy()
        
        # Extract the raw, decoupled data components
        data_matrix = subset_in_memory.X
        obs_df = subset_in_memory.obs
        var_index = subset_in_memory.var.index
        original_labels = subset_in_memory.obs[partition_label].values

        # Yield a single tuple containing all arguments for the worker
        yield (i, data_matrix, obs_df, var_index, original_labels, partition_label, block_label)
        
    adata.file.close()

def run_permutation_task(i, data_matrix, obs_df, var_index, original_labels, partition_label, block_label):
    """
    Lean worker function. Receives raw data, rebuilds an AnnData object,
    runs the permutation, and calculates metrics.
    """
    # Rebuild a clean AnnData object
    subset = ad.AnnData(X=data_matrix, obs=obs_df)
    subset.var.index = var_index
    
    # Perform the core logic: shuffling the labels
    subset.obs[partition_label] = np.random.permutation(original_labels)
    
    # Generate metrics
    Psi, Psi_block_df, Zeta = generate_entropy_metrics(subset, partition_label)
    
    # Format results into Series
    Psi_series = pd.Series(Psi, index=Psi_block_df.index)
    Zeta_series = pd.Series(Zeta, index=Psi_block_df.index)
    
    psi_block_series = None
    if block_label is not None:
        psi_block_series = pd.Series(Psi_block_df[block_label], index=Psi_block_df.index)

    del subset
    gc.collect()
    
    return Psi_series, psi_block_series, Zeta_series


def generate_pvals(
    adata = None,
    h5ad_dir = None,
    partition_label = None,
    entropy_metrics_dir = None,
    is_sampled = True,
    save_dir = None,
    Psi_real = None,
    Psi_block_df_real = None,
    Zeta_real = None,
    sample_id_col = None, 
    category_col = None, 
    condition_col = None,
    seed = 42,
    block_label=None,
    n_iterations=1000,
    n_cpus=1
):
    """
    Calculate empirical p-values for entropy metrics
    from permutation test results.
    
    Entropy metrics generated:
        - Psi : Fraction of infromation explained by partition of choice
        - Psi_block : Specificity of infromation to a block
        - Zeta : Speicifcty to a partition/ distance of Psi_blocks distribution from uniform

    Parameters
    ----------
    adata : AnnData, optional
        An AnnData object already loaded in memory.
        Required if `h5ad_dir` is not provided.
        
    h5ad_dir : str, optional
        Path to the AnnData `.h5ad` file to process.
        Required if `adata` is not provided.

    partition_label : str
        Column in `.obs` used to partition cells for entropy calculations 
        (e.g., "celltype"). Required to run process.
        
    entropy_metrics_dir : str, optional
        Path to csv with entropy metrics to use for .
        Required if `adata` is not provided.
        
    is_sampled : bool, optional
        Boolean to tell the program whether the input entropy metrics are sampled metrics or unsampled metrics. 
        Default is True.
    
    save_dir : str
        Path to directory where results will be saved. Required to run process.
        
    partition_label : str
        Column in `.obs` to partition by when calculating entropy metrics.
        
    Psi_real : pd.Series
        Observed Psi values for each gene.
        
    Psi_block_df_real : pd.Dataframe
        Observed Psi_block values for all blocks in chosen partition. 
        
    Zeta_real : pd.Series
        Observed Zeta values for each gene.
        
    sample_id_col : str
        The column in `.obs` with unique identifiers for each sample or replicate
        (e.g., 'sample_id', 'mouse_id'). Required to run process.
        
    category_col : str
        The column in `.obs` defining the primary groups to balance within
        (e.g., 'disease_status', 'mouse_strain'). Required to run process.
        
    condition_col : str
        The column in `.obs` containing the conditions to balance across within
        each category (e.g., 'sex', 'treatment'). Required to run process.
        
    num_draws : int, optional
        The number of balanced subsets to generate, by default 100.
        
    seed : int, optional
        The random seed for reproducible draws, by default 42.
        
    block_label : str
        Block in partition that the user wishes to calucate p-values for. 
        Default set to None, program will continue generating p-values for only Psi and Zeta. 
        
    n_iterations : int
        Number of iterations to calulate p-vales. Default set to 1000. 
        Note that doing fewer than 1000 iterations is a good choice to get first pass p-values
        but for reliable p-values 1000 iterations is recommended. 
        Larger than 1000 will give you more relibale p-values but will increase runtime significantly. 
        
    n_cpus : int
        Number of cpus to use to perfrom p-value calculation. 
        Default set to 1 assuming no parallel compute power on local machine. 
        User can input -1 to use all available cpus but one. 


    Returns
    -------
    final : pd.DataFrame
        If block is set equal to the default 'none', program outputs dataFrame with columns:
        ['Psi', 'Psi p-value', 'Zeta', 'Zeta p-value']
        
        If block is not equal to 'none', program outputs dataFrame with columns:
        ['Psi', 'Psi p-value', 'Zeta', 'Zeta p-value', Psi_block, 'Psi_block p-value]
        
    """
    # Required argument validation 
    #Validate AnnData

    # 1. Handle the case where the user provides an in-memory object
    if adata is not None:
        if not isinstance(adata, sc.AnnData):
            raise TypeError("The provided `adata` must be an AnnData object.")
        print('Using the provided in-memory AnnData object. Please use `h5ad_dir` argument for backed mode.')

    # 2. Handle the case where the user provides a file path
    elif h5ad_dir is not None:
        adata_path = os.path.expanduser(h5ad_dir)
        if not os.path.exists(adata_path):
            raise FileNotFoundError(f"The file specified by `h5ad_dir` was not found at: {adata_path}")
        
        print(f'Loading AnnData object from {adata_path} in backed mode.')
        adata = sc.read_h5ad(h5ad_dir, backed = 'r')

    # 3. Handle the case where no input is provided
    else:
        raise ValueError("You must provide either an `adata` object or an `h5ad_dir` path.")

     # Validate partition_label
    if partition_label is None:
        raise ValueError("`partition_label` must be provided.")
   
    if partition_label not in adata.obs.columns:
        raise ValueError(
            f"partition_label '{partition_label}' not found in adata.obs columns. "
            f"Available columns: {list(adata.obs.columns)}"
        )
    # Valid block label 
    if block_label is not None:    
        valid_blocks = adata.obs[partition_label].unique()
        if block_label not in valid_blocks:
            raise ValueError(
                f"block_label '{block_label}' not found in adata.obs['{partition_label}']. "
                f"Available block labels: {list(valid_blocks)}"
        )
            
    if entropy_metrics_dir is not None:
        entropy_metrics_dir = os.path.expanduser(entropy_metrics_dir)
        
        if is_sampled:
            all_entropy_metrics = pd.read_csv(os.path.join(entropy_metrics_dir, 
                                                       f"combined_entropy_metrics_{partition_label}.csv"), 
                                          index_col = 0)
            Psi_real = all_entropy_metrics[f"Psi_mean_{partition_label}"]
            Zeta_real = all_entropy_metrics[f"Zeta_mean_{partition_label}"]
            mean_Psi_block_dir = os.path.join(entropy_metrics_dir, "mean_Psi_block_df")
            Psi_block_df_real = pd.read_csv(os.path.join(mean_Psi_block_dir, 
                                                       f"mean_Psi_block_df_{partition_label}.csv"), 
                                          index_col = 0)
        if not is_sampled:
            all_entropy_metrics = pd.read_csv(os.path.join(entropy_metrics_dir, 
                                                       f"entropy_metrics_unsampled_{partition_label}.csv"), 
                                          index_col = 0)
            Psi_real = all_entropy_metrics[f"Psi_unsampled_{partition_label}"]
            Zeta_real = all_entropy_metrics[f"Zeta_unsampled_{partition_label}"]
            unsampled_Psi_block_dir = os.path.join(entropy_metrics_dir, "unsampled_Psi_block_df")
            Psi_block_df_real = pd.read_csv(os.path.join(unsampled_Psi_block_dir, 
                                                       f"Psi_block_unsampled_{partition_label}.csv"), 
                                          index_col = 0)
            
    if entropy_metrics_dir is None and all(
        x is None for x in (Psi_real, Psi_block_df_real, Zeta_real)
        ):
        raise ValueError(
            "You must provide either `entropy_metrics_dir` (path to folder with metrics file) "
            "OR all of `Psi_real`, `Psi_block_df_real`, and `Zeta_real` individually."
        )
    # Validate sampling inputs
    if sample_id_col is None or category_col is None or condition_col is None:
        raise ValueError("You must provide `sample_id_col`, `category_col`, and `condition_col` for sampling of replicates for p-values.")

    
    mask = Psi_real > 0
    Psi_real = Psi_real[mask]
    Psi_block_df_real = Psi_block_df_real[mask]
    Zeta_real = Zeta_real[mask]
    
    # Generate 1000 subsets of replicates to caluclate p-values
    sets, usage = generate_balanced_draws(adata, sample_id_col, category_col, condition_col,num_draws = n_iterations, seed = seed)
    
    #Overwrite adata object from memory if using h5ad_dir to activate backed mode
    if h5ad_dir is not None:
        del(adata)
        gc.collect()
        adata = None
    
    # Function to scramble partition labels to calculate p-values. 
    #Outputs caluclated Psi, Psi_block and Zeta value for each iteration i. 
    
    # Create the generator that will produce the arguments for each task
    task_generator = prepare_permutation_tasks(
        h5ad_dir,
        sets,
        sample_id_col,
        partition_label,
        block_label,
        n_iterations
    )

    print(f"\nStarting {n_iterations} parallel permutations with {n_cpus} workers...")
    results = []
    # Create a pool of worker processes
    with multiprocessing.Pool(processes=n_cpus) as pool:
        # tqdm creates a progress bar for all your tasks
        # imap_unordered is highly efficient
        for result in tqdm(pool.imap_unordered(worker_wrapper, task_generator), total=n_iterations):
            results.append(result)

    Psi_df = pd.DataFrame({f"Psi_{i}": r[0] for i, r in enumerate(results)})
    Zeta_df = pd.DataFrame({f"Zeta_{i}": r[2] for i, r in enumerate(results)})

    Psi_p_values = (Psi_df.ge(Psi_real, axis=0).sum(axis=1) + 1) / (Psi_df.shape[1] + 1)
    Zeta_p_values = (Zeta_df.ge(Zeta_real, axis=0).sum(axis=1) + 1) / (Zeta_df.shape[1] + 1)

    # Define a common index explicitly
    common_index = Psi_block_df_real.index
    
    if block_label is not None:
        psi_block_real = Psi_block_df_real[block_label]
        psi_block_df = pd.DataFrame({f"Psi_block_{i}": r[1] for i, r in enumerate(results)})
        psi_block_p_values = (psi_block_df.ge(psi_block_real, axis=0).sum(axis=1) + 1) / (psi_block_df.shape[1] + 1)
        
        final = pd.DataFrame({
            'Psi': pd.Series(Psi_real, index=common_index),
            'Psi p-value': pd.Series(Psi_p_values, index=common_index), 
            'psi_block': pd.Series(psi_block_real, index=common_index),
            'psi_block p-value': pd.Series(psi_block_p_values, index=common_index), 
            'Zeta': pd.Series(Zeta_real, index=common_index),
            'Zeta p-value': pd.Series(Zeta_p_values, index=common_index), 
        }, index=common_index)

        final.index.name = 'gene_name'
        pval_cols = ['Psi p-value', 'psi_block p-value', 'Zeta p-value']
        
    else:
        final = pd.DataFrame({
            'Psi': pd.Series(Psi_real, index=common_index),
            'Psi p-value': pd.Series(Psi_p_values, index=common_index), 
            'Zeta': pd.Series(Zeta_real, index=common_index),
            'Zeta p-value': pd.Series(Zeta_p_values, index=common_index), 
        }, index=common_index)

        final.index.name = 'gene_name'
        pval_cols = ['Psi p-value', 'Zeta p-value']
        
    # Perfrom global multiple testing correction for all tests and save in FDR columns
   
    # 1) Identify p-value columns (or hardcode your list)
    pval_cols = [c for c in final.columns if c.lower().endswith('p-value')]

    # 2) Collect only valid (non-null, finite) p-values, remembering their (row index, col)
    records = []  # list of (row_index, col_name, pval)
    for col in pval_cols:
        s = final[col]
        mask_valid = s.notna() & np.isfinite(s.values)
        for idx in final.index[mask_valid]:
            records.append((idx, col, float(final.at[idx, col])))

    # 3) If nothing to correct, optionally skip
    if records:
        all_pvals = np.array([r[2] for r in records], dtype=float)

        # 4) Global FDR across ALL p-values (both columns) at once
        _, qvals, _, _ = multipletests(all_pvals, method='fdr_bh')

        # 5) Ensure destination FDR columns exist, initialized to NaN
        for col in pval_cols:
            fdr_col = col.replace('p-value', 'FDR')
            if fdr_col not in final.columns:
                final[fdr_col] = np.nan

        # 6) Write back corrected q-values to the matching rows/columns
        for (idx, col, _), q in zip(records, qvals):
            fdr_col = col.replace('p-value', 'FDR')
            final.at[idx, fdr_col] = q
    else:
        print("No valid p-values found to correct.")

    # 7) Save or return
    if save_dir is not None:
        save_dir = os.path.expanduser(save_dir)
        os.makedirs(save_dir, exist_ok=True)
        if block_label is not None:
            out_path = os.path.join(save_dir, f"pvals_entropy_metrics_{partition_label}_{block_label}.csv")
        else:
            out_path = os.path.join(save_dir, f"pvals_entropy_metrics_{partition_label}.csv")
        final.to_csv(out_path)
        print(f'\nSaved all entropy metrics along with pvalues to {out_path}')
    else:
        return final
    