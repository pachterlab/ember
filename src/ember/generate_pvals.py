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

def _worker_wrapper(args):
    """
    Unpacks a tuple of arguments and calls the real worker function.
    Needed because pool.imap only accepts functions with a single argument.
    """
    return _run_permutation_task(*args)


def _prepare_permutation_tasks(h5ad_dir, sets, sample_id_col, partition_label, block_label, n_iterations):
    """
    Generator that handles all file I/O. It prepares and yields a tuple
    of all arguments needed for each permutation task.
    """
    adata = ad.read_h5ad(h5ad_dir, backed='r')
    
    for i in range(n_iterations):
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

def _run_permutation_task(i, data_matrix, obs_df, var_index, original_labels, partition_label, block_label):
    """
    Worker that accepts raw data and rebuilds a clean AnnData object.
    Runs generate_entropy_metrics on rebuild anndata.
    """
    # Rebuild a clean AnnData object
    subset = ad.AnnData(X=data_matrix, obs=obs_df)
    subset.var.index = var_index
    
    # Shuffling the labels
    subset.obs[partition_label] = np.random.permutation(original_labels)
    
    # Generate metrics
    Psi, Psi_block_df, Zeta = generate_entropy_metrics(subset, partition_label)
    
    Psi_series = pd.Series(Psi, index=Psi_block_df.index)
    Zeta_series = pd.Series(Zeta, index=Psi_block_df.index)
    
    psi_block_series = None
    if block_label is not None:
        psi_block_series = pd.Series(Psi_block_df[block_label], index=Psi_block_df.index)

    del subset
    gc.collect()
    
    return Psi_series, psi_block_series, Zeta_series


def generate_pvals(
    h5ad_dir ,
    partition_label,
    entropy_metrics_dir,
    save_dir,
    sample_id_col, 
    category_col, 
    condition_col,
    block_label=None,
    seed = 42,
    n_iterations=1000,
    n_cpus=1, 
    Psi_real = None,
    Psi_block_df_real = None,
    Zeta_real = None
):
    """
    Calculate empirical p-values for entropy metrics from permutation test results.
    This function can be called manually or accessed through light_ember with 
    partition_pvals = True or block_pvals = True. 

    Manual access useful if wanting to generate p-values for multiple blocks and partitions of 
    interest after initial investigation using entropy metrics. 
    
    Integrated access with light_ember is easier if investigating only a partition or
    a block in a partition. 
    
    Entropy metrics generated:
        - Psi : Fraction of infromation explained by partition of choice
        - Psi_block : Specificity of infromation to a block
        - Zeta : Speicifcty to a partition/ distance of Psi_blocks distribution from uniform

    Parameters
    ----------
    h5ad_dir : str, Required
        Path to the `.h5ad` file to process.
        Data should be log1p and depth normalized before running ember. 
        Remove genes with less than 100 reads before running ember. 
        
    partition_label : str, Required
        Column in `.obs` used to partition cells for entropy calculations 
        (e.g., "celltype", "Genotype", "Age"). Required to run process. 
        If performing calculation on interaction term, first create a column 
        in `.obs` concatnating the two columns of interested with a semicolon (:).
        
    entropy_metrics_dir : str, Required
        Path to csv with entropy metrics to use for generating pvals.
    
    save_dir : str, Required
        Path to directory where results will be saved. 
        
    sample_id_col : str, Required
        The column in `.obs` with unique identifiers for each sample or replicate
        (e.g., 'sample_id', 'mouse_id').
        
    category_col : str, Required
        The column in `.obs` defining the primary group to balance across in order
        to generate a balanced sample of the experiment. (e.g., 'disease_status', 'mouse_strain').
        Refer to readme for further explanation on how to select category and condition columns.
        category_col and condition_col are interchangable.
        If balancing across more than 2 variables, generate interaction terms, create a column
        in `.obs` concatnating the two (or more) columns of interested with a semicolon (:). 
        This way balancing can be done across as many variables as desired. 
        
    condition_col : str, Required
        The column in `.obs` containing the conditions to balance within
        each category to generate a balanced sample of the experiment.  (e.g., 'sex', 'treatment').
        Refer to readme for further explanation on how to select category and condition columns. 
        category_col and condition_col are interchangable. 
        If balancing across more than 2 variables, generate interaction terms, create a column
        in `.obs` concatnating the two (or more) columns of interested with a semicolon (:). 
        This way balancing can be done across as many variables as desired. 
        
    block_label : str, default=None
        Block in partition to calucate p-values for. 
        Default set to None, program will continue generating p-values for only Psi and Zeta. 
    
    seed : int, default=42
        The random seed for reproducible draws, by default 42.

    n_iterations : int, default = 1000
        Number of iterations to calulate p-vales. Default set to 1000. 
        Note that doing fewer than 1000 iterations is a good choice to get first pass p-values
        but for reliable p-values 1000 iterations is recommended. 
        Larger than 1000 will give you more relibale p-values but will increase runtime significantly. 
        
    n_cpus : int, default=1
        Number of cpus to use to perfrom p-value calculation. 
        Default set to 1 assuming no parallel compute power on local machine. 
        User can input -1 to use all available cpus but one.
    
    Psi_real : pd.Series, default=None
        Observed Psi values for each gene. 
        Used by light_ember, not necessary for user use. 
        
    Psi_block_df_real : pd.Dataframe, default = None
        Observed Psi_block values for all blocks in chosen partition. 
        Used by light_ember, not necessary for user use.
        
    Zeta_real : pd.Series, default=None
        Observed Zeta values for each gene.
        Used by light_ember, not necessary for user use.
        

    Returns
    -------
    None
        
    Notes
    -------
    **What to expect inside 'pvals_entropy_metrics.csv'**:
    
    - gene_name: All genes in `.var`
    - Psi: Psi scores averaged over n draws (between 0 and 1) generated by light_ember for each gene in `.var`.
    - Psi p-value: Permutation based empirical p-values for observed Psi scores for each gene in `.var`.
    - Zeta: Zeta scores averaged over n draws (between 0 and 1) generated by light_ember for each gene in `.var`.
    - Zeta p-value: Permutation based empirical p-values for observed Zeta scores for each gene in `.var`.
    - Psi FDR: Multiple testing corrected q-values for Psi scores.  
    - Zeta FDR: Multiple testing corrected q-values for Zeta scores.Correction perfromed to include all p-values generated in a single file (Psi and Zeta). 
    
    if block_pvals = True and a single block_label is given:
    
    - psi_block: psi_block scores (between 0 and 1) generated by light_ember for each gene in `.var`.
    - psi_block p-value: Permutation based empirical p-values for observed psi_block scores for each gene in `.var`.
    - psi_block FDR: Multiple testing corrected q-values for psi_block scores. Correction perfromed to include all p-values generated in a single file (Psi, psi_block and Zeta). 
                     
    """
    #Validate h5ad_dir
    adata_path = os.path.expanduser(h5ad_dir)
    if not os.path.exists(adata_path):
        raise FileNotFoundError(f"The file specified by `h5ad_dir` was not found at: {adata_path}")
        
    print(f'Loading AnnData object from {adata_path} in backed mode.')
    adata = sc.read_h5ad(h5ad_dir, backed = 'r')

    # Validate partition_label
    if partition_label not in adata.obs.columns:
        raise ValueError(
            f"partition_label '{partition_label}' not found in adata.obs columns. "
            f"Available columns: {list(adata.obs.columns)}"
        )
        
    # Validate entropy_metrics_dir
    entropy_path = os.path.expanduser(entropy_metrics_dir)
    if not os.path.exists(entropy_path):
        raise FileNotFoundError(f"The folder specified by `entropy_metrics_dir` was not found at: {entropy_path}")
        
    # Valid block label 
    if block_label is not None:    
        valid_blocks = adata.obs[partition_label].unique()
        if block_label not in valid_blocks:
            raise ValueError(
                f"block_label '{block_label}' not found in adata.obs['{partition_label}']. "
                f"Available block labels: {list(valid_blocks)}"
        )
    # Validate sample_id_col, category_col, and condition_col 
    required_params = {
        "sample_id_col": sample_id_col, 
        "category_col": category_col, 
        "condition_col": condition_col
    }

    required_cols = list(required_params.values())
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]

    if missing_cols:
        raise ValueError(
            f"The following required column(s) were not found in adata.obs: {missing_cols}. "
            f"\nAvailable columns are: {list(adata.obs.columns)}"
        )
        
    #Override entropy_metrics_dir if required arguments are provided. 
    #Functionality for intergration with light_ember
    
    if any(v is not None for v in (Psi_real, Psi_block_df_real, Zeta_real)):
        entropy_metrics_dir = os.path.expanduser(entropy_metrics_dir)
        all_entropy_metrics = pd.read_csv(os.path.join(entropy_metrics_dir, 
                                                   f"entropy_metrics_{partition_label}.csv"), 
                                      index_col = 0)
        Psi_real = all_entropy_metrics[f"Psi_mean_{partition_label}"]
        Zeta_real = all_entropy_metrics[f"Zeta_mean_{partition_label}"]
        mean_Psi_block_dir = os.path.join(entropy_metrics_dir, "Psi_block_df")
        Psi_block_df_real = pd.read_csv(os.path.join(mean_Psi_block_dir, 
                                                   f"mean_Psi_block_df_{partition_label}.csv"), 
                                      index_col = 0)
            
    print(f'\nGenerating p-values for entropy metrics.')
    
    mask = Psi_real > 0
    Psi_real = Psi_real[mask]
    Psi_block_df_real = Psi_block_df_real[mask]
    Zeta_real = Zeta_real[mask]
    
    # Generate 1000 subsets of replicates to caluclate p-values
    sets, usage = generate_balanced_draws(adata, sample_id_col, category_col, condition_col,num_draws = n_iterations, seed = seed)
    del(adata)
    gc.collect()
    
    # Create the generator that will produce the arguments for each task
    task_generator = _prepare_permutation_tasks(
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
        for result in tqdm(pool.imap_unordered(_worker_wrapper, task_generator), total=n_iterations):
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
   
    pval_cols = [c for c in final.columns if c.lower().endswith('p-value')]
    records = [] 
    for col in pval_cols:
        s = final[col]
        mask_valid = s.notna() & np.isfinite(s.values)
        for idx in final.index[mask_valid]:
            records.append((idx, col, float(final.at[idx, col])))

    all_pvals = np.array([r[2] for r in records], dtype=float)

    # Global FDR across ALL p-values (both columns) at once
    _, qvals, _, _ = multipletests(all_pvals, method='fdr_bh')

    for col in pval_cols:
        fdr_col = col.replace('p-value', 'FDR')
        if fdr_col not in final.columns:
            final[fdr_col] = np.nan

    # Write back corrected q-values to the matching rows/columns
    for (idx, col, _), q in zip(records, qvals):
        fdr_col = col.replace('p-value', 'FDR')
        final.at[idx, fdr_col] = q


    save_dir = os.path.expanduser(save_dir)
    os.makedirs(save_dir, exist_ok=True)
    if block_label is not None:
        out_path = os.path.join(save_dir, f"pvals_entropy_metrics_{partition_label}_{block_label}.csv")
    else:
        out_path = os.path.join(save_dir, f"pvals_entropy_metrics_{partition_label}.csv")
    final.to_csv(out_path)
    print(f'\nSaved all entropy metrics along with pvalues to {out_path}')
    