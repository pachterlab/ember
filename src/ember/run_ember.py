import os
import gc
import glob
import shutil
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from tempfile import mkdtemp
from joblib import Parallel, delayed, parallel_backend
from .generate_entropy_metrics import generate_entropy_metrics_optimized
from .generate_pvals import generate_pvals
from .sample_individuals import generate_balanced_sets, aitchison_mean_and_std


def run_ember(
    adata = None,
    adata_dir = None,
    partition_label = None,
    save_dir = None,
    sampling=True,
    individual_var=None,
    sex_var=None,
    balance_by_var=None,
    num_draws=100,
    partition_pvals=True,
    block_pvals=False,
    block_label=None,
    n_pval_iterations=1000,
    n_cpus=1
):
    """
    Runs the ember entropy metrics and p-value generation workflow on an AnnData object.

    This function loads an AnnData `.h5ad` file, optionally performs balanced sampling
    across individuals, computes entropy metrics for the specified partition,
    and generates p-values for Psi and Zeta and optionally Psi_block for a block of choice.
    
    Entropy metrics generated:
        - Psi : Fraction of infromation explained by partition of choice
        - Psi_block : Specificity of infromation to a block
        - Zeta : Speicifcty to a partition/ distance of Psi_blocks distribution from uniform

    Parameters
    ----------
    adata_dir : str, optional
        Path to the AnnData `.h5ad` file to process.
        Required if `adata` is not provided.

    adata : AnnData, optional
        An AnnData object already loaded in memory.
        Required if `adata_dir` is not provided.
        
    partition_label : str
        Column in `.obs` used to partition cells for entropy calculations 
        (e.g., "celltype"). Required to run process. 
        
    save_dir : str
        Path to directory where results will be saved. Required to run process. 

    sampling : bool, default=True
        Whether to perform balanced sampling across individuals before entropy calculation.
        If True, `individual_var`, `sex_var`, and `balance_by_var` must be provided.

    individual_var : str, optional
        Column in `.obs` containing unique individual IDs.
        Required if `sampling=True`.

    sex_var : str, optional
        Column in `.obs` containing sex categories (e.g., "Male", "Female").
        Required if `sampling=True`.

    balance_by_var : str, optional
        Column in `.obs` used to balance individual sampling 
        (e.g., strain or genotype). Required if `sampling=True`.

    num_draws : int, default=100
        Number of balanced draws to perform if `sampling=True`.

    partition_pvals : bool, default=True
        Whether to compute permutation-based p-values for the `partition_label`.

    block_pvals : bool, default=False
        Whether to compute permutation-based p-values for the `block_label`.

    block_label : str, optional
        Column in `.obs` to use for block-based permutation tests.
        Required if `block_pvals=True`.

    n_pval_iterations : int, default=1000
        Number of permutations to use for p-value calculation.

    n_cpus : int, default=1
        Number of CPU cores to use for parallel permutation testing.

    Notes
    -----
    - Results are saved to `save_dir` as CSV files.
    
    """

        

    #Validate AnnData

    # 1. Handle the case where the user provides an in-memory object
    if adata is not None:
        if not isinstance(adata, sc.AnnData):
            raise TypeError("The provided `adata` must be an AnnData object.")
        print('Using the provided in-memory AnnData object. Please use `adata_dir` argument for backed mode.')

    # 2. Handle the case where the user provides a file path
    elif adata_dir is not None:
        adata_path = os.path.expanduser(adata_dir)
        if not os.path.exists(adata_path):
            raise FileNotFoundError(f"The file specified by `adata_dir` was not found at: {adata_path}")
        
        print(f'Loading AnnData object from {adata_path} in backed mode.')
        adata = sc.read_h5ad(adata_dir, backed = 'r')

    # 3. Handle the case where no input is provided
    else:
        raise ValueError("You must provide either an `adata` object or an `adata_dir` path.")

    
     # Validate partition_label
    if partition_label is None:
        raise ValueError("`partition_label` must be provided.")
   
    if partition_label not in adata.obs.columns:
        raise ValueError(
            f"partition_label '{partition_label}' not found in adata.obs columns. "
            f"Available columns: {list(adata.obs.columns)}"
        )
    # Validate save_dir
    if save_dir is None:
        raise ValueError("`save_dir` must be provided.")
    
    
   # Validate required params if sampling or pvals true
    if sampling or partition_pvals or block_pvals:
        if not all([individual_var, sex_var, balance_by_var]):
            raise ValueError(
                "If sampling or generating pvals (requires sampling), you must provide `individual_var`, `sex_var`, and `balance_by_var`."
            )


    # Validate block_pvals requirements
    if block_pvals:
        if not all([block_label]):
            raise ValueError(
                "If `block_pvals=True`, you must provide `block_label`."
            )

        valid_blocks = adata.obs[partition_label].unique()
        if block_label not in valid_blocks:
            raise ValueError(
                f"block_label '{block_label}' not found in adata.obs['{partition_label}']. "
                f"Available block labels: {list(valid_blocks)}"
        )
    
    if sampling:
        # Generate balanced sets
        sets, usage = generate_balanced_sets(
            adata, individual_var, sex_var, balance_by_var, num_draws=num_draws
        )

        # Temporary directory
        temp_dir = mkdtemp(prefix="ember_balanced_")
        print(f"\nSaving intermediate results to temp dir: {temp_dir}")
        
        #Overwrite adata object from memory if using adata_dir to activate backed mode
        if adata_dir is not None:
            del(adata)
            gc.collect()
            adata = None
            
        try:
            

            def run_draw(i,  
                         partition_label, 
                         individual_var, 
                         sets, 
                         temp_dir, 
                         adata = None, 
                         adata_dir = None):
                
                draw_code = f"BALANCED_{i:02d}"
                draw_dir = os.path.join(temp_dir, draw_code)
                os.makedirs(draw_dir, exist_ok=True)
                
                if adata_dir is not None:
                    adata = sc.read_h5ad(adata_dir, backed = 'r')

                # Only load small subsets to memory to maximize speed.
                subset_ids_index = adata.obs[adata.obs[individual_var].isin(sets[i])].index
                subset = adata[subset_ids_index, :].to_memory()
                temp_psi, temp_psi_block, temp_zeta = generate_entropy_metrics_optimized(subset, partition_label)

                # Save entropy metrics
                entropy_df = pd.DataFrame({
                    f"Psi_{partition_label}": temp_psi, 
                    f"Zeta_{partition_label}": temp_zeta
                }, index=subset.var.index)
                entropy_df.to_csv(os.path.join(draw_dir, f"entropy_metrics_{partition_label}.csv"))
                del(subset)
                gc.collect()

                # Save Psi_block
                temp_psi_block.to_csv(os.path.join(draw_dir, f"Psi_block_{partition_label}.csv"))
            
            print(f'\nGenerating entropy metrics for {num_draws} samples. Using {n_cpus} CPUs')
            
            # Run paralellized permutaations if compute power permits. Default set to 1 unless changed in n_cpus. 
            with parallel_backend('loky'):
                results = Parallel(n_jobs=n_cpus, verbose=10)(
                    delayed(run_draw)(
                        i, 
                        partition_label,
                        individual_var,
                        sets,
                        temp_dir, 
                        adata, 
                        adata_dir
                    ) for i in range(num_draws)
                )


            # ========== Aggregation Phase ==========
            # Prep
            save_dir = os.path.expanduser(save_dir)
            os.makedirs(save_dir, exist_ok=True)
            print(f'\nAggregating entropy metrics from {num_draws} samples.')

            all_Psi, all_Zeta, Psi_block_list = [], [], []

            # Load all runs first
            for pair_dir in sorted(glob.glob(os.path.join(temp_dir, "BALANCED_*"))):
                var_df = pd.read_csv(os.path.join(pair_dir, f"entropy_metrics_{partition_label}.csv"), index_col=0)
                sib_df = pd.read_csv(os.path.join(pair_dir, f"Psi_block_{partition_label}.csv"), index_col=0)
                gene_names = var_df.index

                all_Psi.append(var_df[f"Psi_{partition_label}"].values)
                all_Zeta.append(var_df[f"Zeta_{partition_label}"].values)
                Psi_block_list.append(sib_df)

            # Align SIB columns
            common_blocks = sorted(set.intersection(*(set(df.columns) for df in Psi_block_list)))
            Psi_block_list = [df[common_blocks] for df in Psi_block_list]

            # Compute mean/std ONCE across all draws
            Psi_arr = np.stack(all_Psi)
            Zeta_arr = np.stack(all_Zeta)
            #Psi for genes with close to zero counts is set to -1 in generate_entropy_metrics. Masking over these values to avoid miscaulcations. 
            mask = Psi_arr != -1

            # Suppress numpy warnings for empty slices and degrees of freedom <= 0
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)

                aggregate_entropy_df = pd.DataFrame({
                    f"Psi_mean_{partition_label}" : np.nanmean(np.where(mask, Psi_arr, np.nan), axis=0),
                    f"Psi_std_{partition_label}"  : np.nanstd(np.where(mask, Psi_arr, np.nan), axis=0),
                    f"Psi_valid_counts_{partition_label}" : np.sum(mask, axis=0),
                    f"Zeta_mean_{partition_label}": np.nanmean(np.where(mask, Zeta_arr, np.nan), axis=0),
                    f"Zeta_std_{partition_label}" : np.nanstd(np.where(mask, Zeta_arr, np.nan), axis=0),
                    f"Zeta_valid_counts_{partition_label}" : np.sum(mask, axis=0)
                }, index=gene_names)

            # Replace mean/std with NaN when there are no valid counts
            aggregate_entropy_df.loc[
                aggregate_entropy_df[f"Psi_valid_counts_{partition_label}"] == 0,
                [f"Psi_mean_{partition_label}", f"Psi_std_{partition_label}"]
            ] = np.nan

            aggregate_entropy_df.loc[
                aggregate_entropy_df[f"Zeta_valid_counts_{partition_label}"] == 0,
                [f"Zeta_mean_{partition_label}", f"Zeta_std_{partition_label}"]
            ] = np.nan


            # Save Psi_block means/stds
            mean_Psi_block_dir = os.path.join(save_dir, "mean_Psi_block_df")
            os.makedirs(mean_Psi_block_dir, exist_ok=True)
            aitchison_results = aitchison_mean_and_std(Psi_block_list)
            mean_Psi_block = pd.DataFrame(aitchison_results[0], index=gene_names, columns=common_blocks)
            std_Psi_block = pd.DataFrame(aitchison_results[1], index=gene_names, columns=common_blocks)
            mean_Psi_block.to_csv(os.path.join(mean_Psi_block_dir, f"mean_Psi_block_df_{partition_label}.csv"))
            std_Psi_block.to_csv(os.path.join(mean_Psi_block_dir, f"std_Psi_block_df_{partition_label}.csv"))

            # Save final aggregated file
            aggregate_entropy_df.to_csv(os.path.join(save_dir, f"combined_entropy_metrics_{partition_label}.csv"))
            print(f'\nSaved all entropy metrics to {os.path.join(save_dir, f"combined_entropy_metrics_{partition_label}.csv")}')

            if partition_pvals or block_pvals:
                print(f'\nGenerating p-values for entropy metrics.')
            # If partition_pvals == True and block_pvals == False
            if partition_pvals and not block_pvals:
                pvals = generate_pvals(adata = adata, 
                                       partition_label = partition_label, 
                                       Psi_real = aggregate_entropy_df[f"Psi_mean_{partition_label}"], 
                                       Psi_block_df_real = mean_Psi_block, 
                                       Zeta_real = aggregate_entropy_df[f"Zeta_mean_{partition_label}"], 
                                       individual_var = individual_var,
                                       sex_var = sex_var, 
                                       balance_by_var = balance_by_var, 
                                       block_label=block_label, 
                                       n_iterations=n_pval_iterations, 
                                       n_cpus=n_cpus)
                out_path = os.path.join(save_dir, f"pvals_combined_entropy_metrics_{partition_label}.csv")
                
                pvals.to_csv(out_path)
                print(f'\nSaved all entropy metrics along with pvalues to {out_path}')
                

            # Run when block_pvals is True whether or not partition_pvals is  
            elif block_pvals:
                pvals = generate_pvals(adata = adata,
                                       partition_label = partition_label, 
                                       Psi_real = aggregate_entropy_df[f"Psi_mean_{partition_label}"], 
                                       Psi_block_df_real = mean_Psi_block, 
                                       Zeta_real = aggregate_entropy_df[f"Zeta_mean_{partition_label}"], 
                                       individual_var = individual_var, 
                                       sex_var = sex_var, 
                                       balance_by_var = balance_by_var, 
                                       block_label=block_label, 
                                       n_iterations=n_pval_iterations, 
                                       n_cpus=n_cpus)
                out_path = os.path.join(save_dir, f"pvals_combined_entropy_metrics_{partition_label}_{block_label}.csv")
                
                pvals.to_csv(out_path)
                print(f'\nSaved all entropy metrics along with pvalues to {out_path}')

            
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)
            print(f"\nDeleted temp dir: {temp_dir}")



    else:
        # If sampling == False
        print("\nSampling turned off, generating entropy metrics without sampling")
        adata.to_memory()
        temp_psi, temp_psi_block, temp_zeta = generate_entropy_metrics(adata, partition_label)
        # Save entropy metrics
        entropy_df = pd.DataFrame({
            f"Psi_unsampled_{partition_label}": temp_psi, 
            f"Zeta_unsampled_{partition_label}": temp_zeta

        }, index=adata.var.index)
        entropy_df.to_csv(os.path.join(save_dir, f"entropy_metrics_unsampled_{partition_label}.csv"))

        # Save Psi_block
        unsampled_Psi_block_dir = os.path.join(save_dir, "unsampled_Psi_block_df")
        os.makedirs(unsampled_Psi_block_dir, exist_ok=True)
        temp_psi_block.to_csv(os.path.join(unsampled_Psi_block_dir, 
                                           f"Psi_block_unsampled_{partition_label}.csv"))

        print(f'\nSaved all unsampled entropy metrics to {os.path.join(save_dir, f"entropy_metrics_unsampled_{partition_label}.csv")}')
        
        if partition_pvals or block_pvals:
            print(f'\nGenerating p-values for unsampled entropy metrics.')
            
        # If partition_pvals == True and block_pvals == False
        if partition_pvals and not block_pvals:
            pvals = generate_pvals(adata = adata,
                                   partition_label = partition_label, 
                                   Psi_real = entropy_df[f"Psi_unsampled_{partition_label}"], 
                                   Psi_block_df_real = temp_psi_block, 
                                   Zeta_real = entropy_df[f"Zeta_unsampled_{partition_label}"], 
                                   individual_var = individual_var, 
                                   sex_var = sex_var, 
                                   balance_by_var = balance_by_var, 
                                   block_label=block_label, 
                                   n_iterations=n_pval_iterations, 
                                   n_cpus=n_cpus)
            out_path = os.path.join(save_dir, f"pvals_entropy_metrics_unsampled_{partition_label}_{block_label}.csv")
            pvals.to_csv(out_path)
            print(f'\nSaved all unsampled entropy metrics with pvalues to {out_path}')
         
         # Run when block_pvals = True whether or not partition_pvals is   
        elif block_pvals:
            pvals = generate_pvals(adata = adata, 
                                   partition_label = partition_label, 
                                   Psi_real = entropy_df[f"Psi_unsampled_{partition_label}"], 
                                   Psi_block_df_real = temp_psi_block, 
                                   Zeta_real = entropy_df[f"Zeta_unsampled_{partition_label}"], 
                                   individual_var = individual_var, 
                                   sex_var = sex_var, 
                                   balance_by_var = balance_by_var, 
                                   block_label=block_label, 
                                   n_iterations=n_pval_iterations, 
                                   n_cpus=n_cpus)
            out_path = (os.path.join(save_dir, f"pvals_entropy_metrics_unsampled_{partition_label}_{block_label}.csv"))
            pvals.to_csv(out_path)
            print(f'\nSaved all unsampled entropy metrics with pvalues to {out_path}')

        
    print('\nember run complete')
            