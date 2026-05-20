import os
import gc
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import multiprocessing
from tqdm import tqdm
from scipy.sparse import csr_matrix
from scipy import sparse
from .generate_entropy_metrics import generate_entropy_metrics
from .generate_pvals import generate_pvals
from .sample_replicates import generate_balanced_draws, aitchison_mean_and_std
from . import _shm_utils


def _worker_wrapper(args):
    """
    Unpacks a tuple of arguments and calls the real worker function.
    Needed because pool.imap only accepts functions with a single argument.
    """
    return _run_computation_on_subset(*args)


def _prepare_task_args(sets, sample_id_col, obs_full, var_index, partition_label, shm_desc):
    """
    Generator that yields task arguments using shared memory for the X matrix.
    Only the partition-label column and row indices are pickled per task;
    the count matrix is read from shared memory inside the worker.
    """
    sample_col = obs_full[sample_id_col]
    for i in range(len(sets)):
        row_mask = sample_col.isin(sets[i]).values
        row_indices = np.where(row_mask)[0].astype(np.int64)
        labels = obs_full.loc[row_mask, partition_label].values
        yield (i, row_indices, labels, var_index, partition_label, shm_desc)


def _run_computation_on_subset(i, row_indices, labels, var_index, partition_label, shm_desc):
    """
    Worker that reads its X slice from shared memory, builds a minimal AnnData,
    and returns entropy metrics.
    """
    from ._shm_utils import read_shm_csr_rows
    data_matrix = read_shm_csr_rows(shm_desc, row_indices)

    obs_df = pd.DataFrame({partition_label: labels}, index=np.arange(len(labels)).astype(str))
    subset = ad.AnnData(
        X=data_matrix,
        obs=obs_df,
        var=pd.DataFrame(index=var_index),
    )

    temp_psi, temp_psi_block, temp_zeta = generate_entropy_metrics(subset, partition_label)
    del subset
    gc.collect()

    return i, np.asarray(temp_psi), temp_psi_block, np.asarray(temp_zeta)


def light_ember(
    h5ad_dir,
    partition_label,
    save_dir,
    sampling=True,
    sample_id_col=None,
    category_col=None,
    condition_col=None,
    num_draws=100,
    save_draws=False,
    seed=42,
    partition_pvals=True,
    block_pvals=False,
    block_label=None,
    n_pval_iterations=1000,
    n_cpus=1
):
    """
    Runs the ember entropy metrics and p-value generation workflow on an AnnData object.

    This function loads an AnnData `.h5ad` file, optionally performs balanced sampling
    across replicates, computes entropy metrics for the specified partition,
    and generates p-values for Psi and Zeta and optionally Psi_block for a block of choice.

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

    save_dir : str, Required
        Path to directory where results will be saved. Required to run process.

    sampling : bool, default=True
        Whether to perform balanced sampling across replicates before entropy calculation.
        If True, `sample_id_col`, `category_col`, and `condition_col` must be provided.
        Sampling should only be False if fast intermediate results are desired or
        if there are no replicates to sample over.
        If sampling is set to False but either partition_pvals or block_pvals are set to
        True then the sampling=False will be overridden as pval generation requires sampling.

    sample_id_col : str, default = None
        The column in `.obs` with unique identifiers for each sample or replicate
        (e.g., 'sample_id', 'mouse_id').

    category_col : str, default = None
        The column in `.obs` defining the primary group to balance across in order
        to generate a balanced sample of the experiment. (e.g., 'disease_status', 'mouse_strain').
        Refer to readme for further explanation on how to select category and condition columns.
        category_col and condition_col are interchangable.
        If balancing across more than 2 variables, generate interaction terms, create a column
        in `.obs` concatnating the two (or more) columns of interested with a semicolon (:).
        This way balancing can be done across as many variables as desired.

    condition_col : str, default = None
        The column in `.obs` containing the conditions to balance within
        each categoryto generate a balanced sample of the experiment.  (e.g., 'sex', 'treatment').
        Refer to readme for further explanation on how to select category and condition columns.
        category_col and condition_col are interchangable.
        If balancing across more than 2 variables, generate interaction terms, create a column
        in `.obs` concatnating the two (or more) columns of interested with a semicolon (:).
        This way balancing can be done across as many variables as desired.

    num_draws : int, default = 100
        The number of balanced subsets to generate, by default 100.

    save_draws : bool, default=False
        Whether to save intermediate draws to save_dir.

    seed : int,  default = 42
        The random seed for reproducible draws, by default 42.

    partition_pvals : bool, default=True
        Whether to compute permutation-based p-values for the `partition_label`.
        P-values are generated by sampling. If sampling = False and partition_pvals = True,
        the sampling=False will be overwritten.
        Calls generate_pavls, which can be called manually after metric generation as well.

    block_pvals : bool, default=False
        Whether to compute permutation-based p-values for the `block_label`.
        P-values are generated by sampling. If sampling = False and block_pvals = True,
        the sampling=False will be overwritten.
        Calls generate_pavls, which can be called manually after metric generation as well.

    block_label : str,  default = None
        One value in the `.obs` column for partition_label to use for block-based permutation tests.
        Required if `block_pvals=True`.

    n_pval_iterations : int, default=1000
        Number of permutations to use for p-value calculation.

    n_cpus : int, default=1
        Number of CPU cores to use for parallel permutation testing.
        For this script, performance is I/O-bound and may not improve beyond 4-8 cores.'

    Returns
    -------
    None

    Notes
    -----
    - Results are saved to `save_dir` as CSV files.
    - one csv file with all entropy metrics
    - one csv file in a new Psi_block_df folder with psi block values for all blocks in a partition
    - Separate file for pvals
    - Separate files for each partition
    - Alternate file names depending on sampling on or off.

    **What to expect inside 'entropy_metrics.csv'**:

    - gene_name: All genes in `.var`
    - Psi_mean: Psi scores averaged over n draws (between 0 and 1) corresponding to the selected partition for each gene in `.var`.
    - Psi_std: Standard deviation of Psi scores across n draws corresponding to the selected partition for each gene in `.var`.
    - Psi_valid_counts: Number of valid Psi scores observed across n draws. Only use genes for downstream analysis that have valid counts=num_draws. If valid counts is not close to num_draws, increase threshold for filtering genes with low reads beforehand(recommended <100 reads, increase as needed).
    - Zeta_mean: Zeta scores averaged over n draws (between 0 and 1) corresponding to the selected partition for each gene in `.var`.
    - Zeta_std: Standard deviation of Zeta scores across n draws corresponding to the selected partition for each gene in `.var`.
    - Zeta_valid_counts: Number of valid Psi scores observed across n draws. Only use genes for downstream analysis that have valid counts=num_draws. If valid counts is not close to num_draws, increase threshold for filtering genes with low reads beforehand (recommended <100 reads, increase as needed).

    **What to expect inside 'Psi_block_df/'**:

    - mean_Psi_block_df.csv : A dataframe of mean Psi_block scores (between 0 and 1) corresponding to the selected partition for each gene in `.var`. Scores are caluclated for all blocks, each column of the dataframe corresponds to one block.
    - std_Psi_block_df.csv : A dataframe of standard deviations for Psi_block scores corresponding to the selected partition for each gene in `.var`.Scores are caluclated for all blocks, each column of the dataframe corresponds to one block.

    **What to expect inside 'pvals_entropy_metrics.csv'**:

    - gene_name: All genes in `.var`
    - Psi: Psi scores averaged over n draws (between 0 and 1) generated by light_ember for each gene in `.var`.
    - Psi p-value: Permutation based empirical p-values for observed Psi scores for each gene in `.var`.
    - Zeta: Zeta scores averaged over n draws (between 0 and 1) generated by light_ember for each gene in `.var`.
    - Zeta p-value: Permutation based empirical p-values for observed Zeta scores for each gene in `.var`.
    - Psi q-value: Multiple testing corrected q-values for Psi scores.
    - Zeta q-value: Multiple testing corrected q-values for Zeta scores.Correction perfromed to include all p-values generated in a single file (Psi and Zeta).

    If block_pvals = True and a single block_label is given:

    - psi_block: psi_block scores (between 0 and 1) generated by light_ember for each gene in `.var`.
    - psi_block p-value: Permutation based empirical p-values for observed psi_block scores for each gene in `.var`.
    - psi_block q-value: Multiple testing corrected q-values for psi_block scores. Correction perfromed to include all p-values generated in a single file (Psi, psi_block and Zeta).

    """

    adata_path = os.path.expanduser(h5ad_dir)
    if not os.path.exists(adata_path):
        raise FileNotFoundError(f"The file specified by `h5ad_dir` was not found at: {adata_path}")

    adata = None
    shm_handles = None

    try:
        print(f'Loading AnnData object from {adata_path}.')
        adata = ad.read_h5ad(h5ad_dir)

        # Validate partition_label
        if partition_label not in adata.obs.columns:
            raise ValueError(
                f"partition_label '{partition_label}' not found in adata.obs columns. "
                f"Available columns: {list(adata.obs.columns)}"
            )

        # Validate sample_id_col, category_col, and condition_col if sampling or p-value calculation is enabled.
        if sampling or partition_pvals or block_pvals:
            required_params = {
                "sample_id_col": sample_id_col,
                "category_col": category_col,
                "condition_col": condition_col
            }
            if not all(required_params.values()):
                missing_params = [name for name, value in required_params.items() if not value]
                raise ValueError(
                    f"If sampling or generating p-values, you must provide `sample_id_col`, `category_col`, and `condition_col`. "
                    f"You are missing the following parameter(s): {missing_params}"
                )

            required_cols = list(required_params.values())
            missing_cols = [col for col in required_cols if col not in adata.obs.columns]
            if missing_cols:
                raise ValueError(
                    f"The following required column(s) were not found in adata.obs: {missing_cols}. "
                    f"\nAvailable columns are: {list(adata.obs.columns)}"
                )

        # Override sampling = False if pvals are true.
        if (partition_pvals or block_pvals) and not sampling:
            print("\np-value generation requires sampling. sampling=False overriden")
            sampling = True

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
        elif not block_pvals:
            block_label = None

        if sampling:
            sets, usage = generate_balanced_draws(
                adata, sample_id_col, category_col, condition_col, num_draws=num_draws, seed=seed
            )
            print(f'\nGenerating entropy metrics.')

            # Generate pval sets from adata while it is still in memory
            sets_pval = None
            if partition_pvals or block_pvals:
                sets_pval, _ = generate_balanced_draws(
                    adata, sample_id_col, category_col, condition_col,
                    num_draws=n_pval_iterations, seed=seed
                )

            # Build shared memory from adata.X before releasing adata
            _x = adata.X
            if not sparse.issparse(_x):
                _csr = csr_matrix(np.asarray(_x, dtype=np.float32))
            else:
                _csr = csr_matrix(_x.astype(np.float32))
            obs_full = adata.obs.copy()
            var_index_saved = adata.var.index.copy()
            del adata
            adata = None
            gc.collect()

            shm_handles, shm_desc = _shm_utils.create_shm_csr(_csr)
            del _csr
            gc.collect()

            if save_draws:
                draws_dir = os.path.join(os.path.expanduser(save_dir), "balanced_draws")
                os.makedirs(draws_dir, exist_ok=True)
                print(f'\nDraws will be saved to: {draws_dir}.')

            task_generator = _prepare_task_args(
                sets, sample_id_col, obs_full, var_index_saved, partition_label, shm_desc
            )

            print(f"\nComputing {num_draws} iterations in parallel with {n_cpus} workers")
            results = []
            with multiprocessing.Pool(processes=n_cpus) as pool:
                for result in tqdm(pool.imap_unordered(_worker_wrapper, task_generator), total=len(sets)):
                    results.append(result)
                    if save_draws:
                        i, temp_psi, temp_psi_block, temp_zeta = result
                        draw_dir = os.path.join(draws_dir, f"BALANCED_{i:02d}")
                        os.makedirs(draw_dir, exist_ok=True)
                        pd.DataFrame({
                            f"Psi_{partition_label}": temp_psi,
                            f"Zeta_{partition_label}": temp_zeta,
                        }, index=temp_psi_block.index).to_csv(
                            os.path.join(draw_dir, f"entropy_metrics_{partition_label}.csv")
                        )
                        temp_psi_block.to_csv(os.path.join(draw_dir, f"Psi_block_{partition_label}.csv"))

            # ========== Aggregation Phase ==========
            save_dir = os.path.expanduser(save_dir)
            os.makedirs(save_dir, exist_ok=True)
            print(f'\nAggregating entropy metrics from {num_draws} samples.')

            all_Psi = [r[1] for r in results]
            Psi_block_list = [r[2] for r in results]
            all_Zeta = [r[3] for r in results]
            gene_names = results[0][2].index

            common_blocks = sorted(set.intersection(*(set(df.columns) for df in Psi_block_list)))
            Psi_block_list = [df[common_blocks] for df in Psi_block_list]

            Psi_arr = np.stack(all_Psi)
            Zeta_arr = np.stack(all_Zeta)
            mask = Psi_arr != -1

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)

                aggregate_entropy_df = pd.DataFrame({
                    f"Psi_mean_{partition_label}":        np.nanmean(np.where(mask, Psi_arr, np.nan), axis=0),
                    f"Psi_std_{partition_label}":         np.nanstd( np.where(mask, Psi_arr, np.nan), axis=0),
                    f"Psi_valid_counts_{partition_label}": np.sum(mask, axis=0),
                    f"Zeta_mean_{partition_label}":       np.nanmean(np.where(mask, Zeta_arr, np.nan), axis=0),
                    f"Zeta_std_{partition_label}":        np.nanstd( np.where(mask, Zeta_arr, np.nan), axis=0),
                    f"Zeta_valid_counts_{partition_label}": np.sum(mask, axis=0)
                }, index=gene_names)

            aggregate_entropy_df.loc[
                aggregate_entropy_df[f"Psi_valid_counts_{partition_label}"] == 0,
                [f"Psi_mean_{partition_label}", f"Psi_std_{partition_label}"]
            ] = np.nan

            aggregate_entropy_df.loc[
                aggregate_entropy_df[f"Zeta_valid_counts_{partition_label}"] == 0,
                [f"Zeta_mean_{partition_label}", f"Zeta_std_{partition_label}"]
            ] = np.nan

            mean_Psi_block_dir = os.path.join(save_dir, "Psi_block_df")
            os.makedirs(mean_Psi_block_dir, exist_ok=True)
            aitchison_results = aitchison_mean_and_std(Psi_block_list)
            mean_Psi_block = pd.DataFrame(aitchison_results[0], index=gene_names, columns=common_blocks)
            std_Psi_block  = pd.DataFrame(aitchison_results[1], index=gene_names, columns=common_blocks)
            mean_Psi_block.to_csv(os.path.join(mean_Psi_block_dir, f"mean_Psi_block_df_{partition_label}.csv"))
            std_Psi_block.to_csv( os.path.join(mean_Psi_block_dir, f"std_Psi_block_df_{partition_label}.csv"))

            aggregate_entropy_df.to_csv(os.path.join(save_dir, f"entropy_metrics_{partition_label}.csv"))
            print(f'\nSaved all entropy metrics to {os.path.join(save_dir, f"entropy_metrics_{partition_label}.csv")}')

            if partition_pvals or block_pvals:
                generate_pvals(
                    h5ad_dir=h5ad_dir,
                    partition_label=partition_label,
                    entropy_metrics_dir=save_dir,
                    save_dir=save_dir,
                    sample_id_col=sample_id_col,
                    category_col=category_col,
                    condition_col=condition_col,
                    block_label=block_label,
                    n_iterations=n_pval_iterations,
                    n_cpus=n_cpus,
                    Psi_real=aggregate_entropy_df[f"Psi_mean_{partition_label}"],
                    Psi_block_df_real=mean_Psi_block,
                    Zeta_real=aggregate_entropy_df[f"Zeta_mean_{partition_label}"],
                    _merge_into_entropy_file=True,
                    _shm_descriptor=shm_desc,
                    _obs_full=obs_full,
                    _var_index=var_index_saved,
                    _balanced_sets=sets_pval,
                )

        else:
            # sampling=False: adata is still in memory
            print("\nSampling turned off, generating entropy metrics without sampling.")
            temp_psi, temp_psi_block, temp_zeta = generate_entropy_metrics(adata, partition_label)
            var_idx = adata.var.index

            save_dir = os.path.expanduser(save_dir)
            os.makedirs(save_dir, exist_ok=True)

            entropy_df = pd.DataFrame({
                f"Psi_unsampled_{partition_label}": temp_psi,
                f"Zeta_unsampled_{partition_label}": temp_zeta,
            }, index=var_idx)
            entropy_df.to_csv(os.path.join(save_dir, f"entropy_metrics_unsampled_{partition_label}.csv"))

            unsampled_Psi_block_dir = os.path.join(save_dir, "Psi_block_df")
            os.makedirs(unsampled_Psi_block_dir, exist_ok=True)
            temp_psi_block.to_csv(os.path.join(unsampled_Psi_block_dir,
                                               f"Psi_block_unsampled_{partition_label}.csv"))

            print(f'\nSaved all unsampled entropy metrics to {os.path.join(save_dir, f"entropy_metrics_unsampled_{partition_label}.csv")}')

            del adata
            adata = None
            gc.collect()

        print('\nember run complete')

    finally:
        if adata is not None:
            del adata
            gc.collect()
        if shm_handles is not None:
            _shm_utils.release_shm(shm_handles)
