import pandas as pd
import numpy as np
import os
import gc
import multiprocessing
import anndata as ad
from itertools import islice
from tqdm import tqdm
from scipy.sparse import csr_matrix
from scipy import sparse
from statsmodels.stats.multitest import multipletests
from .generate_entropy_metrics import generate_entropy_metrics
from .sample_replicates import generate_balanced_draws
from . import _shm_utils


def _worker_wrapper(args):
    """
    Unpacks a tuple of arguments and calls the real worker function.
    Needed because pool.imap only accepts functions with a single argument.
    """
    return _run_permutation_task(*args)


def _prepare_permutation_tasks(obs_full, var_index, sets, sample_id_col,
                                partition_label, block_label, n_iterations, shm_desc):
    """
    Generator that yields task arguments for each permutation.
    Only row indices (int array), label array, and the shared-memory descriptor are
    pickled per task — the count matrix is read from shared memory inside the worker.
    """
    sample_col = obs_full[sample_id_col]
    for i in range(n_iterations):
        current_set_index = i % len(sets)
        row_mask = sample_col.isin(sets[current_set_index]).values
        row_indices = np.where(row_mask)[0].astype(np.int64)
        labels = obs_full.loc[row_mask, partition_label].values
        yield (i, row_indices, labels, var_index, partition_label, block_label, shm_desc)


def _run_permutation_task(i, row_indices, labels, var_index, partition_label, block_label, shm_desc):
    """
    Worker that reads its X slice from shared memory, shuffles partition labels,
    and returns numpy arrays of entropy metrics.
    """
    from ._shm_utils import read_shm_csr_rows
    data_matrix = read_shm_csr_rows(shm_desc, row_indices)
    shuffled = np.random.permutation(labels)
    obs_df = pd.DataFrame({partition_label: shuffled}, index=np.arange(len(shuffled)).astype(str))
    subset = ad.AnnData(X=data_matrix, obs=obs_df, var=pd.DataFrame(index=var_index))

    Psi, Psi_block_df, Zeta = generate_entropy_metrics(subset, partition_label)

    psi_block_vals = Psi_block_df[block_label].values if block_label is not None else None

    del subset
    gc.collect()

    return np.asarray(Psi), psi_block_vals, np.asarray(Zeta)


def generate_pvals(
    h5ad_dir,
    partition_label,
    entropy_metrics_dir,
    save_dir,
    sample_id_col,
    category_col,
    condition_col,
    block_label=None,
    seed=42,
    n_iterations=1000,
    n_cpus=1,
    Psi_real=None,
    Psi_block_df_real=None,
    Zeta_real=None,
    _merge_into_entropy_file=False,
    _shm_descriptor=None,
    _obs_full=None,
    _var_index=None,
    _balanced_sets=None,
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
    - Psi q-value: Multiple testing corrected q-values for Psi scores.
    - Zeta q-value: Multiple testing corrected q-values for Zeta scores.Correction perfromed to include all p-values generated in a single file (Psi and Zeta).

    if block_pvals = True and a single block_label is given:

    - psi_block: psi_block scores (between 0 and 1) generated by light_ember for each gene in `.var`.
    - psi_block p-value: Permutation based empirical p-values for observed psi_block scores for each gene in `.var`.
    - psi_block q-value: Multiple testing corrected q-values for psi_block scores. Correction perfromed to include all p-values generated in a single file (Psi, psi_block and Zeta).

    """
    own_shm = _shm_descriptor is None
    shm_handles = None
    adata = None

    try:
        if own_shm:
            # ── Standalone path: load adata, validate, create shm ─────────────
            adata_path = os.path.expanduser(h5ad_dir)
            if not os.path.exists(adata_path):
                raise FileNotFoundError(f"The file specified by `h5ad_dir` was not found at: {adata_path}")

            entropy_path_check = os.path.expanduser(entropy_metrics_dir)
            if not os.path.exists(entropy_path_check):
                raise FileNotFoundError(f"The folder specified by `entropy_metrics_dir` was not found at: {entropy_path_check}")

            print(f'Loading AnnData object from {adata_path}.')
            adata = ad.read_h5ad(h5ad_dir)

            if partition_label not in adata.obs.columns:
                raise ValueError(
                    f"partition_label '{partition_label}' not found in adata.obs columns. "
                    f"Available columns: {list(adata.obs.columns)}"
                )

            if block_label is not None:
                valid_blocks = adata.obs[partition_label].unique()
                if block_label not in valid_blocks:
                    raise ValueError(
                        f"block_label '{block_label}' not found in adata.obs['{partition_label}']. "
                        f"Available block labels: {list(valid_blocks)}"
                    )

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

            if any(v is None for v in (Psi_real, Psi_block_df_real, Zeta_real)):
                entropy_metrics_dir_exp = os.path.expanduser(entropy_metrics_dir)
                all_entropy_metrics = pd.read_csv(
                    os.path.join(entropy_metrics_dir_exp, f"entropy_metrics_{partition_label}.csv"),
                    index_col=0
                )
                Psi_real = all_entropy_metrics[f"Psi_mean_{partition_label}"]
                Zeta_real = all_entropy_metrics[f"Zeta_mean_{partition_label}"]
                mean_Psi_block_dir = os.path.join(entropy_metrics_dir_exp, "Psi_block_df")
                Psi_block_df_real = pd.read_csv(
                    os.path.join(mean_Psi_block_dir, f"mean_Psi_block_df_{partition_label}.csv"),
                    index_col=0
                )

            sets, _ = generate_balanced_draws(
                adata, sample_id_col, category_col, condition_col,
                num_draws=n_iterations, seed=seed
            )

            # Build shm before releasing adata
            _x = adata.X
            if not sparse.issparse(_x):
                _csr = csr_matrix(np.asarray(_x, dtype=np.float32))
            else:
                _csr = csr_matrix(_x.astype(np.float32))
            obs_full = adata.obs.copy()
            var_index = adata.var.index.copy()
            del adata
            adata = None
            gc.collect()

            shm_handles, shm_desc = _shm_utils.create_shm_csr(_csr)
            del _csr
            gc.collect()

        else:
            # ── Integrated path: use data and shm from light_ember ────────────
            shm_desc = _shm_descriptor
            obs_full = _obs_full
            var_index = _var_index
            sets = _balanced_sets

            entropy_path_check = os.path.expanduser(entropy_metrics_dir)
            if not os.path.exists(entropy_path_check):
                raise FileNotFoundError(f"The folder specified by `entropy_metrics_dir` was not found at: {entropy_path_check}")

            if partition_label not in obs_full.columns:
                raise ValueError(
                    f"partition_label '{partition_label}' not found in obs columns. "
                    f"Available columns: {list(obs_full.columns)}"
                )

            if block_label is not None:
                valid_blocks = obs_full[partition_label].unique()
                if block_label not in valid_blocks:
                    raise ValueError(
                        f"block_label '{block_label}' not found in partition '{partition_label}'. "
                        f"Available block labels: {list(valid_blocks)}"
                    )

        print(f'\nGenerating p-values for entropy metrics.')

        # Filter to genes with valid observed Psi
        mask = Psi_real > 0
        Psi_real = Psi_real[mask]
        Psi_block_df_real = Psi_block_df_real[mask]
        Zeta_real = Zeta_real[mask]

        gene_idx = Psi_real.index
        n_genes = len(gene_idx)
        # gene_idx is a filtered subset (Psi > 0); get_indexer maps it to positions
        # in the full permutation arrays returned by workers (which cover all genes).
        gene_positions = var_index.get_indexer(gene_idx)

        psi_real_vals = Psi_real.values
        zeta_real_vals = Zeta_real.values
        psi_block_real_vals = (
            Psi_block_df_real[block_label].values if block_label is not None else None
        )

        psi_count_ge = np.zeros(n_genes)
        zeta_count_ge = np.zeros(n_genes)
        psi_block_count_ge = np.zeros(n_genes) if block_label is not None else None

        # Adaptive stopping: check convergence every 100 iters after a 300-iter
        # burn-in (enough iterations for p-value estimates to stabilise before
        # the first comparison). Stop when no p-value shifts by more than 0.002.
        min_iters = 300
        check_interval = 100
        tol = 0.002
        prev_pvals = None
        n_done = 0

        task_gen = _prepare_permutation_tasks(
            obs_full, var_index, sets, sample_id_col,
            partition_label, block_label, n_iterations, shm_desc
        )

        print(f"\nStarting up to {n_iterations} parallel permutations with {n_cpus} workers...")
        with multiprocessing.Pool(processes=n_cpus) as pool:
            with tqdm(total=n_iterations, desc="Permutations") as pbar:
                while n_done < n_iterations:
                    batch_n = min(check_interval, n_iterations - n_done)
                    batch_tasks = list(islice(task_gen, batch_n))
                    if not batch_tasks:
                        break

                    # imap_unordered: result order doesn't matter (we only accumulate counts),
                    # and it avoids head-of-line blocking when one task is slower than others.
                    for psi_arr, psi_block_arr, zeta_arr in pool.imap_unordered(_worker_wrapper, batch_tasks):
                        psi_count_ge  += (psi_arr[gene_positions]  >= psi_real_vals)
                        zeta_count_ge += (zeta_arr[gene_positions] >= zeta_real_vals)
                        if psi_block_count_ge is not None:
                            psi_block_count_ge += (psi_block_arr[gene_positions] >= psi_block_real_vals)
                        n_done += 1

                    pbar.update(len(batch_tasks))

                    if n_done >= min_iters:
                        curr = np.concatenate([
                            (psi_count_ge  + 1) / (n_done + 1),
                            (zeta_count_ge + 1) / (n_done + 1),
                        ])
                        if prev_pvals is not None:
                            max_change = float(np.nanmax(np.abs(curr - prev_pvals)))
                            if max_change < tol:
                                print(
                                    f"\nAdaptive stopping: converged after {n_done} iterations "
                                    f"(max Δp = {max_change:.4f})."
                                )
                                break
                        prev_pvals = curr

        # +1 Laplace pseudo-count keeps all p-values in (0, 1] and prevents
        # zero p-values from dominating FDR correction when n_iterations is small.
        Psi_p_values  = pd.Series((psi_count_ge  + 1) / (n_done + 1), index=gene_idx)
        Zeta_p_values = pd.Series((zeta_count_ge + 1) / (n_done + 1), index=gene_idx)

        common_index = gene_idx

        if block_label is not None:
            psi_block_real = Psi_block_df_real[block_label]
            psi_block_p_values = pd.Series(
                (psi_block_count_ge + 1) / (n_done + 1), index=gene_idx
            )
            final = pd.DataFrame({
                'Psi':                 pd.Series(Psi_real,           index=common_index),
                'Psi p-value':         pd.Series(Psi_p_values,       index=common_index),
                'psi_block':           pd.Series(psi_block_real,     index=common_index),
                'psi_block p-value':   pd.Series(psi_block_p_values, index=common_index),
                'Zeta':                pd.Series(Zeta_real,          index=common_index),
                'Zeta p-value':        pd.Series(Zeta_p_values,      index=common_index),
            }, index=common_index)
        else:
            final = pd.DataFrame({
                'Psi':         pd.Series(Psi_real,     index=common_index),
                'Psi p-value': pd.Series(Psi_p_values, index=common_index),
                'Zeta':        pd.Series(Zeta_real,    index=common_index),
                'Zeta p-value': pd.Series(Zeta_p_values, index=common_index),
            }, index=common_index)

        final.index.name = 'gene_name'

        # Global multiple testing correction across all p-value columns at once
        pval_cols = [c for c in final.columns if c.lower().endswith('p-value')]
        parts = []
        for col in pval_cols:
            s = final[col]
            valid = s[s.notna() & np.isfinite(s.values)]
            parts.append(valid.rename(col))

        stacked = pd.concat(parts, keys=pval_cols)
        _, qvals, _, _ = multipletests(stacked.values, method='fdr_bh')
        qval_series = pd.Series(qvals, index=stacked.index)

        for col in pval_cols:
            fdr_col = col.replace('p-value', 'q-value')
            final[fdr_col] = np.nan
            if col in qval_series.index.get_level_values(0):
                col_qvals = qval_series[col]
                final.loc[col_qvals.index, fdr_col] = col_qvals.values

        save_dir_exp = os.path.expanduser(save_dir)
        os.makedirs(save_dir_exp, exist_ok=True)

        if _merge_into_entropy_file:
            entropy_path = os.path.join(save_dir_exp, f"entropy_metrics_{partition_label}.csv")
            entropy_df = pd.read_csv(entropy_path, index_col=0)
            merge_cols = [c for c in final.columns if 'p-value' in c.lower() or 'q-value' in c.lower()]
            if block_label is not None:
                merge_cols = ['psi_block'] + merge_cols
            # Drop any pre-existing p/q-value columns to avoid join collision on re-runs
            entropy_df = entropy_df.drop(columns=[c for c in merge_cols if c in entropy_df.columns], errors='ignore')
            entropy_df = entropy_df.join(final[merge_cols], how='left')
            entropy_df.to_csv(entropy_path)
            print(f'\nMerged p-values into {entropy_path}')
        else:
            if block_label is not None:
                out_path = os.path.join(save_dir_exp, f"pvals_entropy_metrics_{partition_label}_{block_label}.csv")
            else:
                out_path = os.path.join(save_dir_exp, f"pvals_entropy_metrics_{partition_label}.csv")
            final.to_csv(out_path)
            print(f'\nSaved all entropy metrics along with pvalues to {out_path}')

    finally:
        if adata is not None:
            del adata
            gc.collect()
        if own_shm and shm_handles is not None:
            _shm_utils.release_shm(shm_handles)
