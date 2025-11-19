from collections import defaultdict
import random
import math
import numpy as np
import pandas as pd
from scipy import stats

def generate_balanced_draws(adata, sample_id_col, category_col, condition_col, num_draws=100, seed=42):
    """
    Generate balanced subsets of replicates based on a categorical variable 
    (eg: Genotype, Age, Estrus Stage) and balances draws within each category
    by a condition variable (eg: Sex, Knockout/Wildtype/Overexpressed, Treated/Control). 
 
    For each unique value in `category_col`, this function creates draws where each
    unique value from `condition_col` is represented exactly once. It iterates
    through the available samples for each category-condition pair to generate
    diverse subsets across multiple draws.
    
    Parameters
    ----------
    adata : AnnData
        The annotated data object with sample metadata in `.obs`.
    sample_id_col : str
        The column in `.obs` with unique identifiers for each sample or replicate
        (e.g., 'sample_id', 'mouse_id').
    category_col : str
        The column in `.obs` defining the primary groups to balance within
        (e.g., 'disease_status', 'mouse_strain').
    condition_col : str
        The column in `.obs` containing the conditions to balance across within
        each category (e.g., 'sex', 'treatment').
    num_draws : int, optional
        The number of balanced subsets to generate, by default 100.
    seed : int, optional
        The random seed for reproducible draws, by default 42.

    Returns
    -------
    list[list[str]]
        A list of draws, where each draw is a list of sample IDs.
    dict[str, int]
        A dictionary tracking the selection count for each sample.

    """
    
    df = adata.obs
    
    # Group samples by (category, condition)
    grouped = (
        df.groupby([category_col, condition_col], observed=True)[sample_id_col]
        .unique()
        .apply(list)
        .to_dict()
    )

    all_category_vals = sorted(df[category_col].unique())
    all_condition_vals = sorted(df[condition_col].unique())

    used_count = defaultdict(int)
    balanced_draws = []

    random.seed(seed)

    for _ in range(num_draws):
        current_draw = []
        for category in all_category_vals:
            for condition in all_condition_vals:
                
                choices = grouped.get((category, condition), [])
                if not choices:
                    raise ValueError(
                        f"Cannot create balanced draws. No samples found for the combination "
                        f"({category}, {condition}). Please ensure every category-condition "
                        "pair has at least one sample."
                    )

                # Find all samples with the minimum usage count for this group
                min_count = min(used_count[c] for c in choices)
                least_used_choices = [c for c in choices if used_count[c] == min_count]
                
                # Randomly select from the least-used candidates
                chosen_sample = random.choice(least_used_choices)
                current_draw.append(chosen_sample)
                
                # Update usage count for the chosen sample
                used_count[chosen_sample] += 1
                
        balanced_draws.append(current_draw)

    total_combinations = math.prod([len(v) for v in grouped.values() if v])
    print(f'\nSampled {num_draws} out of {total_combinations} possible unique combinations.')

    return balanced_draws, dict(used_count)


# Function to combine sampled values
def aitchison_mean_and_std(Psi_block_dfs_list):
    """
    Compute the Aitchison mean and geometric standard deviation of compositional 
    DataFrames across blocks.

    Parameters
    ----------
    Psi_block_dfs_list : list of pd.DataFrame
        List of compositional matrices (genes x choice_of_partition).
        
    Returns
    -------
    mean_df : pd.DataFrame
        Aitchison mean of Psi block values.
    var_df : pd.DataFrame
        geometric std of the Psi block values.
    """
    epsilon = 1e-10
    gene_names = Psi_block_dfs_list[0].index
    block_names = Psi_block_dfs_list[0].columns

    with np.errstate(divide='ignore'):
        # === Aitchison mean ===
        data_stack = np.stack([df.values for df in Psi_block_dfs_list], axis=0) + epsilon
        gmean_data = stats.gmean(data_stack, axis=0)
        gmean_df = pd.DataFrame(gmean_data, index=gene_names, columns=block_names)
        row_sums = gmean_df.sum(axis=1).replace(0, np.nan)
        normalized_df = gmean_df.div(row_sums, axis=0).fillna(0)

        # === Geometric std ===
        log_data = np.log(np.where(data_stack > 0, data_stack, epsilon))
        log_std = np.std(log_data, axis=0, ddof=1)
        gsd_data = np.exp(log_std)
        gsd_df = pd.DataFrame(gsd_data, index=gene_names, columns=block_names)

    return normalized_df, gsd_df