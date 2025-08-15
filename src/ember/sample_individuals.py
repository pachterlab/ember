from collections import defaultdict
import random
import math
import numpy as np
import pandas as pd
from scipy import stats

def generate_balanced_sets(adata, individual_var, sex_var, balance_by_var, num_draws=100, seed=42):
    """
    Generate balanced subsets of individuals based on sex and a balance-by variable.
    Each draw includes 1 male and 1 female from each balance_by category.
    Rotates evenly over available individuals within each (category, sex) group.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object with `.obs` containing metadata.
    individual_var : str
        Column in `.obs` that uniquely identifies individuals (e.g., mouse_id, sample_id).
    sex_var : str
        Column in `.obs` that uniquely identifies sex (e.g., Sex, sex).
    balance_by_var : str
        Column in `.obs` to blance draws by (e.g., mouse_strain, genotype, category, disease).
    num_draws : int
        Number of subsets to generate.
    seed : int
        Random seed for reproducibility.
    
    Returns
    -------
    subset_sample_ids : list of lists
        A list of subsets, each containing sample IDs.
    used_count : dict
        Dictionary tracking how often each individual and group combination was used.
        
    """

    df = adata.obs
    # Group individuals by (balance_by, sex)
    grouped = (
        df.groupby([balance_by_var, sex_var])[individual_var]
        .unique()
        .apply(list)
        .to_dict()
    )

    all_balance_vals = sorted(set(df[balance_by_var]))
    all_sex_vals = sorted(set(df[sex_var]))

    used_count = defaultdict(int)
    subset_sample_indivs = []

    random.seed(seed)

    for _ in range(num_draws):
        selected_inds = []

        for b in all_balance_vals:
            for s in all_sex_vals:
                if (b, s) not in grouped or not grouped[(b, s)]:
                    continue

                # Pick the least-used individual from this combination
                choices = sorted(grouped[(b, s)], key=lambda x: used_count[x])
                chosen = choices[0]
                selected_inds.append(chosen)

                # Update usage counts
                used_count[chosen] += 1

        # Append selected individual IDs
        subset_sample_indivs.append(selected_inds)
    print(f'\nSampling {num_draws} out of {math.prod([len(v) for v in grouped.values()])} possible combinations of individuals')

    return subset_sample_indivs, used_count


# Function to combine sampled values
def aitchison_mean_and_std(Psi_block_dfs_list):
    """
    Compute the Aitchison mean of compositional DataFrames along with standard deviation.

    Parameters
    ----------
    Psi_block_dfs_list : list of pd.DataFrame
        List of compositional matrices (genes x choice_of_partition).
        
    Returns
    -------
    mean_df : pd.DataFrame
        Aitchison mean of Psi block values.
    var_df : pd.DataFrame
        Std of the Psi block values.
    """
    epsilon = 1e-10
    gene_names = Psi_block_dfs_list[0].index
    strain_names = Psi_block_dfs_list[0].columns

    # === Aitchison mean ===
    data_stack = np.stack([df.values for df in Psi_block_dfs_list], axis=0) + epsilon
    gmean_data = stats.gmean(data_stack, axis=0)
    gmean_df = pd.DataFrame(gmean_data, index=gene_names, columns=strain_names)
    row_sums = gmean_df.sum(axis=1).replace(0, np.nan)
    normalized_df = gmean_df.div(row_sums, axis=0).fillna(0)

    # === STD ===
    normalized_list = []
    for df in Psi_block_dfs_list:
        norm_df = df.div(df.sum(axis=1), axis=0).fillna(0)
        normalized_list.append(norm_df.values)

    norm_stack = np.stack(normalized_list, axis=0)
    var_data = np.std(norm_stack, axis=0, ddof=1)
   
    var_df = pd.DataFrame(var_data, index=gene_names, columns=strain_names)

    return normalized_df, var_df