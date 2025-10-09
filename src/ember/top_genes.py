import pandas as pd
import os

def highly_specific_to_partition(
    partition_label,
    pvals_dir,
    save_dir,
    psi_thresh = 0.5,
    zeta_thresh = 0.5,
    q_thresh = 0.05
):
    """
    Identifies significant and specific genes from a ember generated 
    p-values/q-values CSV file based on thresholds for Psi, Zeta, and q-values.

    This function reads a CSV file containing Psi and Zeta metrics (and their
    corresponding q-values), filters genes that meet given significance and
    specificity thresholds, and saves the resulting subset to a new CSV file.

    Parameters
    ----------
    pvals_dir : str, Required
        Path to the input CSV file (e.g., 'pvals_entropy_metrics_Age_E16.5.csv').
        The CSV must include the following columns:
        'Psi q-value', 'Zeta q-value', 'Psi', and 'Zeta'.

    save_dir : str, Required
        Directory where the filtered results CSV will be saved.

    partition_label : Required
        Name of partition used to generate entropy metrics, used to label saved csv. 
    
    psi_thresh : float, default = 0.5
        Threshold for Psi values. Only genes with Psi > psi_thresh are kept.

    zeta_thresh : float, Required, default = 0.5
        Threshold for Zeta values. Only genes with Zeta > zeta_thresh are kept.

    q_thresh : float, Required, default = 0.05
        Threshold for q-values. Genes are retained if both
        'Psi q-value' <= q_thresh and 'Zeta q-value' <= q_thresh.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the significant and specific genes that meet
        all threshold criteria. Also saved as "highly_specific_genes_to_{partition_label}.csv"
        in the specified save directory.
    """

    print(f"Finding highly specific genes by {partition_label}.")
    
    # Expand user paths
    pvals_dir = os.path.expanduser(pvals_dir)
    save_dir = os.path.expanduser(save_dir)

    
    # Load the input CSV file
    try:
        pvals = pd.read_csv(pvals_dir, index_col=0)
    except FileNotFoundError:
        print(f"Error: Input file not found at {pvals_dir}")
        return pd.DataFrame()
    except Exception as e:
        print(f"An error occurred while reading the CSV file: {e}")
        return pd.DataFrame()

    # Filter for significant and specific genes
    # Significant: q-values <= q_thresh
    # Specific: Psi > psi_thresh and Zeta > zeta_thresh
    try:
        significant = pvals[
            (pvals['Psi q-value'] <= q_thresh) &
            (pvals['Zeta q-value'] <= q_thresh) &
            (pvals['Psi'] > psi_thresh) &
            (pvals['Zeta'] > zeta_thresh)
        ]
    except KeyError as e:
        print(f"Error: Missing expected column in CSV file: {e}")
        return pd.DataFrame()
    
    # Sort significant genes - highest to lowest specificty
    significant = significant.sort_values(by = ['Zeta', 'Psi'], ascending = [False, False])
    
    # Save the filtered DataFrame to CSV
    os.makedirs(save_dir, exist_ok=True)
    output_path = os.path.join(save_dir, f"highly_specific_genes_to_{partition_label}.csv")
    
    try:
        significant.to_csv(output_path)
    except Exception as e:
        print(f"An error occurred while saving the results: {e}")
        return pd.DataFrame()

    print(f"Filtering complete. {len(significant)} specific genes identified.")
    print(f"Results saved to: {output_path}")

    return significant

def highly_specific_to_block(
    partition_label,
    block_label,
    pvals_dir,
    save_dir,
    psi_thresh = 0.5,
    psi_block_thresh = 0.5,
    q_thresh = 0.05
):
    """
    Identifies significant and specific genes from a ember generated 
    p-values/q-values CSV file based on thresholds for Psi, psi_block, and q-values.
    (Potential marker genes)
    
    This function reads a CSV file containing Psi and psi_block metrics (and their
    corresponding q-values), filters genes that meet given significance and
    specificity thresholds, and saves the resulting subset to a new CSV file.

    Parameters
    ----------
    pvals_dir : str, Required
        Path to the input CSV file (e.g., 'pvals_entropy_metrics_Age_E16.5.csv').
        The CSV must include the following columns:
        'Psi q-value', 'psi_block q-value', 'Psi', and 'psi_block'.

    save_dir : str, Required
        Directory where the filtered results CSV will be saved.

    partition_label : Required
        Name of partition used to generate entropy metrics, used to label saved csv. 
    
    block_label : Required
        Name of block in partition used to generate entropy metrics, used to label saved csv. 
    
    psi_thresh : float, default = 0.5
        Threshold for Psi values. Only genes with Psi > psi_thresh are kept.

    psi_block_thresh : float, Required, default = 0.5
        Threshold for psi_block values. Only genes with psi_block > psi_block_thresh are kept.

    q_thresh : float, Required, default = 0.05
        Threshold for q-values. Genes are retained if both
        'Psi q-value' <= q_thresh and 'psi_block q-value' <= q_thresh.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the significant and specific genes that meet
        all threshold criteria. Also saved as 
        "highly_specific_genes_by_{partition_label}_{block_label}.csv"
        in the specified save directory.
    """

    print(f"Finding highly specific genes by {partition_label} and {block_label}.")
    
    # Expand user paths
    pvals_dir = os.path.expanduser(pvals_dir)
    save_dir = os.path.expanduser(save_dir)

    
    # Load the input CSV file
    try:
        pvals = pd.read_csv(pvals_dir, index_col=0)
    except FileNotFoundError:
        print(f"Error: Input file not found at {pvals_dir}")
        return pd.DataFrame()
    except Exception as e:
        print(f"An error occurred while reading the CSV file: {e}")
        return pd.DataFrame()

    # Filter for significant and specific genes
    # Significant: q-values <= q_thresh
    # Specific: Psi > psi_thresh and psi_block > psi_block_thresh
    try:
        significant = pvals[
            (pvals['Psi q-value'] <= q_thresh) &
            (pvals['psi_block q-value'] <= q_thresh) &
            (pvals['Psi'] > psi_thresh) &
            (pvals['psi_block'] > psi_block_thresh)
        ]
    except KeyError as e:
        print(f"Error: Missing expected column in CSV file: {e}")
        return pd.DataFrame()
    
    # Sort significant genes - highest to lowest specificty
    significant = significant.sort_values(by = ['psi_block', 'Psi'], ascending = [False, False])
    
    # Save the filtered DataFrame to CSV
    os.makedirs(save_dir, exist_ok=True)
    output_path = os.path.join(save_dir, f"highly_specific_genes_by_{partition_label}_{block_label}.csv")
    
    try:
        significant.to_csv(output_path)
    except Exception as e:
        print(f"An error occurred while saving the results: {e}")
        return pd.DataFrame()

    print(f"Filtering complete. {len(significant)} specific genes identified.")
    print(f"Results saved to: {output_path}")

    return significant

def non_specific_to_partition(
    partition_label,
    pvals_dir,
    save_dir,
    psi_thresh = 0.5,
    zeta_thresh = 0.5,
    q_thresh = 0.05
):
    """
    Identifies significant and non-specific genes from a ember generated 
    p-values/q-values CSV file based on thresholds for Psi, Zeta, and q-values.
    (Potential housekeeping genes)

    This function reads a CSV file containing Psi and Zeta metrics (and their
    corresponding q-values), filters genes that meet given significance and
    specificity thresholds, and saves the resulting subset to a new CSV file.

    Parameters
    ----------
    pvals_dir : str, Required
        Path to the input CSV file (e.g., 'pvals_entropy_metrics_Age_E16.5.csv').
        The CSV must include the following columns:
        'Psi q-value', 'Zeta q-value', 'Psi', and 'Zeta'.

    save_dir : str, Required
        Directory where the filtered results CSV will be saved.

    partition_label : Required
        Name of partition used to generate entropy metrics, used to label saved csv. 
    
    psi_thresh : float, default = 0.5
        Threshold for Psi values. Only genes with Psi > psi_thresh are kept.

    zeta_thresh : float, Required, default = 0.5
        Threshold for Zeta values. Only genes with Zeta < zeta_thresh are kept.

    q_thresh : float, Required, default = 0.05
        Threshold for q-values. Genes are retained if both
        'Psi q-value' <= q_thresh and 'Zeta q-value' <= q_thresh.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the significant and specific genes that meet
        all threshold criteria. Also saved as "non_specific_genes_to_{partition_label}.csv"
        in the specified save directory.
    """

    print(f"Finding non-specific genes by {partition_label}.")
    
    # Expand user paths
    pvals_dir = os.path.expanduser(pvals_dir)
    save_dir = os.path.expanduser(save_dir)

    
    # Load the input CSV file
    try:
        pvals = pd.read_csv(pvals_dir, index_col=0)
    except FileNotFoundError:
        print(f"Error: Input file not found at {pvals_dir}")
        return pd.DataFrame()
    except Exception as e:
        print(f"An error occurred while reading the CSV file: {e}")
        return pd.DataFrame()

    # Filter for significant and specific genes
    # Significant: q-values <= q_thresh
    # Specific: Psi > psi_thresh and Zeta > zeta_thresh
    try:
        significant = pvals[
            (pvals['Psi q-value'] <= q_thresh) &
            (pvals['Zeta q-value'] <= q_thresh) &
            (pvals['Psi'] > psi_thresh) &
            (pvals['Zeta'] < zeta_thresh)
        ]
    except KeyError as e:
        print(f"Error: Missing expected column in CSV file: {e}")
        return pd.DataFrame()
    
    # Sort significant genes - highest to lowest specificty
    significant = significant.sort_values(by = ['Zeta', 'Psi'], ascending = [True, False])
    
    # Save the filtered DataFrame to CSV
    os.makedirs(save_dir, exist_ok=True)
    output_path = os.path.join(save_dir, f"non_specific_genes_to_{partition_label}.csv")
    
    try:
        significant.to_csv(output_path)
    except Exception as e:
        print(f"An error occurred while saving the results: {e}")
        return pd.DataFrame()

    print(f"Filtering complete. {len(significant)} non-specific genes identified.")
    print(f"Results saved to: {output_path}")

    return significant

