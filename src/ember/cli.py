import argparse
import sys
from .light_ember import light_ember
from .generate_pvals import generate_pvals
from .plots import plot_partition_specificity, plot_block_specificity, plot_sample_counts, plot_psi_blocks

def main():
    """Main CLI entry point for the ember toolkit."""
    parser = argparse.ArgumentParser(
        prog="ember",
        description="A command-line toolkit for ember: Entropy Metrics for Biological ExploRation."
    )
    
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available sub-commands")

    
    # =================================================================
    # ==                 COMMAND 1: light_ember                      ==
    # =================================================================
    light_ember_parser = subparsers.add_parser(
    "light_ember",
    help="Runs the ember entropy metrics and p-value generation workflow on an AnnData object.",
    formatter_class=argparse.RawTextHelpFormatter,
    description="""\
    Runs the ember entropy metrics and p-value generation workflow on an AnnData object.

    This function loads an AnnData `.h5ad` file, optionally performs balanced sampling
    across replicates, computes entropy metrics for the specified partition,
    and generates p-values for Psi and Zeta and optionally Psi_block for a block of choice.

    Entropy metrics generated:
        - Psi : Fraction of information explained by partition of choice
        - Psi_block : Specificity of information to a block
        - Zeta : Specificity to a partition / distance of Psi_blocks distribution from uniform

    Notes
    -----
    - Results are saved to `save_dir` as CSV files.
    - One CSV file with all entropy metrics.
    - One CSV file in a new Psi_block_df folder with Psi_block values for all blocks in a partition.
    - Separate file for p-values.
    - Separate files for each partition.
    - Alternate file names depending on sampling on or off.
    """
    )

    # --- Required Positional Arguments ---
    light_ember_parser.add_argument(
        "h5ad_dir",
        help="Path to the `.h5ad` file to process. Data should be log1p and depth normalized "
             "before running ember. Remove genes with <100 reads before running ember."
    )
    light_ember_parser.add_argument(
        "partition_label",
        help="Column in `.obs` used to partition cells for entropy calculations "
             "(e.g., 'celltype', 'Genotype', 'Age'). For interaction terms, create a new "
             "column concatenating multiple `.obs` columns with a semicolon (:)."
    )
    light_ember_parser.add_argument(
        "save_dir",
        help="Path to directory where results will be saved."
    )

    # --- Sampling Arguments ---
    sampling_group = light_ember_parser.add_argument_group('Sampling Parameters')
    sampling_group.add_argument(
        "--no_sampling",
        action="store_false",
        dest="sampling",
        help="Disable balanced sampling. Default: True. "
             "Note: If partition_pvals or block_pvals are enabled, sampling will be re-enabled."
    )
    sampling_group.add_argument(
        "--sample_id_col",
        type=str,
        default=None,
        help="Column in `.obs` with unique identifiers for each sample or replicate "
             "(e.g., 'sample_id', 'mouse_id')."
    )
    sampling_group.add_argument(
        "--category_col",
        type=str,
        default=None,
        help="Column in `.obs` defining the primary group to balance across "
             "(e.g., 'disease_status', 'mouse_strain'). Interchangeable with condition_col. "
             "For >2 variables, create interaction terms by concatenating columns with `:`."
    )
    sampling_group.add_argument(
        "--condition_col",
        type=str,
        default=None,
        help="Secondary column in `.obs` to balance sampling across (e.g., 'sex', 'treatment'). "
             "Interchangeable with category_col. Supports interaction terms."
    )
    sampling_group.add_argument(
        "--num_draws",
        type=int,
        default=100,
        help="Number of balanced subsets to generate (default: 100)."
    )
    sampling_group.add_argument(
        "--save_draws",
        action="store_true",
        help="Save intermediate sampled draws to save_dir (default: False)."
    )
    sampling_group.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible draws (default: 42)."
    )

    # --- P-value Arguments ---
    pval_group = light_ember_parser.add_argument_group('P-value Parameters')
    pval_group.add_argument(
        "--no_partition_pvals",
        action="store_false",
        dest="partition_pvals",
        help="Disable permutation p-value calculation for the main partition. Default: True."
    )
    pval_group.add_argument(
        "--block_pvals",
        action="store_true",
        help="Enable permutation p-value calculation for a specific block. Default: False."
    )
    pval_group.add_argument(
        "--block_label",
        type=str,
        default=None,
        help="Specific value in 'partition_label' for block p-values. Required if --block_pvals is set."
    )
    pval_group.add_argument(
        "--n_pval_iterations",
        type=int,
        default=1000,
        help="Number of permutations for p-value calculation (default: 1000)."
    )

    # --- Performance Arguments ---
    perf_group = light_ember_parser.add_argument_group('Performance Parameters')
    perf_group.add_argument(
        "--n_cpus",
        type=int,
        default=1,
        help="Number of CPU cores to use for parallel processing (default: 1). "
             "Performance is I/O-bound and may not improve beyond 4–8 cores."
    )


    # =================================================================
    # ==                 COMMAND 2: generate_pvals                   ==
    # =================================================================
    generate_pvals_parser = subparsers.add_parser(
    "generate_pvals",
    help="Calculate empirical p-values for entropy metrics from permutation test results.",
    formatter_class=argparse.RawTextHelpFormatter,
    description="""\
    Calculate empirical p-values for entropy metrics from permutation test results.

    Entropy metrics generated:
        - Psi : Fraction of information explained by partition of choice
        - Psi_block : Specificity of information to a block
        - Zeta : Specificity to a partition / distance of Psi_blocks distribution from uniform
    """
    )

    # --- Required Positional Arguments ---
    generate_pvals_parser.add_argument(
        "h5ad_dir",
        help="Path to the `.h5ad` file to process. Data should be log1p and depth normalized "
             "before running ember. Remove genes with <100 reads before running ember."
    )
    generate_pvals_parser.add_argument(
        "partition_label",
        help="Column in `.obs` used to partition cells for entropy calculations "
             "(e.g., 'celltype', 'Genotype', 'Age'). For interaction terms, create a new "
             "column concatenating multiple `.obs` columns with a semicolon (:)."
    )
    generate_pvals_parser.add_argument(
        "entropy_metrics_dir",
        help="Path to CSV with entropy metrics to use for generating p-values."
    )
    generate_pvals_parser.add_argument(
        "save_dir",
        help="Path to directory where results will be saved."
    )

    # --- Sampling Arguments ---
    sampling_group = generate_pvals_parser.add_argument_group('Sampling Parameters')
    sampling_group.add_argument(
        "--sample_id_col",
        type=str,
        default=None,
        help="Column in `.obs` with unique identifiers for each sample or replicate "
             "(e.g., 'sample_id', 'mouse_id')."
    )
    sampling_group.add_argument(
        "--category_col",
        type=str,
        default=None,
        help="Column in `.obs` defining the primary group to balance across "
             "(e.g., 'disease_status', 'mouse_strain'). Interchangeable with condition_col. "
             "For >2 variables, create interaction terms by concatenating columns with `:`."
    )
    sampling_group.add_argument(
        "--condition_col",
        type=str,
        required=True,
        help="Column in `.obs` containing the conditions to balance within each category "
             "(e.g., 'sex', 'treatment'). Interchangeable with category_col. Supports interaction terms."
    )

    # --- Block Argument ---
    generate_pvals_parser.add_argument(
        "--block_label",
        type=str,
        default=None,
        help="Block in partition to calculate p-values for. Default: None (Psi and Zeta only)."
    )

    # --- Performance & Iterations ---
    perf_group = generate_pvals_parser.add_argument_group('Performance Parameters')
    perf_group.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible draws (default: 42)."
    )
    perf_group.add_argument(
        "--n_iterations",
        type=int,
        default=1000,
        help="Number of iterations to calculate p-values (default: 1000). "
             "Use fewer for quick runs, more for reliable results."
    )
    perf_group.add_argument(
        "--n_cpus",
        type=int,
        default=1,
        help="Number of CPUs to use for p-value calculation (default: 1). "
             "Set to -1 to use all available cores but one."
    )
    
    # --- Internal-use Arguments ---
    internal_group = generate_pvals_parser.add_argument_group('Internal Arguments (used by light_ember)')
    internal_group.add_argument(
        "--Psi_real",
        type=str,
        default=None,
        help="Observed Psi values for each gene (pd.Series). Not required for user runs."
    )
    internal_group.add_argument(
        "--Psi_block_df_real",
        type=str,
        default=None,
        help="Observed Psi_block values for all blocks in chosen partition (pd.DataFrame). "
             "Not required for user runs."
    )
    internal_group.add_argument(
        "--Zeta_real",
        type=str,
        default=None,
        help="Observed Zeta values for each gene (pd.Series). Not required for user runs."
    )

    

    # =================================================================
    # ==          COMMAND 3: plot_partition_specificity            ==
    # =================================================================
    plot_partition_specificity_parser = subparsers.add_parser(
    "plot_partition_specificity",
    help="Generate a Zeta vs. Psi scatter plot to visualize partition-specific genes.",
    formatter_class=argparse.RawTextHelpFormatter,
    description="""\
    Generate a Zeta vs. Psi scatter plot to visualize partition-specific genes.

    This function reads p-value data, colors genes based on their statistical
    significance for Psi and Zeta scores, and highlights top "marker" and
    "housekeeping" genes. Allows for custom highlighting of a user-provided
    gene list. Font size and color palette can be customized.
    """
    )

    # --- Required Positional Arguments ---
    plot_partition_specificity_parser.add_argument(
        "partition_label",
        help="Label for the partition being plotted, used in the plot title."
    )
    plot_partition_specificity_parser.add_argument(
        "pvals_dir",
        help="Path to input CSV containing p-values and scores (Psi, Zeta, FDRs). "
             "CSV must have gene names as its index."
    )
    plot_partition_specificity_parser.add_argument(
        "save_dir",
        help="Path where the output plot image will be saved."
    )

    # --- Optional Arguments ---
    plot_partition_specificity_parser.add_argument(
        "--highlight_genes",
        nargs="+",
        default=None,
        help="List of gene names to highlight and annotate on the plot (default: None)."
    )
    plot_partition_specificity_parser.add_argument(
        "--fontsize",
        type=int,
        default=18,
        help="Base font size for plot labels and text (default: 18)."
    )
    plot_partition_specificity_parser.add_argument(
        "--custom_palette",
        nargs="+",
        default=None,
        help="List of 7 hex color codes to customize the color scheme. Order:\n"
             "['significant by psi', 'significant by zeta', 'highlight genes', "
             "'significant by both', 'circle markers', 'circle housekeeping genes', "
             "'significant by neither']. Default: None (uses built-in palette)."
    )


    # =================================================================
    # ==            COMMAND 4: plot_block_specificity              ==
    # =================================================================
    plot_block_specificity_parser = subparsers.add_parser(
    "plot_block_specificity",
    help="Generate a psi_block vs. Psi scatter plot to visualize block-specific genes.",
    formatter_class=argparse.RawTextHelpFormatter,
    description="""\
    Generate a psi_block vs. Psi scatter plot to visualize block-specific genes.

    This function reads p-value data, colors genes based on their statistical
    significance for Psi and psi_block scores, and highlights the top genes
    significant in both metrics. Allows for custom highlighting of a user-provided
    gene list. Font size and color palette can be customized.
    """
    )

    # --- Required Positional Arguments ---
    plot_block_specificity_parser.add_argument(
        "partition_label",
        help="Label for the partition, used in the plot title."
    )
    plot_block_specificity_parser.add_argument(
        "block_label",
        help="Label for the block variable (e.g., a cell type or condition)."
    )
    plot_block_specificity_parser.add_argument(
        "pvals_dir",
        help="Path to input CSV containing p-values and scores. "
             "CSV must have gene names as its index."
    )
    plot_block_specificity_parser.add_argument(
        "save_dir",
        help="Path where the output plot image will be saved."
    )

    # --- Optional Arguments ---
    plot_block_specificity_parser.add_argument(
        "--highlight_genes",
        nargs="+",
        default=None,
        help="List of gene names to highlight and annotate on the plot (default: None)."
    )
    plot_block_specificity_parser.add_argument(
        "--fontsize",
        type=int,
        default=18,
        help="Base font size for plot labels and text (default: 18)."
    )
    plot_block_specificity_parser.add_argument(
        "--custom_palette",
        nargs="+",
        default=None,
        help="List of 6 hex color codes to customize the color scheme. Order:\n"
             "['significant by psi', 'significant by psi_block', 'highlight genes', "
             "'significant by both', 'circle markers', 'circle housekeeping genes', "
             "'significant by neither']. Default: None (uses built-in palette)."
    )


    # =================================================================
    # ==               COMMAND 5: plot_sample_counts                 ==
    # =================================================================
    plot_sample_counts_parser = subparsers.add_parser(
    "plot_sample_counts",
    help="Generate a bar plot showing the number of unique individuals per category and condition.",
    formatter_class=argparse.RawTextHelpFormatter,
    description="""\
    Generate a bar plot showing the number of unique individuals per category and condition.

    This function reads an AnnData object from an .h5ad file in backed mode, 
    calculates the number of unique individuals for each combination of a given 
    category and condition, and visualizes these counts as a grouped bar plot.
    Font size can be customized.
    """
    )

    # --- Required Positional Arguments ---
    plot_sample_counts_parser.add_argument(
        "h5ad_dir",
        help="Path to the input AnnData (.h5ad) file."
    )
    plot_sample_counts_parser.add_argument(
        "save_dir",
        help="Path to directory to save the output plot image."
    )
    plot_sample_counts_parser.add_argument(
        "sample_id_col",
        help="Column name in `.obs` that contains unique sample IDs."
    )
    plot_sample_counts_parser.add_argument(
        "category_col",
        help="Column name to use for the primary categories on the x-axis."
    )
    plot_sample_counts_parser.add_argument(
        "condition_col",
        help="Column name to use for grouping the bars (hue)."
    )

    # --- Optional Arguments ---
    plot_sample_counts_parser.add_argument(
        "--fontsize",
        type=int,
        default=18,
        help="Base font size for plot labels and text (default: 18)."
    )


    # =================================================================
    # ==                 COMMAND 6: plot_psi_blocks                  ==
    # =================================================================
    plot_psi_blocks_parser = subparsers.add_parser(
    "plot_psi_blocks",
    help="Generate a bar plot of mean psi block values with error bars for a given gene.",
    formatter_class=argparse.RawTextHelpFormatter,
    description="""\
    Generates and saves a bar plot of mean psi block values with error bars.

    This function reads two CSV files from a specified directory: one for mean
    psi block values and one for standard deviations. It plots the mean values
    for a specific gene as a bar plot with corresponding standard deviation
    error bars. Font size can be customized.
    """
    )

    # --- Required Positional Arguments ---
    plot_psi_blocks_parser.add_argument(
        "gene_name",
        help="Name of the gene (row) to select and plot from the CSV files."
    )
    plot_psi_blocks_parser.add_argument(
        "partition_label",
        help="Partition label used to find the correct files (e.g., 'Genotype')."
    )
    plot_psi_blocks_parser.add_argument(
        "psi_block_df_dir",
        help="Directory containing the mean and std CSV files. Files must be named "
             "'mean_Psi_block_df_{partition_label}.csv' and "
             "'std_Psi_block_df_{partition_label}.csv'."
    )
    plot_psi_blocks_parser.add_argument(
        "save_dir",
        help="Path to directory to save the output plot image."
    )

    # --- Optional Arguments ---
    plot_psi_blocks_parser.add_argument(
        "--fontsize",
        type=int,
        default=18,
        help="Base font size for plot labels and text (default: 18)."
    )


    # =================================================================
    # ==                  DISPATCH LOGIC (Call functions)            ==
    # =================================================================
    args = parser.parse_args()

    kwargs = vars(args)
    command_to_run = kwargs.pop('command')


    if command_to_run == "light_ember":
        light_ember(**kwargs)

    elif command_to_run == "generate_pvals":
        generate_pvals(**kwargs)

    elif command_to_run == "plot_partition_specificity":
        plot_partition_specificity(**kwargs)

    elif command_to_run == "plot_block_specificity":
        plot_block_specificity(**kwargs)

    elif command_to_run == "plot_sample_counts":
        plot_sample_counts(**kwargs)

    elif command_to_run == "plot_psi_blocks":
        plot_psi_blocks(**kwargs)


if __name__ == "__main__":
    main()