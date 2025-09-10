import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from adjustText import adjust_text
import numpy as np
import pandas as pd
import os

def _find_top_genes(df, n_top=15, method='sum', metric_cols=None, fdr_cols=None):
    """Finds top genes based on significance and a ranking metric."""
    if metric_cols is None or fdr_cols is None:
        raise ValueError("metric_cols and fdr_cols must be provided.")

    # Filter for genes significant in both metrics
    significant = df[(df[fdr_cols[0]] < 0.05) & (df[fdr_cols[1]] < 0.05)].copy()

    # Calculate ranking score based on the specified method
    if method == 'sum':
        significant['score'] = significant[metric_cols[0]] + significant[metric_cols[1]]
    elif method == 'difference':
        significant['score'] = significant[metric_cols[0]] - significant[metric_cols[1]]
    else:
        raise ValueError("Method must be 'sum' or 'difference'.")
    
    return significant.sort_values(by='score', ascending=False).index.tolist()[:n_top]


def _annotate_genes(ax, df, genes, x_col, y_col, outline_color, outline=True, text=False, scatter_kwargs=None, text_kwargs=None):
    """Outlines and/or annotates a list of genes on the plot."""
    # Set empty dictionaries as defaults to avoid errors
    if scatter_kwargs is None:
        scatter_kwargs = {}
    if text_kwargs is None:
        text_kwargs = {}

    gene_texts = []
    for gene in genes:
        if gene in df.index:
            x, y = df.at[gene, x_col], df.at[gene, y_col]

            if outline:
                # Use a copy to combine defaults with provided kwargs
                current_scatter_kwargs = scatter_kwargs.copy()
                current_scatter_kwargs.setdefault('facecolor', 'none')
                current_scatter_kwargs.setdefault('edgecolor', outline_color)
                ax.scatter(x, y, **current_scatter_kwargs)

            if text:
                gene_texts.append(ax.text(x, y, gene, **text_kwargs))

    if text and gene_texts:
        adjust_text(
            gene_texts,
            ax=ax,
            arrowprops=dict(arrowstyle='->', color='black', lw=0.5)
        )
        
def plot_partition_specificity(partition_label,
                               pvals_dir,
                               save_dir,
                               highlight_genes=None,
                               fontsize=18,
                               custom_palette=None
                              ):
    """
    Generate a Zeta vs. Psi scatter plot to visualize partition-specific genes.

    This function reads p-value data, colors genes based on their statistical
    significance for Psi and Zeta scores, and highlights top "marker" and
    "housekeeping" genes. It also allows for custom highlighting of a
    user-provided gene list.

    Parameters
    ----------
    partition_label : str
        The label for the partition being plotted, used in the plot title.
    pvals_dir : str
        Path to the input CSV file containing p-values and scores (Psi, Zeta, FDRs).
        The CSV must have gene names as its index column.
    save_dir : str
        Path where the output plot image will be saved.
    highlight_genes : list[str], optional
        A list of gene names to highlight and annotate on the plot, by default None.
    fontsize : int, optional
        The base font size for plot labels and text, by default 18.
    custom_palette : list[str], optional
        A list of 7 hex color codes to customize the plot's color scheme.
        If None, a default palette is used.

    Returns
    -------
    None
    """
    if highlight_genes is None:
        highlight_genes = []
    if custom_palette is None:
        custom_palette = ['#CC6677', '#44AA99', '#AA4499', '#88CCEE', 
                          'black', '#332288', '#DDDDDD']

    # --- Data and Color Setup ---
    pvals_dir = os.path.expanduser(pvals_dir)
    pvals = pd.read_csv(pvals_dir, index_col=0)
    df_plot = pvals[['Psi', 'Psi FDR', 'Zeta', 'Zeta FDR']].copy()
    
    colors = {
        'psi': custom_palette[0], 'zeta': custom_palette[1], 'highlight': custom_palette[2],
        'both': custom_palette[3], 'marker': custom_palette[4], 'housekeeping': custom_palette[5],
        'none': custom_palette[6]
    }

    conditions = [
        (df_plot['Psi FDR'] < 0.05) & (df_plot['Zeta FDR'] < 0.05),
        (df_plot['Psi FDR'] < 0.05),
        (df_plot['Zeta FDR'] < 0.05)
    ]
    color_choices = [colors['both'], colors['psi'], colors['zeta']]
    df_plot['color'] = np.select(conditions, color_choices, default=colors['none'])

    # --- Plotting ---
    fig, ax = plt.subplots(figsize=(12, 10))

    for color in df_plot['color'].unique():
        subset = df_plot[df_plot['color'] == color]
        ax.scatter(subset['Psi'], subset['Zeta'], color=color, alpha=0.6, zorder=1 if color == colors['none'] else 2)

    # --- Gene Highlighting & Annotation ---
    top_markers = _find_top_genes(df_plot, method='sum', metric_cols=['Psi', 'Zeta'], fdr_cols=['Psi FDR', 'Zeta FDR'])
    top_housekeepers = _find_top_genes(df_plot, method='difference', metric_cols=['Psi', 'Zeta'], fdr_cols=['Psi FDR', 'Zeta FDR'])

    _annotate_genes(ax, df_plot, top_markers, 'Psi', 'Zeta', colors['marker'],
                    scatter_kwargs={'s': 120, 'linewidths': 1.5, 'zorder': 4})
    _annotate_genes(ax, df_plot, top_housekeepers, 'Psi', 'Zeta', colors['housekeeping'],
                    scatter_kwargs={'s': 120, 'linewidths': 1.5, 'zorder': 4})
    scatter_args = {'s': 120, 'linewidths': 2.5, 'zorder': 5}
    text_args = {'fontsize': fontsize - 2, 'color': 'black', 'weight': 'bold', 'zorder': 5}
    _annotate_genes(ax, df_plot, highlight_genes, 'Psi', 'Zeta', colors['highlight'],
                    text=True,
                    scatter_kwargs=scatter_args,
                    text_kwargs=text_args)

    # --- Text Boxes & Legend ---
    boxprops_marker = dict(boxstyle='round', facecolor='white', alpha=0.3, edgecolor=colors['marker'])
    ax.text(1.02, 0.95, 'Marker Genes:\n' + '\n'.join(top_markers), transform=ax.transAxes,
            fontsize=fontsize-4, bbox=boxprops_marker, va='top', ha='left')

    boxprops_hk = dict(boxstyle='round', facecolor='white', alpha=0.3, edgecolor=colors['housekeeping'])
    ax.text(1.02, 0.45, 'Housekeeping Genes:\n' + '\n'.join(top_housekeepers), transform=ax.transAxes,
            fontsize=fontsize-4, bbox=boxprops_hk, va='top', ha='left')
    
    handles = [
        mpatches.Patch(color=colors['none'], label='Not Significant'),
        mpatches.Patch(color=colors['psi'], label='Significant Psi Only'),
        mpatches.Patch(color=colors['zeta'], label='Significant Zeta Only'),
        mpatches.Patch(color=colors['both'], label='Significant Both'),
        mpatches.Patch(facecolor='none', edgecolor=colors['highlight'], linewidth=2, label='Highlighted Genes'),
        mpatches.Patch(facecolor='none', edgecolor=colors['marker'], linewidth=2, label='Marker Genes'),
        mpatches.Patch(facecolor='none', edgecolor=colors['housekeeping'], linewidth=2, label='Housekeeping Genes')
    ]
    ax.legend(handles=handles, fontsize=fontsize - 2, loc='lower left', frameon=True)

    # --- Final Formatting ---
    ax.set_xlabel(r'$\Psi$', fontsize=fontsize + 3)
    ax.set_ylabel(r'$\zeta$', fontsize=fontsize + 3)
    ax.set_title(fr'$\zeta$ vs $\Psi$ for {partition_label}', fontsize=fontsize + 3)
    ax.set_xticks(np.arange(0, 1.1, 0.2))
    ax.set_yticks(np.arange(0, 1.1, 0.2))
    ax.tick_params(axis='both', labelsize=fontsize - 2)
    ax.grid(True, linestyle='--', linewidth=0.5)
    fig.tight_layout(rect=[0, 0, 0.85, 1])
    
    save_dir = os.path.expanduser(save_dir)
    out_path = (os.path.join(save_dir, f"partition_specificity_scatterplot_{partition_label}.png"))
    plt.savefig(out_path)
    plt.show()
    print(f'Partition specificity plot saved to {out_path}')
    
    
def plot_block_specificity(partition_label,
                             block_label,
                             pvals_dir,
                             save_dir,
                             highlight_genes=None,
                             fontsize=18,
                             custom_palette=None
                            ):
    """
    Generate a psi_block vs. Psi scatter plot to visualize block-specific genes.

    This function reads p-value data, colors genes based on their statistical
    significance for Psi and psi_block scores, and highlights the top genes
    significant in both metrics. It also allows for custom highlighting of a
    user-provided gene list.

    Parameters
    ----------
    partition_label : str
        The label for the partition, used in the plot title.
    block_label : str
        The label for the block variable (e.g., a cell type or condition).
    pvals_dir : str
        Path to the input CSV file containing p-values and scores.
        The CSV must have gene names as its index column.
    save_dir : str
        Path where the output plot image will be saved.
    highlight_genes : list[str], optional
        A list of gene names to highlight and annotate on the plot, by default None.
    fontsize : int, optional
        The base font size for plot labels and text, by default 18.
    custom_palette : list[str], optional
        A list of 6 hex color codes to customize the plot's color scheme.
        If None, a default palette is used.

    Returns
    -------
    None
    """
    if highlight_genes is None:
        highlight_genes = []
    if custom_palette is None:
        custom_palette = ['#CC6677', '#44AA99', '#AA4499', '#88CCEE', 
                          'black', '#DDDDDD']
        
    # --- Data and Color Setup ---
    pvals_dir = os.path.expanduser(pvals_dir)
    pvals = pd.read_csv(pvals_dir, index_col=0)
    df_plot = pvals[['Psi', 'Psi FDR', 'psi_block', 'psi_block FDR']].copy()

    colors = {
        'psi': custom_palette[0], 'block': custom_palette[1], 'highlight': custom_palette[2],
        'both': custom_palette[3], 'top_genes': custom_palette[4], 'none': custom_palette[5]
    }
    
    conditions = [
        (df_plot['Psi FDR'] < 0.05) & (df_plot['psi_block FDR'] < 0.05),
        (df_plot['Psi FDR'] < 0.05),
        (df_plot['psi_block FDR'] < 0.05)
    ]
    color_choices = [colors['both'], colors['psi'], colors['block']]
    df_plot['color'] = np.select(conditions, color_choices, default=colors['none'])
    
    # --- Plotting ---
    fig, ax = plt.subplots(figsize=(12, 10))

    for color in df_plot['color'].unique():
        subset = df_plot[df_plot['color'] == color]
        ax.scatter(subset['Psi'], subset['psi_block'], color=color, alpha=0.6)

    # --- Gene Highlighting & Annotation ---
    top_genes = _find_top_genes(df_plot, n_top=10, method='sum', metric_cols=['Psi', 'psi_block'], fdr_cols=['Psi FDR', 'psi_block FDR'])
    
    _annotate_genes(ax, df_plot, top_genes, 'Psi', 'psi_block', colors['top_genes'],
                    scatter_kwargs={'s': 140, 'linewidths': 1.5, 'zorder': 4})
    scatter_args = {'s': 120, 'linewidths': 2.5, 'zorder': 5}
    text_args = {'fontsize': fontsize - 2, 'color': 'black', 'weight': 'bold', 'zorder': 5}
    _annotate_genes(ax, df_plot, highlight_genes, 'Psi', 'psi_block', colors['highlight'],
                    text=True,
                    scatter_kwargs=scatter_args,
                    text_kwargs=text_args)

    # --- Text Boxes & Legend ---
    boxprops = dict(boxstyle='round', facecolor='white', alpha=0.3, edgecolor=colors['top_genes'])
    ax.text(1.02, 0.95, 'Top Genes:\n' + '\n'.join(top_genes), transform=ax.transAxes,
            fontsize=fontsize-4, bbox=boxprops, va='top', ha='left')
    
    handles = [
        mpatches.Patch(color=colors['none'], label='Not Significant'),
        mpatches.Patch(color=colors['psi'], label='Significant Psi Only'),
        mpatches.Patch(color=colors['block'], label=f'Significant {block_label} Only'),
        mpatches.Patch(color=colors['both'], label='Significant Both'),
        mpatches.Patch(facecolor='none', edgecolor=colors['highlight'], linewidth=2, label='Highlighted Genes'),
        mpatches.Patch(facecolor='none', edgecolor=colors['top_genes'], linewidth=2, label='Top Genes (Both)')
    ]
    ax.legend(handles=handles, fontsize=fontsize - 2, loc='lower left', frameon=True)

    # --- Final Formatting ---
    ax.set_xlabel(r'$\Psi$', fontsize=fontsize + 3)
    ax.set_ylabel(fr'$\psi_{{{block_label}}}$', fontsize=fontsize + 3)
    ax.set_title(fr'$\psi_{{{block_label}}}$ vs $\Psi$ for {partition_label}', fontsize=fontsize + 3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(True, linestyle='--', linewidth=0.5)
    fig.tight_layout(rect=[0, 0, 0.85, 1])

    save_dir = os.path.expanduser(save_dir)
    out_path = (os.path.join(save_dir, f"block_specificity_scatterplot_{partition_label}_{block_label}.png"))
    plt.savefig(out_path)
    plt.show()
    print(f'Block specificity plot saved to {out_path}')
    
import os
import subprocess
import tempfile
import pandas as pd
import scanpy as sc

def check_r_package(pkg, pachterlab=True):
    """
    Checks for the presence of R and a specified R package.

    This function executes a non-interactive R command to determine if a
    package namespace can be loaded. It raises an error if the R executable
    is not found or if the specified package is not installed.

    Parameters
    ----------
    pkg (str): The name of the R package to check for.
    
    pachterlab (bool, optional): If True, the error message for a missing
        package will include a URL pointing to the pachterlab GitHub
        organization. By default, True.

    Raises
    ----------
    RuntimeError:
        If R is not installed or if the specified R package is not found.
    """
    
    try:
        result = subprocess.run(
            ["R", "--quiet", "-e", f"if (!requireNamespace('{pkg}', quietly=TRUE)) quit(status=1)"],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            if pachterlab:
                raise RuntimeError(f"R package '{pkg}' is not installed. Please see installation instructions at https://github.com/pachterlab/{pkg}.git.")
            else:
                raise RuntimeError(f"R package '{pkg}' is not installed.")
    except FileNotFoundError:
        raise RuntimeError("R is not installed.")
        

def generate_alluvial_plot(
    gene_name,
    output_path,
    columns,
    wompwomp_path,
    adata=None,
    h5ad_dir=None,
    column_weights="umi_total",
    sorting_algorithm="tsp",
    verbose=True,
    dev_mode=False,
    auto_adjust_text=True,
    match_order="None",
    axis_text_size=20,
    save_height=9,
    save_width=12
):
    """
    Generates an alluvial plot by first creating the data and then plotting it.

    This function is a one-stop-shop that first processes an AnnData object
    to calculate UMI counts for specific genes and then immediately passes the 
    resulting data to the `wompwomp`command-line tool to generate theplot.

    Parameters
    ----------
    gene_name : (str or list of str)
        The name of the gene or a list of gene names to analyze.
    output_path : (str)
        Path to save the output plot file (e.g., 'plot.pdf').
    columns : (list of str)
        Column names from `adata.obs` to use for the alluvial plot axes.
    wompwomp_path : (str)
        Path to wompwomp file (e.g., '/home/username/wompwomp/).
    adata : (sc.AnnData, optional)
        An in-memory AnnData object. One of `adata` or `h5ad_dir` must be provided.
    h5ad_dir : (str, optional)
        The file path to an .h5ad file. If provided, it's loaded in backed mode.
    column_weights : (str, optional)
        The column to use for flow weights. Defaults to 'umi_total'.
    sorting_algorithm : (str, optional)
        The sorting algorithm to use ('tsp', 'None', etc.). Defaults to 'tsp'.
    verbose : (bool, optional)
        Enables verbose output by adding the '-v' flag. Defaults to True.
    dev_mode : (bool, optional)
        Enables dev mode by adding the '--dev' flag. Defaults to True.
    auto_adjust_text : (bool, optional)
        Adds the '--auto_adjust_text' flag. Defaults to True.
    match_order : (str, optional)
        Argument for the '--match_order' flag. Defaults to 'None'.
    axis_text_size : (int, optional)
        Font size for axis text. Defaults to 20.
    save_height : (int, optional)
        Height of the saved plot in inches. Defaults to 9.
    save_width : (int, optional)
        Width of the saved plot in inches. Defaults to 12.
    """
    
    # Check for R and the wompwomp package before doing any work
    print("Step 0: Checking for R package 'wompwomp'...")
    check_r_package("wompwomp")
    print("Check successful.")

    temp_file_path = None
    try:
        # --- Part 1: Generate the DataFrame ---
        print("\nStep 1: Generating data for the alluvial plot...")

        # 1a. Load AnnData object
        if h5ad_dir is not None:
            adata_path = os.path.expanduser(h5ad_dir)
            if not os.path.exists(adata_path):
                raise FileNotFoundError(f"File not found at: {adata_path}")
            print(f'Loading AnnData from {adata_path} in backed mode.')
            adata = sc.read_h5ad(h5ad_dir, backed='r')
        elif adata is not None:
            if not isinstance(adata, sc.AnnData):
                raise TypeError("`adata` must be an AnnData object.")
            print('Using the provided in-memory AnnData object.')
        else:
            raise ValueError("You must provide either an `adata` object or an `h5ad_dir` path.")

        # 1b. Calculate UMI totals for the specified gene(s)
        alluvial_df = adata.obs[columns].copy()
        gene_adata = adata[:, gene_name].to_memory()
        alluvial_df["umi_total"] = gene_adata.X.sum(axis=1).A1
        alluvial_df['umi_total'] = alluvial_df['umi_total'] * 1000
        rounded_df = alluvial_df.round({"umi_total": 0})
        print("Data generation complete.")

        # --- Part 2: Plot the DataFrame ---
        print("\nStep 2: Plotting the generated data...")

        # 2a. Save the generated DataFrame to a temporary file
        temp_f = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        temp_file_path = temp_f.name
        rounded_df.to_csv(temp_file_path, index=False)
        temp_f.close()

        
        # 2b. Build and run the wompwomp command
        os.environ['RETICULATE_CONDA'] = '/home/nikki/miniconda3/bin/conda'
        wompwomp_executable_path = os.path.join(wompwomp_path, 'exec/wompwomp')
        wompwomp_command = [
            wompwomp_executable_path, "plot_alluvial",
            "--df", temp_file_path,
            "--graphing_columns", *columns,
            "--column_weights", column_weights,
            "--sorting_algorithm", sorting_algorithm,
            "-o", output_path,
            "--match_order", match_order,
            "--axis_text_size", str(axis_text_size),
            "--save_height", str(save_height),
            "--save_width", str(save_width)
        ]

        if verbose: wompwomp_command.append("-v")
        if dev_mode: wompwomp_command.append("--dev")
        if auto_adjust_text: wompwomp_command.append("--auto_adjust_text")

        print(f"Running command: {' '.join(wompwomp_command)}")
        subprocess.run(wompwomp_command, check=True)
        print(f"Plot saved successfully to {output_path}")

    finally:
        # --- Part 3: Clean up the temporary file ---
        if temp_file_path and os.path.exists(temp_file_path):
            print(f"Cleaning up temporary file: {temp_file_path}")
            os.remove(temp_file_path)