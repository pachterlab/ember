import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from adjustText import adjust_text
import numpy as np
import pandas as pd
import seaborn as sns
import anndata as ad
import os
import gc

def _find_top_genes(df, n_top=15, method='sum', metric_cols=None, fdr_cols=None):
    """
    Worker function that finds top genes based on significance and a ranking metric.
    """
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
    """
    Worker function that outlines and/or annotates a list of genes on the plot.
    """
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
    "housekeeping" genes. Only interpret genes that are significant by both Psi
    and Zeta since those are genes that have reliable scores after permutation 
    testing. 
    Allows for custom highlighting of a user-provided gene list.
    Fontsize and color pallette can be customized. 

    Parameters
    ----------
    partition_label : str, Required.
        The label for the partition being plotted, used in the plot title.
        
    pvals_dir : str, Required.
        Path to the input CSV file containing p-values and scores (Psi, Zeta, FDRs).
        The CSV must have gene names as its index column.
        
    save_dir : str, Required.
        Path where the output plot image will be saved.
        
    highlight_genes : list[str], default=None.
        A list of gene names to highlight and annotate on the plot, by default None.
        
    fontsize : int, default=18.
        The base font size for plot labels and text, by default 18.
        
    custom_palette : list[str], default=None.
        A list of 7 hex color codes to customize the plot's color scheme.
        If None, a default palette is used.
        Please provide list in this order 
        ['significant by psi',
        'significant by zeta',
        'highlight genes',
        'significant by both',
        'cirlce markers',
        'circle housekeeping genes',
        'significant by neither']

    Returns
    -------
    None
    """
    if highlight_genes is None:
        highlight_genes = []
    if custom_palette is None:
        custom_palette = ['#CC6677', '#44AA99', '#AA4499', '#88CCEE', 
                          'black', '#332288', '#DDDDDD']

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

    # Plotting
    fig, ax = plt.subplots(figsize=(12, 10))

    for color in df_plot['color'].unique():
        subset = df_plot[df_plot['color'] == color]
        ax.scatter(subset['Psi'], subset['Zeta'], color=color, alpha=0.6, zorder=1 if color == colors['none'] else 2)

    # Gene highlighting and annotation
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

    # Text boxes and legend
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

    # Final formatting
    ax.set_xlabel(r'$\Psi$', fontsize=fontsize + 3)
    ax.set_ylabel(r'$\zeta$', fontsize=fontsize + 3)
    ax.set_title(fr'$\zeta$ vs $\Psi$ for {partition_label}', fontsize=fontsize + 3)
    ax.set_xticks(np.arange(0, 1.1, 0.2))
    ax.set_yticks(np.arange(0, 1.1, 0.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.tick_params(axis='both', labelsize=fontsize - 2)
    ax.grid(True, linestyle='--', linewidth=0.5)
    fig.tight_layout(rect=[0, 0, 0.85, 1])
    
    save_dir = os.path.expanduser(save_dir)
    os.makedirs(save_dir, exist_ok=True)
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
    significant in both metrics. Only interpret genes that are significant by both Psi
    and psi_block since those are genes that have reliable scores after permutation 
    testing. 
    Allows for custom highlighting of a user-provided gene list.
    Fontsize and color pallette can be customized. 

    Parameters
    ----------
    partition_label : str, Required.
        The label for the partition, used in the plot title.
        
    block_label : str, Required.
        The label for the block variable (e.g., a cell type or condition).
        
    pvals_dir : str, Required.
        Path to the input CSV file containing p-values and scores.
        The CSV must have gene names as its index column.
        
    save_dir : str, Required.
        Path where the output plot image will be saved.
        
    highlight_genes : list[str], default=None.
        A list of gene names to highlight and annotate on the plot, by default None.
        
    fontsize : int, default = 18. 
        The base font size for plot labels and text, by default 18.
        
    custom_palette : list[str], default=None.
        A list of 6 hex color codes to customize the plot's color scheme.
        If None, a default palette is used.
        Provide list of colors int his order:
        ['significant by psi',
        'significant by psi_block',
        'highlight genes',
        'significant by both',
        'cirlce markers',
        'circle housekeeping genes',
        'significant by neither']

    Returns
    -------
    None
    """
    if highlight_genes is None:
        highlight_genes = []
    if custom_palette is None:
        custom_palette = ['#CC6677', '#44AA99', '#AA4499', '#88CCEE', 
                          'black', '#DDDDDD']
        
    # Setup
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
    
    # Plotting
    fig, ax = plt.subplots(figsize=(12, 10))

    for color in df_plot['color'].unique():
        subset = df_plot[df_plot['color'] == color]
        ax.scatter(subset['Psi'], subset['psi_block'], color=color, alpha=0.6)

    # Gene highlighting and annotation ---
    top_genes = _find_top_genes(df_plot, n_top=10, method='sum', metric_cols=['Psi', 'psi_block'], fdr_cols=['Psi FDR', 'psi_block FDR'])
    
    _annotate_genes(ax, df_plot, top_genes, 'Psi', 'psi_block', colors['top_genes'],
                    scatter_kwargs={'s': 140, 'linewidths': 1.5, 'zorder': 4})
    scatter_args = {'s': 120, 'linewidths': 2.5, 'zorder': 5}
    text_args = {'fontsize': fontsize - 2, 'color': 'black', 'weight': 'bold', 'zorder': 5}
    _annotate_genes(ax, df_plot, highlight_genes, 'Psi', 'psi_block', colors['highlight'],
                    text=True,
                    scatter_kwargs=scatter_args,
                    text_kwargs=text_args)

    # Text boxes and legend
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

    # Final formatting
    ax.set_xlabel(r'$\Psi$', fontsize=fontsize + 3)
    ax.set_ylabel(fr'$\psi_{{{block_label}}}$', fontsize=fontsize + 3)
    ax.set_title(fr'$\psi_{{{block_label}}}$ vs $\Psi$ for {partition_label}', fontsize=fontsize + 3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.tick_params(axis='both', labelsize=fontsize - 2)
    ax.grid(True, linestyle='--', linewidth=0.5)
    fig.tight_layout(rect=[0, 0, 0.85, 1])

    save_dir = os.path.expanduser(save_dir)
    os.makedirs(save_dir, exist_ok=True)
    out_path = (os.path.join(save_dir, f"block_specificity_scatterplot_{partition_label}_{block_label}.png"))
    plt.savefig(out_path)
    plt.show()
    print(f'Block specificity plot saved to {out_path}')
    

def plot_sample_counts(
    h5ad_dir,
    save_dir,
    sample_id_col,
    category_col,
    condition_col, 
    fontsize = 18):
    """
    Generate a bar plot showing the number of unique individuals 
    per category and condition.

    This function reads an AnnData object from an .h5ad file in backed mode, 
    calculates the number of unique individuals for each combination of a given 
    category and condition, and visualizes these counts as a grouped bar plot.
    Fontsize can be customized. 

    Parameters
    ----------
    h5ad_dir : str, Required
        Path to the input AnnData (.h5ad) file.
        
    save_dir : str, Required
        Path to directory to save the output plot image.
        
    sample_id_col : str, Required
        The column name in adata.obs that contains unique sample IDs.
        
    category_col : str, Required
        The column name to use for the primary categories on the x-axis.
        
    condition_col : str, Required
        The column name to use for grouping the bars (hue).
        
    fontsize : int, default = 18. 
        The base font size for plot labels and text, by default 18.

    Returns
    -------
    None
    
    """
    
    #Validate h5ad_dir
    adata_path = os.path.expanduser(h5ad_dir)
    if not os.path.exists(adata_path):
        raise FileNotFoundError(f"The file specified by `h5ad_dir` was not found at: {adata_path}")
    
    #Read in adata in backed mode
    adata = ad.read_h5ad(h5ad_dir, backed='r')
    
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

    df = pd.DataFrame({
        'id': adata.obs[sample_id_col].astype(str),
        'category': adata.obs[category_col].astype(str),
        'condition': adata.obs[condition_col].astype(str)
    })

    # Calculate unique counts
    counts = df.groupby(['category', 'condition'])['id'].nunique().reset_index(name='count')

    # Create the plot
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(12, 7))
    sns.barplot(data=counts, x='category', y='count', hue='condition')

    # Customize the plot
    # Set labels and title
    plt.title(f'Number of Samples by {category_col} and {condition_col}', fontsize = fontsize+2)
    plt.xlabel(category_col, fontsize = fontsize)
    plt.ylabel('Count', fontsize = fontsize)
    plt.legend(title=condition_col.title(), fontsize = fontsize-2, title_fontsize= fontsize)


    unique_categories = counts['category'].nunique()
        
    plt.xticks(ticks=range(unique_categories), labels=adata.obs[category_col].unique().astype(str), fontsize=fontsize)

    # Set y-axis ticks for clarity
    max_count = counts['count'].max()
    plt.ylim(0, max_count + 1)
    plt.yticks(range(int(max_count) + 2), fontsize=fontsize)

    plt.tight_layout()
    
    save_dir = os.path.expanduser(save_dir)
    os.makedirs(save_dir, exist_ok=True)
    out_path = (os.path.join(save_dir, f"summary_barplot.png"))

    # Save and show the plot
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {out_path}")

    plt.tight_layout() 
    plt.show()
    del(adata)
    gc.collect()
    
    

def plot_psi_blocks(
    gene_name,
    partition_label,
    psi_block_df_dir,
    save_dir, 
    fontsize = 18):
    """
    Generates and saves a bar plot of mean psi block values with error bars.

    This function reads two CSV files from a specified directory: one for mean
    psi block values and one for standard deviations. It plots the mean values
    for a specific gene as a bar plot with corresponding standard deviation
    error bars.
    Fontsize can be customized. 

    Parameters
    ----------
    gene_name : str, Required
        The name of the gene (row) to select and plot from the CSV files.
        
    partition_label : str, Required
        The partition label used to find the correct files (e.g., 'Genotype').
        
    psi_block_df_dir : str, Required
        Path to the directory containing the mean and std CSV files. Files must
        be named 'mean_Psi_block_df_{partition_label}.csv' and 
        'std_Psi_block_df_{partition_label}.csv'.
        
    save_dir : str, Required
        Path to directory to save the output plot image.
        
    fontsize : int, default=18. 
        The base font size for plot labels and text, by default 18.

    Returns
    -------
    None
    """
    
    psi_block_df_dir = os.path.expanduser(psi_block_df_dir)
    save_dir = os.path.expanduser(save_dir)
    
    # Construct paths for the mean and standard deviation files
    psi_mean_path = os.path.join(psi_block_df_dir, f'mean_Psi_block_df_{partition_label}.csv')
    psi_std_path = os.path.join(psi_block_df_dir, f'std_Psi_block_df_{partition_label}.csv')
    
    # Load mean and std data from CSV files
    try:
        df_mean = pd.read_csv(psi_mean_path, index_col=0)
        df_std = pd.read_csv(psi_std_path, index_col=0)
    except FileNotFoundError as e:
        print(f"Error: A required data file was not found. {e}")
        return
    except Exception as e:
        print(f"An error occurred while reading the CSV files: {e}")
        return

    # Select the data for the specified gene from both dataframes
    try:
        gene_mean = df_mean.loc[gene_name]
        gene_std = df_std.loc[gene_name]
    except KeyError:
        print(f"Error: Gene '{gene_name}' not found in the data files.")
        return
        
    # Sort the data alphabetically by the block names to ensure alignment
    gene_mean = gene_mean.sort_index()
    gene_std = gene_std.reindex(gene_mean.index) # Ensure std is in the same order as mean

    # Create the plot
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(12, 7))
    
    # Use plt.bar to include error bars (yerr)
    plt.bar(
        x=gene_mean.index,
        height=gene_mean.values,
        yerr=gene_std.values,
        color='grey',
        capsize=4  
    )

    # Customize the plot
    plt.title(f'Mean $\psi$_block values for {gene_name} in {partition_label}', fontsize=fontsize+2)
    plt.xlabel('Block', fontsize=fontsize)
    plt.ylabel(r'$\psi_{block}$', fontsize=fontsize)
    plt.xticks(rotation=45, ha='right', fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.grid(axis='x')
    
    # Set the y-axis limits to be between 0 and 1
    plt.ylim(0, 1)
    
    plt.tight_layout()
    
    save_dir = os.path.expanduser(save_dir)
    os.makedirs(save_dir, exist_ok=True)
    out_path = (os.path.join(save_dir, f"psi_blocks_barplot_{gene_name}_{partition_label}.png"))

    # Save and show the plot
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"Plot successfully saved to {out_path}")

    plt.show()
    
