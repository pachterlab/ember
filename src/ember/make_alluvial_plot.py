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
        wompwomp_command = [
            "./exec/wompwomp", "plot_alluvial",
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