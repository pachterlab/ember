<p align="center">
  <img alt="How to use ember" src="./images/ember_description_sc.png" width="80%">
</p>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ember Function Documentation</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji";
            line-height: 1.6;
            color: #333;
            max-width: 800px;
            margin: 20px auto;
            padding: 0 15px;
        }
        h2, h3 {
            border-bottom: 1px solid #eaecef;
            padding-bottom: 0.3em;
            margin-top: 24px;
            margin-bottom: 16px;
        }
        code {
            font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, Courier, monospace;
            background-color: #f6f8fa;
            padding: 0.2em 0.4em;
            margin: 0;
            font-size: 85%;
            border-radius: 3px;
        }
        pre {
            background-color: #f6f8fa;
            border-radius: 5px;
            padding: 16px;
            overflow: auto;
        }
        pre code {
            padding: 0;
            margin: 0;
            font-size: 100%;
            background-color: transparent;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 1em;
            margin-bottom: 1em;
            display: block;
            overflow: auto;
        }
        th, td {
            border: 1px solid #dfe2e5;
            padding: 8px 12px;
            text-align: left;
        }
        th {
            background-color: #f6f8fa;
            font-weight: 600;
        }
        ul {
            padding-left: 20px;
        }
    </style>
</head>
<body>

    <h2><code>light_ember()</code></h2>

    <p>
        This is the main workflow function in the <code>ember</code> package. It loads single-cell data, optionally performs balanced sampling across experimental replicates, computes entropy metrics (&Psi; and &Zeta;), and generates permutation-based p-values to assess statistical significance.
    </p>

    <h3>Usage Example</h3>
<pre><code>import ember

ember.light_ember(
    h5ad_dir='path/to/my_data.h5ad',
    save_dir='path/to/results/',
    partition_label='celltype',
    sampling=True,
    sample_id_col='sample_id',
    category_col='disease_status',
    condition_col='sex',
    n_cpus=8
)
</code></pre>

    <h3>Parameters</h3>
    <table>
        <thead>
            <tr>
                <th>Parameter</th>
                <th>Type</th>
                <th>Description</th>
                <th>Default</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td><strong><code>adata</code></strong></td>
                <td>AnnData</td>
                <td>An AnnData object already loaded in memory. Use this or <code>h5ad_dir</code>.</td>
                <td><code>None</code></td>
            </tr>
            <tr>
                <td><strong><code>h5ad_dir</code></strong></td>
                <td>str</td>
                <td>Path to the <code>.h5ad</code> file to process. Use this or <code>adata</code>.</td>
                <td><code>None</code></td>
            </tr>
            <tr>
                <td><strong><code>partition_label</code></strong></td>
                <td>str</td>
                <td><strong>Required.</strong> Column in <code>.obs</code> to group cells by (e.g., 'celltype').</td>
                <td><code>None</code></td>
            </tr>
            <tr>
                <td><strong><code>save_dir</code></strong></td>
                <td>str</td>
                <td><strong>Required.</strong> Path to the directory where results will be saved.</td>
                <td><code>None</code></td>
            </tr>
            <tr>
                <td><strong><code>sampling</code></strong></td>
                <td>bool</td>
                <td>If <code>True</code>, performs balanced sampling. The next three parameters are then required.</td>
                <td><code>True</code></td>
            </tr>
            <tr>
                <td><code>sample_id_col</code></td>
                <td>str</td>
                <td>Column with unique sample/replicate IDs (e.g., 'sample_id').</td>
                <td><code>None</code></td>
            </tr>
            <tr>
                <td><code>category_col</code></td>
                <td>str</td>
                <td>Column defining primary groups to balance within (e.g., 'disease_status').</td>
                <td><code>None</code></td>
            </tr>
            <tr>
                <td><code>condition_col</code></td>
                <td>str</td>
                <td>Column with conditions to balance across within each category (e.g., 'sex').</td>
                <td><code>None</code></td>
            </tr>
            <tr>
                <td><code>num_draws</code></td>
                <td>int</td>
                <td>Number of balanced subsets to generate if sampling.</td>
                <td><code>100</code></td>
            </tr>
            <tr>
                <td><code>seed</code></td>
                <td>int</td>
                <td>Random seed for reproducible sampling.</td>
                <td><code>42</code></td>
            </tr>
            <tr>
                <td><code>partition_pvals</code></td>
                <td>bool</td>
                <td>If <code>True</code>, computes p-values for the main partition.</td>
                <td><code>True</code></td>
            </tr>
            <tr>
                <td><code>block_pvals</code></td>
                <td>bool</td>
                <td>If <code>True</code>, computes p-values for a specific block within the partition.</td>
                <td><code>False</code></td>
            </tr>
            <tr>
                <td><code>block_label</code></td>
                <td>str</td>
                <td>A specific value from <code>partition_label</code> to use for block p-value tests. Required if <code>block_pvals=True</code>.</td>
                <td><code>None</code></td>
            </tr>
            <tr>
                <td><code>n_pval_iterations</code></td>
                <td>int</td>
                <td>Number of permutations for p-value calculation.</td>
                <td><code>1000</code></td>
            </tr>
            <tr>
                <td><code>n_cpus</code></td>
                <td>int</td>
                <td>Number of CPU cores to use for parallel processing.</td>
                <td><code>1</code></td>
            </tr>
        </tbody>
    </table>

    <h3>Output Files</h3>
    <p>After running, the <code>save_dir</code> will contain the following results:</p>
    <ul>
        <li>A main CSV file with the entropy metrics for each feature (e.g., <code>sampled_ember_metrics.csv</code>).</li>
        <li>A folder named <code>Psi_block_df/</code> containing CSVs with the detailed &Psi;<sub>block</sub> values for each feature.</li>
        <li>A separate CSV file with the calculated p-values for &Psi; and &Zeta;.</li>
    </ul>

</body>
</html>