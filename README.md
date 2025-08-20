<p align="center">
  <img alt="How to use ember" src="./images/ember_description_sc.png" width="80%">
</p>


## `light_ember()`

This is the main workflow function in the `ember` package. It loads single-cell data, optionally performs balanced sampling across experimental replicates, computes entropy metrics ($\Psi$ and $\Zeta$), and generates permutation-based p-values to assess statistical significance.

***
### Usage Example

```python
import ember

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
````

-----

### Parameters

| Parameter | Type | Description | Default |
| :--- | :--- | :--- | :--- |
| **`adata`** | AnnData | An AnnData object already loaded in memory. Use this or `h5ad_dir`. | `None` |
| **`h5ad_dir`** | str | Path to the `.h5ad` file to process. Use this or `adata`. | `None` |
| **`partition_label`** | str | **Required.** Column in `.obs` to group cells by (e.g., 'celltype'). | `None` |
| **`save_dir`** | str | **Required.** Path to the directory where results will be saved. | `None` |
| **`sampling`** | bool | If `True`, performs balanced sampling. The next three parameters are then required. | `True` |
| `sample_id_col` | str | Column with unique sample/replicate IDs (e.g., 'sample\_id'). | `None` |
| `category_col` | str | Column defining primary groups to balance within (e.g., 'disease\_status'). | `None` |
| `condition_col` | str | Column with conditions to balance across within each category (e.g., 'sex'). | `None` |
| `num_draws` | int | Number of balanced subsets to generate if sampling. | `100` |
| `seed` | int | Random seed for reproducible sampling. | `42` |
| `partition_pvals` | bool | If `True`, computes p-values for the main partition. | `True` |
| `block_pvals` | bool | If `True`, computes p-values for a specific block within the partition. | `False` |
| `block_label` | str | A specific value from `partition_label` to use for block p-value tests. Required if `block_pvals=True`. | `None` |
| `n_pval_iterations` | int | Number of permutations for p-value calculation. | `1000` |
| `n_cpus` | int | Number of CPU cores to use for parallel processing. | `1` |

-----

### Output Files

After running, the `save_dir` will contain the following results:

  * A main CSV file with the entropy metrics for each feature `entropy_metrics_<partition_label>.csv/`.
  * A folder named `mean_Psi_block_df/` or `sampled_Psi_block_df/` containing CSVs with the detailed $\\Psi\_{block}$ values for each feature.
  * A separate CSV file with the calculated p-values for $\\Psi$ and $\\Zeta$ `pvals_entropy_metrics_<partition_label>_<block>.csv/`.

<!-- end list -->

```
```