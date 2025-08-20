<p align="center">
  <img alt="How to use ember" src="./images/ember_description_sc.png" width="100%">
</p>


## `light_ember()`

This is the main workflow function in `ember`. It loads with an in-memory anndata object or reads in an h5ad file from `h5ad_dir` and computes entropy metrics (&Psi;, &Zeta; and &Psi;<sub>block</sub> ) for a given partition label, optionally performs balanced sampling across experimental replicates, and generates permutation-based p-values to assess statistical significance of entropy metrics.

***

### Installation
```python
!pip install git+https://<enter github PAT>@github.com/pachterlab/ember.git

````

### Usage Example

```python
import ember
from ember.light_ember import light_ember
from ember.generate_pvals import generate_pvals

light_ember(
    h5ad_dir = 'subset_dev_kidney.h5ad',
    partition_label = 'Age',
    save_dir = '/content/output',
    sampling=True,
    sample_id_col= 'sample_id',
    category_col='Age',
    condition_col='sex',
    num_draws=50,
    partition_pvals=False,
    block_pvals=False,
    n_cpus=2
)
````

-----

### Parameters

| Parameter | Type | Description | Default |
| :--- | :--- | :--- | :--- |
| **`adata`** | AnnData | **Either adata or h5ad_dir required.** An AnnData object already loaded in memory. Use this or `h5ad_dir`. In-memory mode - faster with sufficient RAM| `None` |
| **`h5ad_dir`** | str | **Either adata or h5ad_dir required.** Path to the `.h5ad` file to process. Use this or `adata`. Memory efficient mode - slower but will run with lower compute power.| `None` |
| **`partition_label`** | str | **Required.** Column in `.obs` to group cells by (e.g., 'celltype, genotype'). | `None` |
| **`save_dir`** | str | **Required.** Path to the directory where results will be saved. | `None` |
| `sampling` | bool | If `True`, performs balanced sampling. The next three parameters are then required. | `True` |
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
  * A separate CSV file with the calculated p-values for &Psi;, &Zeta; and &Psi;<sub>block</sub> `pvals_entropy_metrics_<partition_label>_<block>.csv/`.

<!-- end list -->

```
```