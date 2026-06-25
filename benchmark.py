"""
Speed comparison: old vs new ember implementations.
Benchmarks each of the four optimizations independently on synthetic data.
"""

import time
import tempfile
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.sparse import random as sparse_random, csr_matrix
from scipy import sparse
from statsmodels.stats.multitest import multipletests


# ── Synthetic data ────────────────────────────────────────────────────────────

def make_adata(n_cells=4000, n_genes=8000, n_blocks=6, density=0.05, seed=0):
    rng = np.random.default_rng(seed)
    X = sparse_random(n_cells, n_genes, density=density, format="csr",
                      data_rvs=lambda s: rng.exponential(scale=2, size=s))
    X = X.astype(np.float32)
    block_labels = rng.choice([f"Block{i}" for i in range(n_blocks)], size=n_cells)
    obs = pd.DataFrame({"partition": block_labels}, index=[f"cell{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"gene{i}" for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


# ── Helper (shared) ───────────────────────────────────────────────────────────

def safe_log2_sparse(mat):
    log_mat = mat.copy()
    log_mat.data = np.log2(log_mat.data)
    log_mat.data[~np.isfinite(log_mat.data)] = 0
    return log_mat

def safe_divide_sparse(numerator, denominator):
    if sparse.issparse(numerator):
        denom_safe = np.copy(denominator)
        denom_safe[np.isclose(denom_safe, 0)] = np.inf
        result = numerator.multiply(1.0 / denom_safe)
        result.eliminate_zeros()
        return result
    else:
        with np.errstate(divide='ignore', invalid='ignore'):
            result = np.divide(numerator, denominator)
            result = np.where(np.isclose(denominator, 0), 0, result)
        return result


# ══════════════════════════════════════════════════════════════════════════════
# Benchmark 1 – Disk I/O round-trip vs in-memory
# ══════════════════════════════════════════════════════════════════════════════

def bench_disk_vs_memory(n_draws=60, n_genes=8000, n_blocks=6, reps=3):
    """Measure write-N-CSVs-then-read-all vs hold results in memory."""
    rng = np.random.default_rng(42)
    gene_index = [f"gene{i}" for i in range(n_genes)]
    blocks = [f"Block{i}" for i in range(n_blocks)]

    # Pre-generate draw results (numpy arrays + DataFrames)
    draw_results = []
    for i in range(n_draws):
        psi = rng.random(n_genes).astype(np.float32)
        zeta = rng.random(n_genes).astype(np.float32)
        psi_block = pd.DataFrame(rng.random((n_genes, n_blocks)).astype(np.float32),
                                 index=gene_index, columns=blocks)
        draw_results.append((i, psi, psi_block, zeta))

    # Old: write CSVs then read back
    times_old = []
    for _ in range(reps):
        with tempfile.TemporaryDirectory() as tmpdir:
            t0 = time.perf_counter()
            for i, psi, psi_block, zeta in draw_results:
                draw_dir = os.path.join(tmpdir, f"BALANCED_{i:02d}")
                os.makedirs(draw_dir)
                pd.DataFrame({"Psi": psi, "Zeta": zeta}, index=gene_index).to_csv(
                    os.path.join(draw_dir, "entropy.csv"))
                psi_block.to_csv(os.path.join(draw_dir, "psi_block.csv"))

            all_psi, all_zeta, block_list = [], [], []
            for i in range(n_draws):
                draw_dir = os.path.join(tmpdir, f"BALANCED_{i:02d}")
                df = pd.read_csv(os.path.join(draw_dir, "entropy.csv"), index_col=0)
                sb = pd.read_csv(os.path.join(draw_dir, "psi_block.csv"), index_col=0)
                all_psi.append(df["Psi"].values)
                all_zeta.append(df["Zeta"].values)
                block_list.append(sb)
            times_old.append(time.perf_counter() - t0)

    # New: accumulate directly from in-memory results
    times_new = []
    for _ in range(reps):
        t0 = time.perf_counter()
        all_psi = [r[1] for r in draw_results]
        block_list = [r[2] for r in draw_results]
        all_zeta = [r[3] for r in draw_results]
        times_new.append(time.perf_counter() - t0)

    return np.mean(times_old), np.mean(times_new)


# ══════════════════════════════════════════════════════════════════════════════
# Benchmark 2 – Block loop: adata[mask].X vs counts[mask]
# ══════════════════════════════════════════════════════════════════════════════

def bench_block_loop(n_cells=4000, n_genes=8000, n_blocks=6, reps=5):
    adata = make_adata(n_cells, n_genes, n_blocks)
    blocks = adata.obs["partition"].unique()
    counts = adata.X
    if not sparse.issparse(counts):
        counts = csr_matrix(counts)
    total_counts_per_gene = np.asarray(counts.sum(axis=0)).ravel()

    def compute_block_entropy_old(adata, blocks, total_counts_per_gene):
        n_genes = adata.shape[1]
        n_blocks = len(blocks)
        E_W = np.zeros(n_genes)
        Psi_block_num = np.zeros((n_genes, n_blocks))
        for idx, block in enumerate(blocks):
            mask = adata.obs["partition"] == block
            block_counts = adata[mask, :].X              # OLD: AnnData subsetting
            block_sum = np.asarray(block_counts.sum(axis=0)).ravel()
            q_j = safe_divide_sparse(block_counts, block_sum)
            log_q_j = safe_log2_sparse(q_j)
            entropy = -q_j.multiply(log_q_j).sum(axis=0).A1
            p_c_j = np.divide(block_sum, total_counts_per_gene,
                               out=np.zeros_like(block_sum),
                               where=~np.isclose(total_counts_per_gene, 0))
            weighted_entropy = entropy * p_c_j
            E_W += weighted_entropy
            Psi_block_num[:, idx] = weighted_entropy

    def compute_block_entropy_new(counts, obs_partition, blocks, total_counts_per_gene):
        n_genes = counts.shape[1]
        n_blocks = len(blocks)
        E_W = np.zeros(n_genes)
        Psi_block_num = np.zeros((n_genes, n_blocks))
        for idx, block in enumerate(blocks):
            mask = obs_partition == block             # NEW: numpy array mask
            block_counts = counts[mask, :]            # NEW: direct sparse slice
            block_sum = np.asarray(block_counts.sum(axis=0)).ravel()
            q_j = safe_divide_sparse(block_counts, block_sum)
            log_q_j = safe_log2_sparse(q_j)
            entropy = -q_j.multiply(log_q_j).sum(axis=0).A1
            p_c_j = np.divide(block_sum, total_counts_per_gene,
                               out=np.zeros_like(block_sum),
                               where=~np.isclose(total_counts_per_gene, 0))
            weighted_entropy = entropy * p_c_j
            E_W += weighted_entropy
            Psi_block_num[:, idx] = weighted_entropy

    obs_partition = adata.obs["partition"].values

    times_old, times_new = [], []
    for _ in range(reps):
        t0 = time.perf_counter()
        compute_block_entropy_old(adata, blocks, total_counts_per_gene)
        times_old.append(time.perf_counter() - t0)

        t0 = time.perf_counter()
        compute_block_entropy_new(counts, obs_partition, blocks, total_counts_per_gene)
        times_new.append(time.perf_counter() - t0)

    return np.mean(times_old), np.mean(times_new)


# ══════════════════════════════════════════════════════════════════════════════
# Benchmark 3 – q-value assignment: row-by-row loop vs vectorized
# ══════════════════════════════════════════════════════════════════════════════

def bench_qvalue_loop(n_genes=15000, reps=5):
    rng = np.random.default_rng(7)
    gene_index = [f"gene{i}" for i in range(n_genes)]
    psi_pvals  = rng.uniform(0, 1, n_genes)
    zeta_pvals = rng.uniform(0, 1, n_genes)
    pval_cols  = ["Psi p-value", "Zeta p-value"]

    def make_final():
        return pd.DataFrame({
            "Psi p-value":  psi_pvals,
            "Zeta p-value": zeta_pvals,
        }, index=gene_index)

    def qval_old(final, pval_cols):
        records = []
        for col in pval_cols:
            s = final[col]
            mask_valid = s.notna() & np.isfinite(s.values)
            for idx in final.index[mask_valid]:
                records.append((idx, col, float(final.at[idx, col])))
        all_pvals = np.array([r[2] for r in records], dtype=float)
        _, qvals, _, _ = multipletests(all_pvals, method='fdr_bh')
        for col in pval_cols:
            fdr_col = col.replace('p-value', 'q-value')
            if fdr_col not in final.columns:
                final[fdr_col] = np.nan
        for (idx, col, _), q in zip(records, qvals):
            fdr_col = col.replace('p-value', 'q-value')
            final.at[idx, fdr_col] = q

    def qval_new(final, pval_cols):
        parts = []
        for col in pval_cols:
            s = final[col]
            valid = s[s.notna() & np.isfinite(s.values)]
            parts.append(valid.rename(col))
        stacked = pd.concat(parts, keys=pval_cols)
        _, qvals, _, _ = multipletests(stacked.values, method='fdr_bh')
        qval_series = pd.Series(qvals, index=stacked.index)
        for col in pval_cols:
            fdr_col = col.replace('p-value', 'q-value')
            final[fdr_col] = np.nan
            if col in qval_series.index.get_level_values(0):
                col_qvals = qval_series[col]
                final.loc[col_qvals.index, fdr_col] = col_qvals.values

    times_old, times_new = [], []
    for _ in range(reps):
        f = make_final()
        t0 = time.perf_counter()
        qval_old(f, pval_cols)
        times_old.append(time.perf_counter() - t0)

        f = make_final()
        t0 = time.perf_counter()
        qval_new(f, pval_cols)
        times_new.append(time.perf_counter() - t0)

    return np.mean(times_old), np.mean(times_new)


# ══════════════════════════════════════════════════════════════════════════════
# Benchmark 4 – E_T cache path: adata.var write vs always compute
# ══════════════════════════════════════════════════════════════════════════════

def bench_et_cache(n_cells=4000, n_genes=8000, reps=10):
    adata = make_adata(n_cells, n_genes)
    counts = adata.X
    if not sparse.issparse(counts):
        counts = csr_matrix(counts)
    total_counts_per_gene = np.asarray(counts.sum(axis=0)).ravel()

    def compute_et_old(adata, counts, total_counts_per_gene):
        if 'E_T' not in adata.var.columns:
            p_i = safe_divide_sparse(counts, total_counts_per_gene)
            log_p_i = safe_log2_sparse(p_i)
            entropy_per_gene = p_i.multiply(log_p_i).sum(axis=0).A1
            E_T = -entropy_per_gene
            E_T[np.isclose(total_counts_per_gene, 0)] = -1
            adata.var['E_T'] = E_T
        else:
            E_T = adata.var['E_T'].values
        return E_T

    def compute_et_new(counts, total_counts_per_gene):
        p_i = safe_divide_sparse(counts, total_counts_per_gene)
        log_p_i = safe_log2_sparse(p_i)
        entropy_per_gene = p_i.multiply(log_p_i).sum(axis=0).A1
        E_T = -entropy_per_gene
        E_T[np.isclose(total_counts_per_gene, 0)] = -1
        return E_T

    times_old, times_new = [], []
    for _ in range(reps):
        # Rebuild adata each time so old path always takes the "compute" branch
        fresh = make_adata(n_cells, n_genes)
        c = fresh.X if sparse.issparse(fresh.X) else csr_matrix(fresh.X)
        tcp = np.asarray(c.sum(axis=0)).ravel()

        t0 = time.perf_counter()
        compute_et_old(fresh, c, tcp)
        times_old.append(time.perf_counter() - t0)

        t0 = time.perf_counter()
        compute_et_new(c, tcp)
        times_new.append(time.perf_counter() - t0)

    return np.mean(times_old), np.mean(times_new)


# ══════════════════════════════════════════════════════════════════════════════
# Run & plot
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Running benchmarks… (this may take a minute)")

    print("  [1/4] Disk I/O round-trip vs in-memory…")
    io_old, io_new = bench_disk_vs_memory()

    print("  [2/4] Block loop: adata[mask].X vs counts[mask]…")
    bl_old, bl_new = bench_block_loop(n_cells=10000)

    print("  [3/4] q-value loop: row-by-row vs vectorized…")
    qv_old, qv_new = bench_qvalue_loop()

    print("  [4/4] E_T cache path overhead…")
    et_old, et_new = bench_et_cache(n_cells=10000)

    print("  [5/5] Full entropy math: sparse intermediate matrices vs COO…")
    from src.ember_py.generate_entropy_metrics import generate_entropy_metrics as entropy_new
    import anndata as ad
    _adata = make_adata(10000, 8000, 6)
    # old entropy: inline re-implementation using sparse intermediates
    from scipy.sparse import csr_matrix as _csr
    from scipy import sparse as _sparse
    def _safe_log2(m):
        m2 = m.copy(); m2.data = np.log2(m2.data); m2.data[~np.isfinite(m2.data)] = 0; return m2
    def _safe_div(num, den):
        d = np.copy(den); d[np.isclose(d,0)] = np.inf
        r = num.multiply(1.0/d); r.eliminate_zeros(); return r
    def entropy_old(adata, pl):
        c = adata.X
        if not _sparse.issparse(c): c = _csr(c)
        tot = np.asarray(c.sum(axis=0)).ravel()
        pi = _safe_div(c, tot); lpi = _safe_log2(pi)
        ET = -pi.multiply(lpi).sum(axis=0).A1; ET[np.isclose(tot,0)] = -1
        blks = adata.obs[pl].unique(); ng = adata.shape[1]; nb = len(blks)
        op = adata.obs[pl].values; EW = np.zeros(ng); PBN = np.zeros((ng,nb))
        for idx, blk in enumerate(blks):
            mask = op == blk; bc = c[mask,:]
            bs = np.asarray(bc.sum(axis=0)).ravel()
            q = _safe_div(bc,bs); lq = _safe_log2(q)
            ent = -q.multiply(lq).sum(axis=0).A1
            pcj = np.divide(bs,tot,out=np.zeros_like(bs),where=~np.isclose(tot,0))
            we = ent*pcj; EW+=we; PBN[:,idx]=we
        with np.errstate(invalid='ignore',divide='ignore'):
            Psi = np.where(ET>0, EW/ET, -1)
        with np.errstate(divide='ignore',invalid='ignore'):
            Pb = np.divide(PBN,EW[:,None],where=EW[:,None]!=0); Pb[~np.isfinite(Pb)]=0
        Pb_df = pd.DataFrame(Pb, index=adata.var.index, columns=blks)
        ez = -np.nansum(Pb_df*np.log2(Pb_df.where(Pb_df>0)),axis=1)
        return Psi, Pb_df, 1-(ez/np.log2(nb))
    REPS = 6
    t_e_old = np.mean([( t0:=time.perf_counter(), entropy_old(_adata,"partition"), time.perf_counter()-t0)[2] for _ in range(REPS)])
    t_e_new = np.mean([( t0:=time.perf_counter(), entropy_new(_adata,"partition"), time.perf_counter()-t0)[2] for _ in range(REPS)])

    labels   = ["Disk I/O\nround-trip", "Block loop\n(per draw)", "q-value\nassignment", "E_T cache\npath", "Entropy math\n(full, per draw)"]
    old_vals = [io_old, bl_old, qv_old, et_old, t_e_old]
    new_vals = [io_new, bl_new, qv_new, et_new, t_e_new]
    speedups = [o / n for o, n in zip(old_vals, new_vals)]

    # ── Overall light_ember estimate (100 draws, single-threaded) ─────────────
    # bench_disk_vs_memory used 60 draws; scale I/O cost to 100.
    # Per-draw compute uses the full entropy math timing (t_e_*).
    N_DRAWS = 100
    IO_SCALE = N_DRAWS / 60

    overall_segments = {
        "Per-draw compute":  (N_DRAWS * t_e_old,  N_DRAWS * t_e_new),
        "Disk I/O overhead": (io_old * IO_SCALE,  io_new * IO_SCALE),
        "q-value calc":      (qv_old,              qv_new),
    }
    total_old = sum(v[0] for v in overall_segments.values())
    total_new = sum(v[1] for v in overall_segments.values())
    overall_speedup = total_old / total_new

    # ── Figure: 3-panel layout ────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(19, 5))
    fig.patch.set_facecolor("#f8f8f8")

    OLD_COLOR = "#e07b54"
    NEW_COLOR = "#5b8db8"
    x = np.arange(len(labels))
    w = 0.35

    # Panel 1: absolute times per operation
    ax = axes[0]
    ax.set_facecolor("#f8f8f8")
    bars_old = ax.bar(x - w/2, old_vals, w, color=OLD_COLOR, zorder=3)
    bars_new = ax.bar(x + w/2, new_vals, w, color=NEW_COLOR, zorder=3)

    for bar, val in zip(bars_old, old_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(old_vals)*0.01,
                f"{val:.2f}s", ha="center", va="bottom", fontsize=8.5, color=OLD_COLOR, fontweight="bold")
    for bar, val in zip(bars_new, new_vals):
        lbl = f"{val*1000:.1f}ms" if val < 0.1 else f"{val:.2f}s"
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(old_vals)*0.01,
                lbl, ha="center", va="bottom", fontsize=8.5, color=NEW_COLOR, fontweight="bold")

    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Wall time (seconds)", fontsize=11)
    ax.set_title("Absolute runtime per operation", fontsize=12, fontweight="bold")
    ax.legend(handles=[
        mpatches.Patch(color=OLD_COLOR, label="Old"), mpatches.Patch(color=NEW_COLOR, label="New"),
    ], fontsize=10)
    ax.set_ylim(0, max(old_vals) * 1.25)
    ax.yaxis.grid(True, linestyle="--", alpha=0.5, zorder=0); ax.set_axisbelow(True)
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)

    # Panel 2: speedup per operation
    ax2 = axes[1]
    ax2.set_facecolor("#f8f8f8")
    bar_colors = [NEW_COLOR if s >= 1 else OLD_COLOR for s in speedups]
    bars_su = ax2.bar(x, speedups, 0.55, color=bar_colors, zorder=3)
    ax2.axhline(1.0, color="black", linewidth=1.2, linestyle="--", zorder=4, label="No change (1×)")

    for bar, su in zip(bars_su, speedups):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                 f"{su:.1f}×", ha="center", va="bottom", fontsize=11, fontweight="bold")

    ax2.set_xticks(x); ax2.set_xticklabels(labels, fontsize=10)
    ax2.set_ylabel("Speedup (old / new)", fontsize=11)
    ax2.set_title("Speedup per operation", fontsize=12, fontweight="bold")
    ax2.set_ylim(0, max(speedups) * 1.25)
    ax2.legend(fontsize=10)
    ax2.yaxis.grid(True, linestyle="--", alpha=0.5, zorder=0); ax2.set_axisbelow(True)
    for sp in ["top", "right"]: ax2.spines[sp].set_visible(False)

    # Panel 3: overall light_ember stacked bar
    ax3 = axes[2]
    ax3.set_facecolor("#f8f8f8")

    seg_colors = ["#7fbfad", OLD_COLOR, "#a78cc9"]
    bottoms = [0.0, 0.0]
    seg_names = list(overall_segments.keys())
    for seg_name, color in zip(seg_names, seg_colors):
        v_old, v_new = overall_segments[seg_name]
        ax3.bar([0], [v_old], 0.5, bottom=[bottoms[0]], color=color,
                label=seg_name, zorder=3, alpha=0.9)
        ax3.bar([1], [v_new], 0.5, bottom=[bottoms[1]], color=color,
                zorder=3, alpha=0.9)
        bottoms[0] += v_old
        bottoms[1] += v_new

    # Total time labels
    for xpos, total, color in [(0, total_old, OLD_COLOR), (1, total_new, NEW_COLOR)]:
        ax3.text(xpos, total + total_old * 0.02, f"{total:.1f}s",
                 ha="center", va="bottom", fontsize=12, fontweight="bold", color=color)

    ax3.annotate(
        f"Overall\n{overall_speedup:.1f}× faster",
        xy=(0.5, (total_old + total_new) / 2),
        fontsize=12, fontweight="bold", ha="center", va="center",
        color="#333333",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#aaaaaa", alpha=0.85),
    )

    ax3.set_xticks([0, 1])
    ax3.set_xticklabels(["Old", "New"], fontsize=12)
    ax3.set_ylabel("Estimated wall time (seconds)", fontsize=11)
    ax3.set_title(f"Full light_ember run\n({N_DRAWS} draws, single-threaded)", fontsize=12, fontweight="bold")
    ax3.legend(fontsize=9, loc="upper right")
    ax3.set_ylim(0, total_old * 1.25)
    ax3.yaxis.grid(True, linestyle="--", alpha=0.5, zorder=0); ax3.set_axisbelow(True)
    for sp in ["top", "right"]: ax3.spines[sp].set_visible(False)

    fig.suptitle(
        "ember speed comparison: old vs new implementation\n"
        "(synthetic data: 10 000 cells × 8 000 genes × 6 blocks)",
        fontsize=12, y=1.02
    )
    plt.tight_layout()
    out_path = os.path.join(os.path.dirname(__file__), "assets", "speed_comparison.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\nFigure saved to {out_path}")

    print("\nResults summary:")
    for label, o, n, s in zip(labels, old_vals, new_vals, speedups):
        print(f"  {label.replace(chr(10),' '):35s}  old={o:.3f}s  new={n:.4f}s  speedup={s:.2f}x")
    print(f"\n  Overall light_ember ({N_DRAWS} draws):  old≈{total_old:.1f}s  new≈{total_new:.1f}s  speedup={overall_speedup:.2f}x")
