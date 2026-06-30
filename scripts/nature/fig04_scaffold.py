"""
fig04_scaffold.py — Figure 4: Scaffold Landscape (4-panel)
Outputs: outputs/figures/nature/fig04_scaffold.{png,pdf}
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Liberation Sans', 'DejaVu Sans', 'Arial'],
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 3,    'ytick.major.size': 3,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,     'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

COLORS = {
    'blue':       '#0077BB',
    'orange':     '#EE7733',
    'teal':       '#009988',
    'navy':       '#004488',
    'grey':       '#BBBBBB',
    'dark_grey':  '#555555',
    'light_grey': '#E8E8E8',
    'red':        '#CC3311',
}

CANON = {
    'n_compounds': 3093,
    'n_scaffolds': 1244,
    'n_series': 375,
    'n_singletons': 869,
    'largest_series': 174,
    'n_in_series': 2224,
}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/figures/nature_v2'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FIGURE 4 — SCAFFOLD LANDSCAPE")
print("=" * 60)

df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
scaffold_sum = pd.read_csv(ROOT / 'outputs/tables/05_scaffold_summary.csv')
scaffold_sum = scaffold_sum.reset_index()
scaffold_sum['scaffold_rank'] = scaffold_sum.index + 1
scaffold_sum = scaffold_sum.sort_values('n_compounds', ascending=False).reset_index(drop=True)
scaffold_sum['scaffold_rank'] = scaffold_sum.index + 1

assert len(scaffold_sum) == CANON['n_scaffolds']
print(f"  Scaffolds: {len(scaffold_sum)} ✓")
print(f"  Largest: {scaffold_sum['n_compounds'].max()} ✓")

fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
ax_a, ax_b, ax_c, ax_d = axes.flatten()

# ── Panel a: Top 30 scaffold bar chart ───────────────────────────────────────
print("\n[Panel a] Top 30 scaffold bar ...")
top30 = scaffold_sum.head(30)
colors_a = [COLORS['navy']] + [COLORS['blue']] * 29
ax_a.bar(top30['scaffold_rank'], top30['n_compounds'],
          color=colors_a, alpha=0.85, edgecolor='white', lw=0.3)
ax_a.set_xlabel('Scaffold rank (by series size)', fontsize=9)
ax_a.set_ylabel('Series size (n compounds)', fontsize=9)
ax_a.set_xlim(0.3, 30.7)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)

# Stats box
stats_txt = (f'1,244 unique scaffolds\n375 series (≥2 cpds)\n'
             f'869 singletons\nLargest: 174 cpds')
ax_a.text(0.97, 0.97, stats_txt, transform=ax_a.transAxes, fontsize=7.5,
          ha='right', va='top',
          bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                    edgecolor=COLORS['grey'], lw=0.8))

# Arrow to rank 1
ax_a.annotate('Rank 1\n(174 cpds)', xy=(1, top30.iloc[0]['n_compounds']),
               xytext=(5, top30.iloc[0]['n_compounds'] - 15),
               arrowprops=dict(arrowstyle='->', color=COLORS['navy'], lw=1.0),
               fontsize=7, color=COLORS['navy'])

ax_a.text(-0.12, 1.04, 'a', transform=ax_a.transAxes, fontsize=11, fontweight='bold')

# ── Panel b: Series size distribution ────────────────────────────────────────
print("[Panel b] Series size distribution ...")
series_only = scaffold_sum[scaffold_sum['n_compounds'] >= 2]['n_compounds']
print(f"  Series: {len(series_only)} | range: {series_only.min()}–{series_only.max()}")

bins_b = np.logspace(np.log10(2), np.log10(200), 20)
ax_b.hist(series_only, bins=bins_b, color=COLORS['teal'], alpha=0.8, edgecolor='white', lw=0.3)
ax_b.set_xscale('log')
ax_b.set_xticks([2, 5, 10, 20, 50, 100])
ax_b.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax_b.set_xlabel('Series size (log scale)', fontsize=9)
ax_b.set_ylabel('Number of scaffolds', fontsize=9)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

med_s = series_only.median()
ax_b.text(0.97, 0.97, f'Range: 2–174\nMedian series size: {med_s:.0f}',
          transform=ax_b.transAxes, fontsize=7.5, ha='right', va='top',
          bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                    edgecolor=COLORS['grey'], lw=0.8))
ax_b.text(-0.13, 1.04, 'b', transform=ax_b.transAxes, fontsize=11, fontweight='bold')

# ── Panel c: Lorenz curve ────────────────────────────────────────────────────
print("[Panel c] Lorenz curve ...")
vals = np.sort(scaffold_sum['n_compounds'].values)
n = len(vals)
cum_vals = np.cumsum(vals)
lorenz_x = np.concatenate([[0], np.arange(1, n + 1) / n])
lorenz_y = np.concatenate([[0], cum_vals / cum_vals[-1]])

# Gini coefficient: G = 1 - 2 * AUC
auc = np.trapz(lorenz_y, lorenz_x)
gini = 1 - 2 * auc
print(f"  Gini coefficient: {gini:.3f}")

ax_c.plot(lorenz_x, lorenz_y, color=COLORS['blue'], lw=1.5, label='Lorenz curve')
ax_c.plot([0, 1], [0, 1], color=COLORS['grey'], lw=1.0, ls='--', label='Perfect equality')
ax_c.fill_between(lorenz_x, lorenz_y, lorenz_x, alpha=0.15, color=COLORS['blue'])

ax_c.text(0.55, 0.25, f'Gini = {gini:.3f}', transform=ax_c.transAxes,
          fontsize=11, fontweight='bold', color=COLORS['navy'],
          ha='center', va='center')

ax_c.set_xlabel('Cumulative fraction of scaffolds', fontsize=9)
ax_c.set_ylabel('Cumulative fraction of compounds', fontsize=9)
ax_c.set_xlim(0, 1)
ax_c.set_ylim(0, 1)
ax_c.set_aspect('equal')
ax_c.legend(fontsize=7.5, loc='upper left', framealpha=0.7, edgecolor='none')
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
ax_c.text(-0.13, 1.04, 'c', transform=ax_c.transAxes, fontsize=11, fontweight='bold')

# ── Panel d: t-SNE scaffold membership ──────────────────────────────────────
print("[Panel d] t-SNE scaffold membership ...")
tsne_coords = np.load(ROOT / 'data/interim/tsne_coords_3093.npy')
tsne_iks = np.load(ROOT / 'data/interim/tsne_inchikeys_3093.npy', allow_pickle=True)
tsne_df = pd.DataFrame({'inchi_key': tsne_iks, 'tx': tsne_coords[:, 0], 'ty': tsne_coords[:, 1]})
df_tsne = df.merge(tsne_df, on='inchi_key', how='inner')

rank1 = df_tsne[df_tsne['scaffold_rank'] == 1]
rank1_iks = set(rank1['inchi_key'])
series_mem = df_tsne[(df_tsne['scaffold_series_size'] >= 2) & (~df_tsne['inchi_key'].isin(rank1_iks))]
# Compounds with NaN scaffold_series_size (no Murcko scaffold) counted as singletons
singletons = df_tsne[~df_tsne['inchi_key'].isin(rank1_iks | set(series_mem['inchi_key']))]

print(f"  Singletons: {len(singletons)} | Series: {len(series_mem)} | Rank-1: {len(rank1)}")

ax_d.scatter(singletons['tx'], singletons['ty'], s=3, alpha=0.3, c=COLORS['grey'],
              marker='.', rasterized=True, zorder=1)
ax_d.scatter(series_mem['tx'], series_mem['ty'], s=4, alpha=0.5, c=COLORS['blue'],
              marker='.', rasterized=True, zorder=2)
ax_d.scatter(rank1['tx'], rank1['ty'], s=6, alpha=0.8, c=COLORS['navy'],
              marker='o', rasterized=True, zorder=3)

# Star at centroid of rank-1 cluster
if len(rank1) > 0:
    cx, cy = rank1['tx'].mean(), rank1['ty'].mean()
    ax_d.scatter([cx], [cy], s=200, c=COLORS['navy'], marker='*', zorder=5)

ax_d.set_xlabel('t-SNE 1', fontsize=9)
ax_d.set_ylabel('t-SNE 2', fontsize=9)
ax_d.set_xticks([])
ax_d.set_yticks([])
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)

leg_handles = [
    mpatches.Patch(color=COLORS['grey'], label=f'Singletons (n=872)'),
    mpatches.Patch(color=COLORS['blue'], label=f'Series members (n={len(series_mem)})'),
    mpatches.Patch(color=COLORS['navy'], label=f'Rank-1 scaffold (n={len(rank1)})'),
]
ax_d.legend(handles=leg_handles, fontsize=7, loc='lower left', framealpha=0.7, edgecolor='none')
ax_d.text(-0.12, 1.04, 'd', transform=ax_d.transAxes, fontsize=11, fontweight='bold')

# ── Save ──────────────────────────────────────────────────────────────────────
import matplotlib.ticker
for ext in ['png', 'pdf']:
    outpath = OUT / f'fig04_scaffold.{ext}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"  Saved: {outpath}")

plt.close(fig)
print("\nFigure 4 complete.")
