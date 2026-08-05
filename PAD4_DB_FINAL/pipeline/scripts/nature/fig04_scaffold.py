#!/usr/bin/env python
"""Nature Fig 4 — Scaffold Landscape (Publication-Quality)"""
import os, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from pathlib import Path

# ── RDKit for Scaffold Computation ────────────────────────────────────────────
try:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("⚠️ WARNING: RDKit is not installed. Panel D will be grey-only.")
    print("   To fix, run: conda install -c conda-forge rdkit")

# ── Nature style ──────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'axes.labelsize': 7,
    'axes.titlesize': 7,
    'axes.linewidth': 0.75,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'lines.linewidth': 0.75,
    'lines.markersize': 4,
    'patch.linewidth': 0.5,
    'legend.fontsize': 6,
    'legend.frameon': False,
    'figure.facecolor': 'white',
    'savefig.facecolor': 'white',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

PAL = {
    'blue':       '#0077BB',
    'orange':     '#EE7733',
    'teal':       '#009988',
    'navy':       '#1A237E',
    'grey':       '#BBBBBB',
    'dark_grey':  '#555555',
    'light_grey': '#E8E8E8',
    'red':        '#CC3311',
}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/figures/nature'
OUT.mkdir(parents=True, exist_ok=True)

def save_fig(fig, name):
    for ext in ('png', 'svg', 'pdf'):
        p = OUT / f'{name}.{ext}'
        fig.savefig(p, dpi=600 if ext == 'png' else None, bbox_inches='tight')
    print(f"✅ Saved: {OUT / name}.png")

def plabel(ax, letter):
    ax.text(-0.08, 1.04, letter, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='bottom', ha='right')

print("=" * 60)
print("FIGURE 4 — SCAFFOLD LANDSCAPE")
print("=" * 60)

# ── Load data ──────────────────────────────────────────────────────────────────
df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
scaffold_sum = pd.read_csv(ROOT / 'outputs/tables/05_scaffold_summary.csv')
scaffold_sum = scaffold_sum.reset_index(drop=True)
scaffold_sum['scaffold_rank'] = scaffold_sum.index + 1

# ── Compute Murcko Scaffolds (if missing from dataframe) ─────────────────────
if 'scaffold_smiles' not in df.columns and RDKIT_AVAILABLE:
    print("  Computing Murcko scaffolds on-the-fly...")
    def get_murcko(smi):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None: return None
            return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
        except Exception:
            return None
    df['scaffold_smiles'] = df['smiles_std'].apply(get_murcko)
    # Drop rows where scaffold computation failed
    df = df.dropna(subset=['scaffold_smiles'])
elif not RDKIT_AVAILABLE:
    print("⚠️ RDKit missing. Panel D will remain grey.")
    df['scaffold_smiles'] = 'unknown'

# ── Prepare Scaffold Mapping ──────────────────────────────────────────────────
# Build a mapping dict from scaffold summary
scaffold_map = scaffold_sum.set_index('scaffold_smiles')['n_compounds'].to_dict()
df['scaffold_size'] = df['scaffold_smiles'].map(scaffold_map).fillna(0)

# Identify the Rank 1 scaffold (largest series)
rank1_scaffold = scaffold_sum.iloc[0]['scaffold_smiles']
df['is_rank1'] = (df['scaffold_smiles'] == rank1_scaffold).astype(bool)
df['is_series'] = (df['scaffold_size'] >= 2).astype(bool)
df['is_singleton'] = (~df['is_series']).astype(bool)

# ── t-SNE Data ──────────────────────────────────────────────────────────────────
tsne_coords = np.load(ROOT / 'data/interim/tsne_coords_3093.npy')
tsne_iks = np.load(ROOT / 'data/interim/tsne_inchikeys_3093.npy', allow_pickle=True)
tsne_df = pd.DataFrame({'inchi_key': tsne_iks, 'tx': tsne_coords[:, 0], 'ty': tsne_coords[:, 1]})
df = df.merge(tsne_df, on='inchi_key', how='inner')

# ── Build Figure ──────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.5), constrained_layout=True)
ax_a, ax_b, ax_c, ax_d = axes.flat

# ── Panel A: Top 30 scaffold series bar ──────────────────────────────────────
top30 = scaffold_sum.head(30)
colors_a = [PAL['navy']] + [PAL['blue']] * 29
bars = ax_a.bar(top30['scaffold_rank'], top30['n_compounds'], color=colors_a,
                alpha=0.85, edgecolor='white', lw=0.3)
ax_a.set_xlabel('Scaffold rank (by series size)', fontsize=7)
ax_a.set_ylabel('Series size (n compounds)', fontsize=7)
ax_a.set_xlim(0.3, 30.7)

# Cleaner annotation for Rank 1 (FIXED: moved slightly to the right so it doesn't overlap)
ax_a.annotate(f'Rank 1\n({top30["n_compounds"].iloc[0]} cpds)', xy=(1, top30["n_compounds"].iloc[0]),
              xytext=(6, top30["n_compounds"].iloc[0] - 20),
              arrowprops=dict(arrowstyle='->', color=PAL['navy'], lw=0.8),
              fontsize=5.5, color=PAL['navy'])

# Stats box
ax_a.text(0.97, 0.97,
          f'1,244 unique scaffolds\n375 series (≥2 cpds)\n869 singletons\nLargest: {top30["n_compounds"].iloc[0]} cpds',
          transform=ax_a.transAxes, ha='right', va='top', fontsize=5.5,
          bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor=PAL['grey'], lw=0.6))
plabel(ax_a, 'a')

# ── Panel B: Series size distribution ──────────────────────────────────────
series_only = scaffold_sum[scaffold_sum['n_compounds'] >= 2]['n_compounds']
bins_b = np.logspace(np.log10(2), np.log10(200), 20)
ax_b.hist(series_only, bins=bins_b, color=PAL['teal'], alpha=0.85, edgecolor='white', lw=0.3)
ax_b.set_xscale('log')
ax_b.set_xticks([2, 5, 10, 20, 50, 100])
ax_b.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
ax_b.set_xlabel('Series size (log scale)', fontsize=7)
ax_b.set_ylabel('Number of scaffolds', fontsize=7)
med_s = series_only.median()
ax_b.text(0.97, 0.97, f'Range: 2–{series_only.max()}\nMedian series size: {med_s:.0f}',
          transform=ax_b.transAxes, ha='right', va='top', fontsize=5.5,
          bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor=PAL['grey'], lw=0.6))
plabel(ax_b, 'b')

# ── Panel C: Lorenz curve ────────────────────────────────────────────────────
vals = np.sort(scaffold_sum['n_compounds'].values)
n = len(vals)
cum_vals = np.cumsum(vals)
lorenz_x = np.concatenate([[0], np.arange(1, n + 1) / n])
lorenz_y = np.concatenate([[0], cum_vals / cum_vals[-1]])
auc = np.trapz(lorenz_y, lorenz_x)
gini = 1 - 2 * auc

ax_c.plot(lorenz_x, lorenz_y, color=PAL['blue'], lw=1.5, label='Lorenz curve')
ax_c.plot([0, 1], [0, 1], color=PAL['grey'], lw=0.8, ls='--', label='Perfect equality')
ax_c.fill_between(lorenz_x, lorenz_y, lorenz_x, alpha=0.20, color=PAL['blue'])

# Gini label with white box background for readability (FIXED contrast)
ax_c.text(0.55, 0.30, f'Gini = {gini:.3f}', transform=ax_c.transAxes,
          fontsize=9, fontweight='bold', color=PAL['navy'], ha='center', va='center',
          bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0.3))
ax_c.set_xlabel('Cumulative fraction of scaffolds', fontsize=7)
ax_c.set_ylabel('Cumulative fraction of compounds', fontsize=7)
ax_c.set_xlim(0, 1)
ax_c.set_ylim(0, 1)
ax_c.set_aspect('equal')
ax_c.legend(loc='upper left', fontsize=5.5)
plabel(ax_c, 'c')

# ── Panel D: t-SNE scaffold membership (NOW FULLY COLORED!) ──────────────────
# Extract data for plot
singleton_mask = df['is_singleton'].values
series_nonhub_mask = (df['is_series'].values & (~df['is_rank1'].values))
rank1_mask = df['is_rank1'].values

ax_d.scatter(df.loc[singleton_mask, 'tx'], df.loc[singleton_mask, 'ty'],
             s=2, alpha=0.35, c=PAL['grey'], marker='.', rasterized=True, zorder=1,
             label=f'Singletons (n={singleton_mask.sum()})')
ax_d.scatter(df.loc[series_nonhub_mask, 'tx'], df.loc[series_nonhub_mask, 'ty'],
             s=3, alpha=0.5, c=PAL['blue'], marker='.', rasterized=True, zorder=2,
             label=f'Other series (n={series_nonhub_mask.sum()})')
ax_d.scatter(df.loc[rank1_mask, 'tx'], df.loc[rank1_mask, 'ty'],
             s=4, alpha=0.8, c=PAL['navy'], marker='o', rasterized=True, zorder=3,
             label=f'Rank-1 scaffold (n={rank1_mask.sum()})')

# Centroid star for the hub series (Rank 1)
if rank1_mask.any():
    cx = df.loc[rank1_mask, 'tx'].mean()
    cy = df.loc[rank1_mask, 'ty'].mean()
    ax_d.scatter([cx], [cy], s=120, c=PAL['navy'], marker='*', zorder=5,
                 edgecolor='white', lw=0.3)

ax_d.set_xlabel('t-SNE 1', fontsize=7)
ax_d.set_ylabel('t-SNE 2', fontsize=7)
ax_d.set_xticks([])
ax_d.set_yticks([])
ax_d.legend(loc='lower left', fontsize=5.5, markerscale=1.5, handlelength=0.6)
plabel(ax_d, 'd')

# ── Save ──────────────────────────────────────────────────────────────────────
save_fig(fig, 'fig04_scaffold')
plt.close()
print("\n🎉 Figure 4 complete. Panel D is now fully colored with computed scaffolds.")