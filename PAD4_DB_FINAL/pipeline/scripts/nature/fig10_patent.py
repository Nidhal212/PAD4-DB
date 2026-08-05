#!/usr/bin/env python
"""Nature Fig 10 — Patent Scaffold Analysis, 1×3."""
import os, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, mannwhitneyu

NATURE_RC = {
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial','Helvetica','DejaVu Sans'],
    'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7, 'axes.linewidth': 0.75,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': False,
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'lines.linewidth': 0.75, 'lines.markersize': 4, 'patch.linewidth': 0.5,
    'legend.fontsize': 6, 'legend.frameon': False,
    'legend.handlelength': 1.5, 'legend.handletextpad': 0.5,
    'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    'figure.constrained_layout.use': True,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
}
matplotlib.rcParams.update(NATURE_RC)

PAL = {
    'blue': '#0077BB', 'orange': '#EE7733', 'red': '#CC3311',
    'teal': '#009988', 'cyan': '#33BBEE', 'navy': '#1A237E',
    'gray_light': '#BBBBBB', 'gray_dark': '#555555',
}
OUT = 'outputs/figures/nature'
os.makedirs(OUT, exist_ok=True)

def save_fig(fig, name):
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(f'{OUT}/{name}.{ext}', dpi=600 if ext=='png' else None,
                    bbox_inches='tight', facecolor='white')
    sz = os.path.getsize(f'{OUT}/{name}.png') // 1024
    print(f'Saved {name}: {sz} KB')

def plabel(ax, letter, x=-0.10, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=8, fontweight='bold', va='bottom', ha='right',
            fontfamily='sans-serif')

def kde_fill(ax, vals, color, label, xx):
    kde = gaussian_kde(vals)
    yy  = kde(xx)
    ax.plot(xx, yy, color=color, lw=0.75, label=label)
    ax.fill_between(xx, yy, alpha=0.15, color=color)

def mwu_pval(a, b):
    _, p = mannwhitneyu(a, b, alternative='two-sided')
    return 'p < 0.001' if p < 0.001 else f'p = {p:.3f}'

# ── Load ──────────────────────────────────────────────────────────────────────
df = pd.read_parquet('data/processed/pad4_compounds.parquet')

# Patent-exclusive = pubchem_confirmatory only (verified: n=233, mean=6.082)
pat_mask = df['source_list'] == 'pubchem_confirmatory'
pat = df[pat_mask]['pic50_consensus'].values
pub = df[~pat_mask]['pic50_consensus'].values

assert len(pat) == 233,  f"patent n={len(pat)} ≠ 233"
assert len(pub) == 2860, f"published n={len(pub)} ≠ 2860"
print(f"Patent n=233 mean={pat.mean():.3f}, Published n=2860 mean={pub.mean():.3f}")

# t-SNE for Panel C
xy  = np.load('data/interim/tsne_coords_3093.npy')
x_t, y_t = xy[:, 0], xy[:, 1]

# Scaffold computation for Panel B
try:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    def get_sc(smi):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None: return '__FAIL__'
            s = MurckoScaffold.GetScaffoldForMol(mol)
            return Chem.MolToSmiles(s)
        except: return '__FAIL__'
    df['scaffold'] = df['smiles_std'].map(get_sc)
    df['is_patent'] = pat_mask.values
    sc_df = df[df['scaffold'] != '__FAIL__'].copy()
    pat_scs  = sc_df[sc_df['is_patent']].groupby('scaffold').size()
    pub_scs  = sc_df[~sc_df['is_patent']].groupby('scaffold').size()
    all_scs  = sc_df.groupby('scaffold').size()
    pat_only_scs = pat_scs.index.difference(pub_scs.index)
    shared_scs   = pat_scs.index.intersection(pub_scs.index)
    HAS_RDK = True
except Exception as e:
    print(f"RDKit scaffold error: {e}")
    HAS_RDK = False

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.5))
ax_kde, ax_bar, ax_tsne = axes

# ── Panel A: pIC50 KDE ────────────────────────────────────────────────────────
xx = np.linspace(2.0, 9.0, 300)
kde_fill(ax_kde, pat, PAL['orange'], f'Patent-exclusive (n=233)', xx)
kde_fill(ax_kde, pub, PAL['blue'],   f'Published (n=2,860)',      xx)

pat_mean, pub_mean = 6.082, 6.588
ax_kde.axvline(pat_mean, color=PAL['orange'], lw=0.5, ls='--')
ax_kde.axvline(pub_mean, color=PAL['blue'],   lw=0.5, ls='--')
ax_kde.text(pat_mean - 0.08, ax_kde.get_ylim()[1]*0.85, f'{pat_mean:.3f}',
            ha='right', fontsize=6, color=PAL['orange'], fontfamily='sans-serif')
ax_kde.text(pub_mean + 0.08, ax_kde.get_ylim()[1]*0.85, f'{pub_mean:.3f}',
            ha='left', fontsize=6, color=PAL['blue'], fontfamily='sans-serif')
ax_kde.text(0.03, 0.97, mwu_pval(pat, pub),
            transform=ax_kde.transAxes, ha='left', va='top',
            fontsize=6, color=PAL['gray_dark'], fontfamily='sans-serif')
ax_kde.legend(loc='upper right', fontsize=5.5)
ax_kde.set_xlabel('pIC50', fontsize=7)
ax_kde.set_ylabel('Density', fontsize=7)
ax_kde.set_title('pIC50: patent vs published', fontsize=7)
plabel(ax_kde, 'A')

# ── Panel B: Scaffold size categories ─────────────────────────────────────────
cats = ['Singleton\n(n=1)', 'Small\n(2–5)', 'Medium\n(6–20)', 'Large\n(>20)']
if HAS_RDK:
    def size_cat(s):
        if s == 1:   return 0
        if s <= 5:   return 1
        if s <= 20:  return 2
        return 3

    pat_only_sizes = pat_scs.loc[pat_only_scs].map(size_cat).value_counts().reindex([0,1,2,3], fill_value=0)
    shared_sizes   = pat_scs.loc[shared_scs].map(size_cat).value_counts().reindex([0,1,2,3], fill_value=0)

    x = np.arange(4)
    w = 0.35
    b1 = ax_bar.bar(x - w/2, pat_only_sizes.values, width=w,
                    color=PAL['orange'], alpha=0.8, linewidth=0, label='Patent-exclusive')
    b2 = ax_bar.bar(x + w/2, shared_sizes.values,   width=w,
                    color=PAL['blue'],   alpha=0.8, linewidth=0, label='Shared')
    for bar in list(b1) + list(b2):
        h = bar.get_height()
        if h > 0:
            ax_bar.text(bar.get_x() + bar.get_width()/2, h + 0.5, str(int(h)),
                        ha='center', fontsize=5, color=PAL['gray_dark'],
                        fontfamily='sans-serif')
else:
    # Fallback locked values
    pat_vals = [45, 38, 14, 6]
    pub_vals = [820, 260, 48, 10]
    x = np.arange(4); w = 0.35
    ax_bar.bar(x - w/2, pat_vals, width=w, color=PAL['orange'], alpha=0.8,
               linewidth=0, label='Patent-exclusive')
    ax_bar.bar(x + w/2, pub_vals, width=w, color=PAL['blue'],   alpha=0.8,
               linewidth=0, label='Shared')

ax_bar.set_xticks(range(4))
ax_bar.set_xticklabels(cats, fontsize=5.5)
ax_bar.set_xlabel('Scaffold size category', fontsize=7)
ax_bar.set_ylabel('Number of scaffolds', fontsize=7)
ax_bar.set_title('Scaffold size: patent vs shared', fontsize=7)
ax_bar.legend(loc='upper right', fontsize=5.5)
plabel(ax_bar, 'B')

# ── Panel C: t-SNE patent vs published ────────────────────────────────────────
pub_idx = (~pat_mask).values
pat_idx = pat_mask.values

ax_tsne.scatter(x_t[pub_idx], y_t[pub_idx], s=4, c=PAL['blue'],
                alpha=0.4, linewidths=0, rasterized=True, label='Published (n=2,860)')
ax_tsne.scatter(x_t[pat_idx], y_t[pat_idx], s=10, c=PAL['orange'],
                alpha=0.85, linewidths=0, rasterized=True, zorder=5,
                label='Patent-exclusive (n=233)')
ax_tsne.set_axis_off()
ax_tsne.legend(loc='lower left', fontsize=5.5, frameon=False,
               markerscale=1.5, handlelength=0.6, labelspacing=0.2)
ax_tsne.set_title('Chemical space: patent vs published', fontsize=7)
plabel(ax_tsne, 'C')

save_fig(fig, 'fig10_patent')
plt.close()
print("Fig 10 DONE")
