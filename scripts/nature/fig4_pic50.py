#!/usr/bin/env python
"""Nature Fig 4 — pIC50 Distribution (4 panels)."""
import os, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, gaussian_kde
from scipy.stats import gaussian_kde

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
        p = f'{OUT}/{name}.{ext}'
        fig.savefig(p, dpi=600 if ext == 'png' else None,
                    bbox_inches='tight', facecolor='white')
    sz = os.path.getsize(f'{OUT}/{name}.png') / 1024
    print(f"Saved {name}: {sz:.0f} KB")

def plabel(ax, letter, x=-0.08, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=8, fontweight='bold', va='bottom', ha='right',
            fontfamily='sans-serif')

def mwu_label(a, b):
    _, p = mannwhitneyu(a, b, alternative='two-sided')
    if p < 0.001: return 'p < 0.001'
    if p < 0.01:  return f'p = {p:.3f}'
    return f'p = {p:.2f}'

def kde_curve(ax, vals, color, label, bw=None, xx=None):
    if xx is None:
        xx = np.linspace(vals.min() - 0.2, vals.max() + 0.2, 300)
    kde = gaussian_kde(vals, bw_method=bw)
    ax.plot(xx, kde(xx), color=color, lw=0.9, label=label)
    ax.fill_between(xx, kde(xx), alpha=0.12, color=color)

# ── Load ──────────────────────────────────────────────────────────────────────
df = pd.read_parquet('data/processed/pad4_compounds.parquet')
pic = df['pic50_consensus'].values.astype(float)

pubchem_only = df[df['source_list'] == 'pubchem_confirmatory']['pic50_consensus'].values.astype(float)
multi_source = df[df['n_sources'] >= 2]['pic50_consensus'].values.astype(float)

mech_colors = {
    'enzymatic':           PAL['gray_light'],
    'enzymatic_confirmed': PAL['teal'],
    'fp_ic50':             PAL['cyan'],
    'covalent':            PAL['red'],
}
mech_labels = {
    'enzymatic':           'Enzymatic',
    'enzymatic_confirmed': 'Enzymatic confirmed',
    'fp_ic50':             'FP-based IC50',
    'covalent':            'Covalent',
}

# Cliff compounds for degree scatter
cliffs = pd.read_parquet('data/processed/activity_cliffs.parquet')
severe = cliffs[cliffs['cliff_tier'] == 'severe']
from collections import Counter
deg_counter = Counter(severe['inchi_key_a'].tolist() + severe['inchi_key_b'].tolist())
deg_df = pd.DataFrame({'inchi_key': list(deg_counter.keys()),
                        'degree': list(deg_counter.values())})
deg_df = deg_df.merge(df[['inchi_key','pic50_consensus']], on='inchi_key', how='left')

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.5))
ax_h, ax_src, ax_vio, ax_deg = axes.flat

# ── Panel A: Global histogram + KDE ──────────────────────────────────────────
bins = np.arange(pic.min() - 0.25, pic.max() + 0.5, 0.25)
ax_h.hist(pic, bins=bins, color=PAL['blue'], alpha=0.4, linewidth=0, density=True)

ax_tw = ax_h.twinx()
ax_tw.spines['top'].set_visible(False)
ax_tw.spines['right'].set_linewidth(0.5)
xx = np.linspace(pic.min() - 0.5, pic.max() + 0.5, 300)
kde = gaussian_kde(pic)
ax_tw.plot(xx, kde(xx), color=PAL['navy'], lw=0.9, label='KDE')
ax_tw.set_ylabel('Density', fontsize=7)
ax_tw.tick_params(labelsize=6, width=0.4)

# Annotations
ax_h.axvline(np.median(pic), color=PAL['red'], lw=0.75, ls='--', label='Median')
ax_h.set_xlabel('pIC50', fontsize=7)
ax_h.set_ylabel('Count density', fontsize=7)
ax_h.set_title('pIC50 distribution (N=3,093)', fontsize=7)
ax_h.text(0.97, 0.96, f'Median={np.median(pic):.2f}\nMean={pic.mean():.2f}\nSD={pic.std():.2f}',
          transform=ax_h.transAxes, ha='right', va='top',
          fontsize=5.5, color=PAL['gray_dark'], fontfamily='sans-serif')
plabel(ax_h, 'A')

# ── Panel B: Source KDE comparison ────────────────────────────────────────────
xx2 = np.linspace(2.0, 9.0, 300)
kde_curve(ax_src, pubchem_only, PAL['orange'], f'PubChem-only (n={len(pubchem_only)})', xx=xx2)
kde_curve(ax_src, multi_source, PAL['blue'],   f'Multi-source (n={len(multi_source)})', xx=xx2)
pval = mwu_label(pubchem_only, multi_source)
ax_src.text(0.03, 0.97, f'Mann–Whitney U\n{pval}',
            transform=ax_src.transAxes, ha='left', va='top',
            fontsize=5.5, color=PAL['gray_dark'], fontfamily='sans-serif')
ax_src.legend(loc='upper right', fontsize=5.5, frameon=False)
ax_src.set_xlabel('pIC50', fontsize=7)
ax_src.set_ylabel('Density', fontsize=7)
ax_src.set_title('pIC50 by source coverage', fontsize=7)
plabel(ax_src, 'B')

# ── Panel C: Violin by mechanism ──────────────────────────────────────────────
mechs = ['enzymatic', 'enzymatic_confirmed', 'fp_ic50', 'covalent']
mech_data = [df[df['mechanism_class'] == m]['pic50_consensus'].values for m in mechs]
mech_ns   = [len(d) for d in mech_data]

vp = ax_vio.violinplot(mech_data, positions=range(len(mechs)),
                        widths=0.6, showmedians=False, showextrema=False)
for body, mech in zip(vp['bodies'], mechs):
    body.set_facecolor(mech_colors[mech])
    body.set_edgecolor('none')
    body.set_alpha(0.65)

# Strip only for small groups; median line for all
for i, (d, mech) in enumerate(zip(mech_data, mechs)):
    if mech != 'enzymatic':  # skip strip for n=2079 (too dense, occludes violin)
        jit = np.random.default_rng(42).uniform(-0.10, 0.10, size=len(d))
        # Use darker shade so points are distinguishable against violin fill
        strip_clr = '#006666' if mech == 'fp_ic50' else mech_colors[mech]
        ax_vio.scatter(i + jit, d, s=1.2, color=strip_clr,
                       alpha=0.8, linewidths=0, zorder=3)
    med = np.median(d)
    ax_vio.hlines(med, i - 0.25, i + 0.25, colors='black', linewidths=1.0, zorder=4)

ax_vio.set_xticks(range(len(mechs)))
ax_vio.set_xticklabels([f'{mech_labels[m]}\n(n={n:,})' for m, n in zip(mechs, mech_ns)],
                        fontsize=5.5, multialignment='center')
ax_vio.set_xlabel('Assay mechanism', fontsize=7)
ax_vio.set_ylabel('pIC50', fontsize=7)
ax_vio.set_title('pIC50 by assay mechanism', fontsize=7)
plabel(ax_vio, 'C')

# ── Panel D: Degree vs pIC50 — hub compounds colored ─────────────────────────
HUB_A = frozenset({'SMADULGDNOCLOP-GISFHXKWSA-N', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N'})
HUB_B = frozenset({'UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'DVCKJOQIVOGXEI-XMMPIXPASA-N'})
mask_a = deg_df['inchi_key'].isin(HUB_A)
mask_b = deg_df['inchi_key'].isin(HUB_B)
mask_n = ~(mask_a | mask_b)

ax_deg.scatter(deg_df.loc[mask_n, 'pic50_consensus'], deg_df.loc[mask_n, 'degree'],
               s=10, c=PAL['blue'], alpha=0.65, linewidths=0.3, edgecolors='white', zorder=3)
ax_deg.scatter(deg_df.loc[mask_a, 'pic50_consensus'], deg_df.loc[mask_a, 'degree'],
               s=40, c=PAL['navy'], marker='*', linewidths=0, zorder=5,
               label='Hub Class A')
ax_deg.scatter(deg_df.loc[mask_b, 'pic50_consensus'], deg_df.loc[mask_b, 'degree'],
               s=40, c=PAL['red'],  marker='*', linewidths=0, zorder=5,
               label='Hub Class B')
ax_deg.legend(loc='upper right', fontsize=5.5, frameon=False)

# Annotate top-degree compounds
top4 = deg_df.nlargest(4, 'degree')
for _, row in top4.iterrows():
    ax_deg.annotate(f"{int(row['degree'])}", (row['pic50_consensus'], row['degree']),
                    xytext=(3, 2), textcoords='offset points',
                    fontsize=5.0, color=PAL['navy'], fontfamily='sans-serif')

ax_deg.set_xlabel('pIC50', fontsize=7)
ax_deg.set_ylabel('Severe cliff degree', fontsize=7)
ax_deg.set_title('Cliff degree vs pIC50 (n=99 severe cliff compounds)', fontsize=7)
ax_deg.text(0.02, 0.02, 'Numbers = cliff degree (severe pairs)',
            transform=ax_deg.transAxes, fontsize=5.0,
            color=PAL['gray_dark'], va='bottom', fontfamily='sans-serif')
plabel(ax_deg, 'D')

save_fig(fig, 'fig4_pic50')
plt.close()
print("Fig 4 DONE")
