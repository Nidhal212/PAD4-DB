#!/usr/bin/env python
"""Nature Fig 11 — Source Independence, 1×3."""
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

# Use threshold=0.6 (locked: True=528, False=2565)
TRUE_THRESH = 0.6
true_mask = df['source_independence_score'] >= TRUE_THRESH
false_mask = ~true_mask
assert true_mask.sum()  == 528,  f"true_multi {true_mask.sum()} ≠ 528"
assert false_mask.sum() == 2565, f"redundant {false_mask.sum()} ≠ 2565"
print("Source independence locked numbers verified ✓")

pic_true  = df.loc[true_mask,  'pic50_consensus'].values
pic_false = df.loc[false_mask, 'pic50_consensus'].values
score_all = df['source_independence_score'].values

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.5))
ax_hist, ax_dot, ax_kde = axes

# ── Panel A: Score histogram (log y) ─────────────────────────────────────────
ax_hist.hist(score_all, bins=40, color=PAL['blue'], alpha=0.6,
             linewidth=0.3, edgecolor='white')
ax_hist.set_yscale('log')
ax_hist.axvline(TRUE_THRESH, color=PAL['red'], lw=0.75, ls='--')
ymax = ax_hist.get_ylim()[1]
ax_hist.fill_betweenx([1, ymax], TRUE_THRESH, 1.0, color=PAL['red'], alpha=0.06)

ax_hist.text(TRUE_THRESH - 0.02, ymax * 0.4,
             f'<{TRUE_THRESH}: n=2,565\n(pipeline redundancy)',
             ha='right', fontsize=6, color=PAL['gray_dark'],
             fontfamily='sans-serif', multialignment='right')
ax_hist.text(TRUE_THRESH + 0.02, ymax * 0.4,
             f'≥{TRUE_THRESH}: n=528\n(true multi-source)',
             ha='left', fontsize=6, color=PAL['red'],
             fontfamily='sans-serif', multialignment='left')

ax_hist.set_xlabel('Source independence score', fontsize=7)
ax_hist.set_ylabel('Count', fontsize=7)
ax_hist.set_title('Independence score distribution', fontsize=7)
plabel(ax_hist, 'A')

# ── Panel B: Dot plot by source combination ───────────────────────────────────
# Locked values from spec
combo_data = [
    ('BDB only',        1.0, 95),
    ('PC only',         1.0, 233),
    ('ChEMBL only',     1.0, 10),
    ('ChEMBL + PC',     0.7, 23),
    ('BDB + ChEMBL',    0.6, 167),
    ('BDB + PC',        0.5, 1199),
    ('BDB+ChEMBL+PC',   0.3, 1366),
]
# Sort by score descending
combo_data.sort(key=lambda x: x[1], reverse=True)
labels = [c[0] for c in combo_data]
scores = [c[1] for c in combo_data]
ns     = [c[2] for c in combo_data]

# Dot sizes proportional to log(n)
dot_sizes = [max(20, np.log10(n)*40) for n in ns]
dot_colors = [PAL['orange'] if s >= 0.6 else PAL['blue'] for s in scores]

y_pos = range(len(combo_data))
for i, (sc, sz, clr, n) in enumerate(zip(scores, dot_sizes, dot_colors, ns)):
    ax_dot.plot([0, sc], [i, i], color=PAL['gray_light'], lw=0.5, zorder=1)
    ax_dot.scatter([sc], [i], s=sz, c=clr, zorder=3, linewidths=0)
    ax_dot.text(sc + 0.015, i, f'n={n:,}', va='center',
                fontsize=5, color=PAL['gray_dark'], fontfamily='sans-serif')

ax_dot.axvline(TRUE_THRESH, color=PAL['red'], lw=0.5, ls='--', alpha=0.7)
ax_dot.set_yticks(list(y_pos))
ax_dot.set_yticklabels(labels, fontsize=5.5)
ax_dot.set_xlim(-0.05, 1.25)
ax_dot.set_xlabel('Source independence score', fontsize=7)
ax_dot.set_title('Score by source combination', fontsize=7)

handles = [mpatches.Patch(facecolor=PAL['orange'], edgecolor='none', label='Score ≥ 0.6'),
           mpatches.Patch(facecolor=PAL['blue'],   edgecolor='none', label='Score < 0.6')]
ax_dot.legend(handles=handles, loc='lower right', fontsize=5.5, frameon=False)
plabel(ax_dot, 'B')

# ── Panel C: pIC50 KDE by independence class ──────────────────────────────────
xx = np.linspace(2.0, 9.0, 300)
kde_fill(ax_kde, pic_true,  PAL['orange'], f'True multi-source (n=528)',       xx)
kde_fill(ax_kde, pic_false, PAL['blue'],   f'Pipeline redundancy (n=2,565)',   xx)

ax_kde.text(0.03, 0.97, mwu_pval(pic_true, pic_false),
            transform=ax_kde.transAxes, ha='left', va='top',
            fontsize=6, color=PAL['gray_dark'], fontfamily='sans-serif')
ax_kde.legend(loc='upper right', fontsize=5.5)
ax_kde.set_xlabel('pIC50', fontsize=7)
ax_kde.set_ylabel('Density', fontsize=7)
ax_kde.set_title('pIC50 by source independence', fontsize=7)
plabel(ax_kde, 'C')

save_fig(fig, 'fig11_independence')
plt.close()
print("Fig 11 DONE")
