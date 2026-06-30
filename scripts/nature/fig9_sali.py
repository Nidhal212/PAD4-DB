#!/usr/bin/env python
"""Nature Fig 9 — SALI Analysis, 1×3."""
import os, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

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

# ── Load ──────────────────────────────────────────────────────────────────────
pairs = pd.read_parquet('data/processed/activity_pairs_with_sali.parquet')
pairs_valid = pairs[pairs['sali'].notna() & (pairs['sali'] > 0.01)].copy()

tier_colors = {
    'severe':   PAL['red'],
    'moderate': PAL['orange'],
    'broad':    PAL['blue'],
    'non_cliff': PAL['gray_light'],
}

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.5))
ax_a, ax_b, ax_c = axes

# ── Panel A: Top 20 SALI pairs ────────────────────────────────────────────────
top20 = pairs.nlargest(20, 'sali').copy().reset_index(drop=True)
top20['tier'] = top20['cliff_tier'].fillna('non_cliff')
top20['label'] = [f'Pair {i+1}' for i in range(len(top20))]
bar_clrs = [tier_colors.get(t, PAL['gray_light']) for t in top20['tier']]

ax_a.barh(range(len(top20))[::-1], top20['sali'].values,
          color=bar_clrs, linewidth=0, height=0.7)
for i, (sali_val, tier) in enumerate(zip(top20['sali'][::-1].values,
                                         top20['tier'][::-1].values)):
    ax_a.text(sali_val + 0.5, i, f'{sali_val:.2f}', va='center', fontsize=5,
              color=PAL['gray_dark'], fontfamily='sans-serif')

ax_a.set_yticks(range(len(top20)))
ax_a.set_yticklabels([f'Pair {20-i}' for i in range(20)], fontsize=5.5)
ax_a.set_xlabel('SALI', fontsize=7)
ax_a.set_title('Top 20 SALI pairs', fontsize=7)

# Annotate pair 1 (highest SALI)
ax_a.annotate('SALI=65.88',
              xy=(top20['sali'].iloc[0], 19),
              xytext=(top20['sali'].iloc[0] - 20, 17),
              fontsize=6, color=PAL['red'], fontfamily='sans-serif',
              arrowprops=dict(arrowstyle='->', lw=0.5, color=PAL['red']))

# Legend
handles = [mpatches.Patch(facecolor=c, edgecolor='none', label=t.capitalize())
           for t, c in tier_colors.items() if t != 'non_cliff']
ax_a.legend(handles=handles, loc='lower right', fontsize=5.0,
            handlelength=0.7, labelspacing=0.2, frameon=False)
ax_a.text(0.02, -0.12, 'Pair identities in Supplementary Table S4.',
          transform=ax_a.transAxes, fontsize=5, color=PAL['gray_dark'],
          ha='left', va='top', style='italic', fontfamily='sans-serif')
plabel(ax_a, 'A')

# ── Panel B: SALI vs ΔpIC50 (sampled, Tanimoto color) ────────────────────────
sample = pairs_valid.sample(n=min(20000, len(pairs_valid)), random_state=42)
sc = ax_b.scatter(sample['delta_pic50'], sample['sali'],
                  s=2, c=sample['tanimoto'], cmap='viridis_r',
                  vmin=0.6, vmax=1.0, alpha=0.3, linewidths=0, rasterized=True)
cb = fig.colorbar(sc, ax=ax_b, pad=0.02, shrink=0.85, aspect=25)
cb.set_label('Tanimoto', fontsize=6)
cb.ax.tick_params(labelsize=5, width=0.4)

# Overlay SALI > 20 (n=19)
high_sali = pairs[pairs['sali'] > 20]
ax_b.scatter(high_sali['delta_pic50'], high_sali['sali'],
             s=20, c=PAL['red'], linewidths=0.3, edgecolors='white',
             zorder=5, label=f'SALI > 20 (n={len(high_sali)})')

ax_b.axhline(10, color=PAL['gray_dark'], lw=0.5, ls='--', alpha=0.6)
ax_b.set_yscale('log')
ax_b.set_xlabel('ΔpIC50', fontsize=7)
ax_b.set_ylabel('SALI (log scale)', fontsize=7)
ax_b.set_title('SALI vs ΔpIC50', fontsize=7)
ax_b.legend(loc='upper left', fontsize=5.5)
plabel(ax_b, 'B')

# ── Panel C: Cumulative SALI distribution (ECDF) ─────────────────────────────
all_sali = pairs['sali'].dropna().values
sorted_s  = np.sort(all_sali)
ecdf_y    = np.arange(1, len(sorted_s)+1) / len(sorted_s)
ax_c.plot(sorted_s, ecdf_y, color=PAL['blue'], lw=0.75)
ax_c.set_xlim(0, 70)
ax_c.set_ylim(0, 1)

ax_c.axvline(10, color=PAL['orange'], lw=0.5, ls='--')
ax_c.axvline(20, color=PAL['red'],    lw=0.5, ls='--')

frac_10 = (all_sali <= 10).mean()
frac_20 = (all_sali <= 20).mean()
ax_c.text(10.5, frac_10 + 0.04, f'335 pairs\n(>10)',
          fontsize=6, color=PAL['orange'], fontfamily='sans-serif',
          multialignment='center')
ax_c.text(20.5, frac_20 + 0.04, f'19 pairs\n(>20)',
          fontsize=6, color=PAL['red'], fontfamily='sans-serif',
          multialignment='center')

ax_c.set_xlabel('SALI', fontsize=7)
ax_c.set_ylabel('Cumulative fraction of pairs', fontsize=7)
ax_c.set_title('SALI cumulative distribution', fontsize=7)
plabel(ax_c, 'C')

save_fig(fig, 'fig9_sali')
plt.close()
print("Fig 9 DONE")
