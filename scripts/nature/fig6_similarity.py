#!/usr/bin/env python
"""Nature Fig 6 — Similarity Landscape (4 panels)."""
import os, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
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
    'teal': '#009988', 'navy': '#1A237E',
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

# ── Load data ─────────────────────────────────────────────────────────────────
pairs = pd.read_parquet('data/processed/activity_pairs_with_sali.parquet')
cliffs = pd.read_parquet('data/processed/activity_cliffs.parquet')

# Merge SALI onto cliffs
if 'sali' not in cliffs.columns:
    if 'sali' in pairs.columns:
        cliffs = cliffs.merge(
            pairs[['inchi_key_a','inchi_key_b','sali']].drop_duplicates(),
            on=['inchi_key_a','inchi_key_b'], how='left'
        )

tan_all  = pairs['tanimoto'].values
dpic_all = pairs['delta_pic50'].values
sali_all = pairs['sali'].dropna().values if 'sali' in pairs.columns else np.array([])

tier_colors = {'severe': PAL['red'], 'moderate': PAL['orange'], 'broad': PAL['blue']}

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.5))
ax_th, ax_sc, ax_sl, ax_sas = axes.flat

# ── Panel A: Tanimoto histogram ───────────────────────────────────────────────
bins_tan = np.linspace(0.6, 1.0, 41)
ax_th.hist(tan_all, bins=bins_tan, color=PAL['blue'], alpha=0.6, linewidth=0, density=False)
ax_th.axvline(0.8, color=PAL['red'], lw=0.8, ls='--', label='Cliff threshold (0.8)')
ax_th.set_xlabel('Tanimoto similarity (ECFP4)', fontsize=7)
ax_th.set_ylabel('Pair count', fontsize=7)
ax_th.set_title(f'Tanimoto distribution\n({len(tan_all):,} pairs, sim≥0.6)', fontsize=7)
ax_th.legend(fontsize=5.5)
ax_th.text(0.97, 0.97, f'sim≥0.8: 12,071 pairs\n(3.4% of sim≥0.6)',
           transform=ax_th.transAxes, ha='right', va='top',
           fontsize=5.5, color=PAL['gray_dark'], fontfamily='sans-serif')
plabel(ax_th, 'A')

# ── Panel B: ΔpIC50 vs Tanimoto scatter (cliff tiers) ─────────────────────────
for tier, clr in tier_colors.items():
    sub = cliffs[cliffs['cliff_tier'] == tier]
    ax_sc.scatter(sub['tanimoto'], sub['delta_pic50'], s=4, c=clr,
                  alpha=0.65, linewidths=0, label=f'{tier.capitalize()} (n={len(sub)})',
                  rasterized=True)
ax_sc.axhline(1.0, color=PAL['gray_light'], lw=0.5, ls='--')
ax_sc.axhline(1.5, color=PAL['gray_light'], lw=0.5, ls=':')
ax_sc.axhline(2.0, color=PAL['gray_dark'],  lw=0.6, ls='--')
ax_sc.axvline(0.8, color=PAL['gray_dark'],  lw=0.6, ls='--')
ax_sc.set_xlabel('Tanimoto similarity (ECFP4)', fontsize=7)
ax_sc.set_ylabel('ΔpIC50', fontsize=7)
ax_sc.set_title('Activity cliff landscape', fontsize=7)
ax_sc.legend(loc='upper left', fontsize=5.5)

# Quadrant labels with patheffects
pefx = [pe.withStroke(linewidth=1.5, foreground='white')]
ax_sc.text(0.82, 0.6, 'broad', fontsize=5.5, color=PAL['blue'],
           path_effects=pefx, fontfamily='sans-serif')
plabel(ax_sc, 'B')

# ── Panel C: SALI histogram — linear x, log y ────────────────────────────────
if len(sali_all) > 0:
    sali_pos = sali_all[sali_all > 0]
    bins_sali = np.linspace(0, min(sali_pos.max(), 70), 50)
    ax_sl.hist(sali_pos, bins=bins_sali, color=PAL['teal'], alpha=0.7, linewidth=0)
    ax_sl.set_yscale('log')
    ax_sl.axvline(10, color=PAL['red'], lw=0.8, ls='--', label='SALI > 10 (n=335)')
    ax_sl.axvline(20, color=PAL['orange'], lw=0.8, ls=':', label='SALI > 20 (n=19)')
    ax_sl.legend(fontsize=5.5)
    ax_sl.text(0.97, 0.97, 'Max SALI: 65.88\nSALI>10: 335\nSALI>20: 19',
               transform=ax_sl.transAxes, ha='right', va='top',
               fontsize=5.5, color=PAL['gray_dark'], fontfamily='sans-serif')
    ax_sl.set_xlabel('SALI', fontsize=7)
    ax_sl.set_ylabel('Pair count (log scale)', fontsize=7)  # explicitly set
    ax_sl.set_title(f'SALI distribution (n={len(sali_pos):,} SALI > 0)', fontsize=7)
else:
    ax_sl.text(0.5, 0.5, 'SALI not available', transform=ax_sl.transAxes,
               ha='center', va='center', fontsize=7, color=PAL['gray_dark'])
plabel(ax_sl, 'C')

# ── Panel D: SALI landscape (SALI vs Tanimoto, colored by cliff tier) ────────
if 'sali' in pairs.columns:
    p_sali = pairs[pairs['sali'].notna() & (pairs['sali'] > 0.01)].copy()  # remove near-zero noise floor
    cliff_tier_map = cliffs[['inchi_key_a','inchi_key_b','cliff_tier']].drop_duplicates()
    p_sali = p_sali.merge(cliff_tier_map, on=['inchi_key_a','inchi_key_b'], how='left')
    # handle suffix if tanimoto or other cols conflicted
    ct_col = [c for c in p_sali.columns if 'cliff_tier' in c]
    if ct_col:
        p_sali['cliff_tier'] = p_sali[ct_col[0]].fillna('non_cliff')
    else:
        p_sali['cliff_tier'] = 'non_cliff'
    tier_cols_d = {
        'non_cliff': PAL['gray_light'],
        'broad':     PAL['blue'],
        'moderate':  PAL['orange'],
        'severe':    PAL['red'],
    }
    for tier in ['non_cliff', 'broad', 'moderate', 'severe']:
        sub = p_sali[p_sali['cliff_tier'] == tier]
        lbl = tier.replace('_', ' ').capitalize()
        ax_sas.scatter(sub['tanimoto'], sub['sali'],
                       s=1.5, c=tier_cols_d[tier], alpha=0.45,
                       linewidths=0, rasterized=True, label=f'{lbl} (n={len(sub):,})')
    ax_sas.set_yscale('log')
    ax_sas.axhline(10, color=PAL['gray_dark'], lw=0.5, ls='--', alpha=0.6)
    ax_sas.axvline(0.8, color=PAL['gray_dark'], lw=0.5, ls='--', alpha=0.6)
    ax_sas.set_xlabel('Tanimoto similarity (ECFP4)', fontsize=7)
    ax_sas.set_ylabel('SALI (log scale)', fontsize=7)
    ax_sas.set_title('SALI landscape', fontsize=7)
    ax_sas.legend(loc='upper left', fontsize=5.0, frameon=False,
                  markerscale=3, handlelength=0.6, labelspacing=0.2)
else:
    ax_sas.text(0.5, 0.5, 'SALI data not available', transform=ax_sas.transAxes,
                ha='center', va='center', fontsize=7, color=PAL['gray_dark'])
plabel(ax_sas, 'D')

save_fig(fig, 'fig6_similarity')
plt.close()
print("Fig 6 DONE")
