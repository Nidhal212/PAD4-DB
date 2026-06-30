"""
fig06_mmp.py — Figure 6: MMP Analysis (4-panel)
Outputs: outputs/figures/nature/fig06_mmp.{png,pdf}
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
    'cyan':       '#33BBEE',
    'red':        '#CC3311',
    'navy':       '#004488',
    'grey':       '#BBBBBB',
    'dark_grey':  '#555555',
    'light_grey': '#E8E8E8',
}

HUB_IKS = {
    'A1': 'SMADULGDNOCLOP-GISFHXKWSA-N',
    'A2': 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
    'B1': 'UDCDEKJNAMHBFH-HSZRJFAPSA-N',
    'B2': 'DVCKJOQIVOGXEI-XMMPIXPASA-N',
}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/figures/nature_v2'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FIGURE 6 — MMP ANALYSIS")
print("=" * 60)

df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
mmp = pd.read_csv(ROOT / 'outputs/mmp/mmp_pairs_cliff99.csv')
disc = pd.read_csv(ROOT / 'outputs/mmp/mmp_discontinuity_scores.csv')

print(f"  MMP pairs: {len(mmp)}")
print(f"  MMP cliff_tier: {mmp['cliff_tier'].value_counts().to_dict()}")
print(f"  MMP mmp_type: {mmp['mmp_type'].value_counts().to_dict()}")

disc_joined = disc.merge(df[['inchi_key', 'pIC50']], on='inchi_key', how='left')

fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
ax_a, ax_b, ax_c, ax_d = axes.flatten()

# ── Panel a: MMP pairs by cliff tier and mmp_type ─────────────────────────────
print("\n[Panel a] MMP pairs by tier/type ...")
tiers = ['severe', 'moderate', 'broad']
mmp_types = ['single_atom_change', 'small_substituent', 'medium_substituent']
type_colors = [COLORS['blue'], COLORS['teal'], COLORS['orange']]
type_labels = ['Single atom change', 'Small substituent', 'Medium substituent']

# Filter to cliff tiers only
mmp_cliff = mmp[mmp['cliff_tier'].isin(tiers)]
grouped = mmp_cliff.groupby(['cliff_tier', 'mmp_type']).size().unstack(fill_value=0)
print(f"  Grouped:\n{grouped}")

x = np.arange(len(tiers))
width = 0.25

for i, (mtype, color, label) in enumerate(zip(mmp_types, type_colors, type_labels)):
    vals = [grouped.get(mtype, pd.Series()).get(t, 0) for t in tiers]
    bars = ax_a.bar(x + (i - 1) * width, vals, width=width * 0.9,
                     color=color, alpha=0.85, label=label, edgecolor='white', lw=0.3)
    for bar, val in zip(bars, vals):
        if val > 0:
            ax_a.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                       str(int(val)), ha='center', va='bottom', fontsize=7)

ax_a.set_xticks(x)
ax_a.set_xticklabels(['Severe\n(n=85)', 'Moderate\n(n=25)', 'Broad\n(n=2)'], fontsize=8)
ax_a.set_ylabel('Number of MMP pairs', fontsize=9)
ax_a.legend(fontsize=7.5, framealpha=0.7, edgecolor='none')
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)

# Annotation
ax_a.annotate('85/94 severe\ncliffs MMP-\nvalidated',
               xy=(0, grouped.get('single_atom_change', pd.Series()).get('severe', 0)),
               xytext=(0.35, 0.85), textcoords='axes fraction',
               arrowprops=dict(arrowstyle='->', color=COLORS['red'], lw=1.0),
               fontsize=8, color=COLORS['red'])

ax_a.text(-0.13, 1.04, 'a', transform=ax_a.transAxes, fontsize=11, fontweight='bold')

# ── Panel b: MMP discontinuity score vs pIC50 ─────────────────────────────────
print("[Panel b] Discontinuity score vs pIC50 ...")
hub_iks_set = set(HUB_IKS.values())
non_hub = disc_joined[~disc_joined['inchi_key'].isin(hub_iks_set)]
hub_a = disc_joined[disc_joined['inchi_key'].isin({HUB_IKS['A1'], HUB_IKS['A2']})]
hub_b = disc_joined[disc_joined['inchi_key'].isin({HUB_IKS['B1'], HUB_IKS['B2']})]

ax_b.scatter(non_hub['pIC50'], non_hub['discontinuity_score'],
              s=40, c=COLORS['grey'], alpha=0.5, zorder=2)
ax_b.scatter(hub_a['pIC50'], hub_a['discontinuity_score'],
              s=250, c=COLORS['navy'], marker='*', zorder=5, edgecolors='white', lw=0.5)
ax_b.scatter(hub_b['pIC50'], hub_b['discontinuity_score'],
              s=200, c=COLORS['red'], marker='D', zorder=5, edgecolors='white', lw=0.5)

for lbl, ik in HUB_IKS.items():
    row = disc_joined[disc_joined['inchi_key'] == ik]
    if len(row) > 0:
        x, y = row['pIC50'].values[0], row['discontinuity_score'].values[0]
        color = COLORS['navy'] if lbl.startswith('A') else COLORS['red']
        ax_b.annotate(lbl, (x, y), xytext=(6, 4), textcoords='offset points',
                       fontsize=8, fontweight='bold', color=color,
                       arrowprops=dict(arrowstyle='-', color=color, lw=0.8))

ax_b.set_xlabel('pIC50', fontsize=9)
ax_b.set_ylabel('Discontinuity score', fontsize=9)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.text(-0.13, 1.04, 'b', transform=ax_b.transAxes, fontsize=11, fontweight='bold')

# ── Panel c: DeltapIC50: all severe vs MMP-validated ─────────────────────────
print("[Panel c] DeltapIC50 distribution ...")
# Load all 94 severe activity cliffs (incl. those not in MMP pairs file)
cliffs_all = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
severe_all_94 = cliffs_all[cliffs_all['cliff_tier'] == 'severe']
# MMP-validated subset: 85 canonical severe cliffs confirmed by MMP
severe_mmp = mmp[mmp['is_canonical_severe_cliff'] == True]
n_non = len(severe_all_94) - len(severe_mmp)
print(f"  All severe: {len(severe_all_94)} | MMP-validated: {len(severe_mmp)} | Non-validated: {n_non}")

delta_all = severe_all_94['delta_pic50'].abs()
delta_mmp = severe_mmp['delta_pic50'].abs()

bins_c = np.linspace(delta_all.min() * 0.9, delta_all.max() * 1.05, 16)
ax_c.hist(delta_all, bins=bins_c, alpha=0.4, color=COLORS['blue'],
           label=f'All severe (n={len(severe_all_94)})')
ax_c.hist(delta_mmp, bins=bins_c, alpha=0.7, color=COLORS['orange'],
           label=f'MMP-validated (n={len(severe_mmp)})')

if n_non > 0:
    ax_c.text(0.97, 0.97, f'{n_non} cliffs\nnot MMP-\nvalidated',
               transform=ax_c.transAxes, fontsize=7.5, ha='right', va='top',
               color=COLORS['red'],
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                         edgecolor=COLORS['red'], lw=0.8))

ax_c.set_xlabel('|ΔpIC50|', fontsize=9)
ax_c.set_ylabel('Number of pairs', fontsize=9)
ax_c.legend(fontsize=7.5, framealpha=0.7, edgecolor='none')
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
ax_c.text(-0.13, 1.04, 'c', transform=ax_c.transAxes, fontsize=11, fontweight='bold')

# ── Panel d: Top 10 shared MMP cores ─────────────────────────────────────────
print("[Panel d] Top 10 shared MMP cores ...")
core_counts = mmp.groupby('shared_core').size().sort_values(ascending=False).head(10)
print(f"  Top cores:\n{core_counts}")

y_labels_d = [f'Core {i+1}' for i in range(len(core_counts))]
colors_d = [COLORS['navy']] + [COLORS['blue']] * 2 + [COLORS['grey']] * 7
y_pos_d = range(len(core_counts))

ax_d.barh(list(y_pos_d), core_counts.values,
           color=colors_d[:len(core_counts)], height=0.6,
           edgecolor='white', lw=0.3)

for i, val in enumerate(core_counts.values):
    ax_d.text(val + 0.2, i, str(int(val)), va='center', ha='left', fontsize=7.5)

ax_d.set_yticks(list(y_pos_d))
ax_d.set_yticklabels(y_labels_d, fontsize=8)
ax_d.text(0.98, 0.02, '† Core labels rank shared MMP scaffolds\nby number of pairs',
           transform=ax_d.transAxes, fontsize=6, ha='right', va='bottom',
           color=COLORS['dark_grey'], style='italic')
ax_d.set_xlabel('Number of MMP pairs', fontsize=9)
ax_d.invert_yaxis()
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)
ax_d.text(-0.22, 1.04, 'd', transform=ax_d.transAxes, fontsize=11, fontweight='bold')

# ── Save ──────────────────────────────────────────────────────────────────────
for ext in ['png', 'pdf']:
    outpath = OUT / f'fig06_mmp.{ext}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"  Saved: {outpath}")

plt.close(fig)
print("\nFigure 6 complete.")
