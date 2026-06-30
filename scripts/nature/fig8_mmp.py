#!/usr/bin/env python
"""Nature Fig 8 — MMP Analysis, 2×2."""
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
mmp   = pd.read_csv('outputs/mmp/mmp_pairs_cliff99.csv')
disc  = pd.read_csv('outputs/mmp/mmp_discontinuity_scores.csv')
df    = pd.read_parquet('data/processed/pad4_compounds.parquet')
cliffs = pd.read_parquet('data/processed/activity_cliffs.parquet')

# ── Verification ──────────────────────────────────────────────────────────────
assert len(mmp) == 707, f"MMP pairs {len(mmp)} ≠ 707"
sev = mmp[mmp['cliff_tier']=='severe']
assert len(sev) == 85, f"MMP severe {len(sev)} ≠ 85"
assert (sev['mmp_type']=='single_atom_change').sum() == 49
assert (sev['mmp_type']=='small_substituent').sum()  == 28
assert (sev['mmp_type']=='medium_substituent').sum() == 8
print("All MMP locked numbers verified ✓")

HUB_LABELS = {
    'SMADULGDNOCLOP-GISFHXKWSA-N': 'A1',
    'RAVBZQAQTVGKIV-XBPDSQQVSA-N': 'A2',
    'UDCDEKJNAMHBFH-HSZRJFAPSA-N': 'B1',
    'DVCKJOQIVOGXEI-XMMPIXPASA-N': 'B2',
}
HUB_A = {'SMADULGDNOCLOP-GISFHXKWSA-N','RAVBZQAQTVGKIV-XBPDSQQVSA-N'}
HUB_B = {'UDCDEKJNAMHBFH-HSZRJFAPSA-N','DVCKJOQIVOGXEI-XMMPIXPASA-N'}

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.5))
ax_a, ax_b, ax_c, ax_d = axes.flat

# ── Panel A: Stacked horizontal bar by cliff tier ─────────────────────────────
tier_order  = ['non_cliff', 'broad', 'moderate', 'severe']
mmp_types   = ['single_atom_change', 'small_substituent', 'medium_substituent']
type_colors = {
    'single_atom_change': PAL['blue'],
    'small_substituent':  PAL['teal'],
    'medium_substituent': PAL['orange'],
}
type_labels = {
    'single_atom_change': 'Single R-group change',
    'small_substituent':  'Small substituent',
    'medium_substituent': 'Medium substituent',
}

pivot = (mmp.groupby(['cliff_tier','mmp_type']).size()
           .unstack(fill_value=0).reindex(tier_order, fill_value=0))
for mt in mmp_types:
    if mt not in pivot.columns:
        pivot[mt] = 0

y_pos = np.arange(len(tier_order))
lefts = np.zeros(len(tier_order))
for mt in mmp_types:
    vals = pivot[mt].values
    ax_a.barh(y_pos, vals, left=lefts, color=type_colors[mt],
              height=0.55, linewidth=0, label=type_labels[mt])
    lefts += vals

totals = pivot[mmp_types].sum(axis=1).values
for i, (tot, y) in enumerate(zip(totals, y_pos)):
    ax_a.text(tot + 3, y, str(int(tot)), va='center', fontsize=6,
              color=PAL['gray_dark'], fontfamily='sans-serif')

# Annotate severe bar
sev_idx = tier_order.index('severe')
ax_a.annotate('85/94 severe\ncliffs validated',
              xy=(totals[sev_idx], sev_idx),
              xytext=(totals[sev_idx] + 40, sev_idx),
              fontsize=6, color=PAL['red'], fontfamily='sans-serif',
              va='center', multialignment='left',
              arrowprops=dict(arrowstyle='->', lw=0.5, color=PAL['red']))

ax_a.set_yticks(y_pos)
ax_a.set_yticklabels(['Non-cliff','Broad','Moderate','Severe'], fontsize=6)
ax_a.set_xlabel('MMP pair count', fontsize=7)
ax_a.set_title('MMP pairs by cliff tier and type', fontsize=7)
ax_a.legend(loc='lower right', fontsize=5.5, ncol=1,
            handlelength=0.8, labelspacing=0.2)
plabel(ax_a, 'A')

# ── Panel B: Discontinuity score vs pIC50 ────────────────────────────────────
disc_valid = disc.dropna(subset=['discontinuity_score'])
hub_a_m = disc_valid['inchi_key'].isin(HUB_A)
hub_b_m = disc_valid['inchi_key'].isin(HUB_B)
non_hub  = ~(hub_a_m | hub_b_m)

ax_b.scatter(disc_valid.loc[non_hub, 'pic50_consensus'],
             disc_valid.loc[non_hub, 'discontinuity_score'],
             s=25, c=PAL['gray_light'], alpha=0.7, linewidths=0, zorder=2)
ax_b.scatter(disc_valid.loc[hub_a_m, 'pic50_consensus'],
             disc_valid.loc[hub_a_m, 'discontinuity_score'],
             s=150, c=PAL['navy'], marker='*',
             edgecolors='black', linewidths=0.5, zorder=5, label='Hub A')
ax_b.scatter(disc_valid.loc[hub_b_m, 'pic50_consensus'],
             disc_valid.loc[hub_b_m, 'discontinuity_score'],
             s=150, c=PAL['red'], marker='*',
             edgecolors='black', linewidths=0.5, zorder=5, label='Hub B')

# Labels for all 4 hubs — manual offsets
hub_offsets = {
    'SMADULGDNOCLOP-GISFHXKWSA-N': (-0.3,  0.06),
    'RAVBZQAQTVGKIV-XBPDSQQVSA-N': (-0.3, -0.10),
    'UDCDEKJNAMHBFH-HSZRJFAPSA-N': ( 0.1,  0.06),
    'DVCKJOQIVOGXEI-XMMPIXPASA-N': ( 0.1, -0.10),
}
for ik, lbl in HUB_LABELS.items():
    row = disc_valid[disc_valid['inchi_key']==ik]
    if len(row):
        dx, dy = hub_offsets[ik]
        clr = PAL['navy'] if ik in HUB_A else PAL['red']
        ax_b.annotate(lbl, (row.iloc[0]['pic50_consensus'],
                             row.iloc[0]['discontinuity_score']),
                      xytext=(row.iloc[0]['pic50_consensus']+dx,
                              row.iloc[0]['discontinuity_score']+dy),
                      fontsize=6, fontweight='bold', color=clr,
                      fontfamily='sans-serif',
                      arrowprops=dict(arrowstyle='->', lw=0.4, color=clr,
                                      shrinkA=4, shrinkB=3))

ax_b.set_xlabel('pIC50', fontsize=7)
ax_b.set_ylabel('Mean MMP ΔpIC50 (discontinuity score)', fontsize=7)
ax_b.set_title('MMP discontinuity score vs pIC50', fontsize=7)
ax_b.legend(loc='upper right', fontsize=6, frameon=False)
plabel(ax_b, 'B')

# ── Panel C: ΔpIC50 histogram — Tanimoto vs MMP severe cliffs ─────────────────
tan_sev_dp  = cliffs[cliffs['cliff_tier']=='severe']['delta_pic50'].values
mmp_sev_dp  = mmp[mmp['cliff_tier']=='severe']['delta_pic50'].values
bins = np.linspace(2.0, 3.5, 16)

ax_c.hist(tan_sev_dp, bins=bins, color=PAL['blue'],   alpha=0.5, linewidth=0.3,
          edgecolor='white', label=f'Tanimoto cliffs (n=94)')
ax_c.hist(mmp_sev_dp, bins=bins, color=PAL['orange'], alpha=0.5, linewidth=0.3,
          edgecolor='white', label=f'MMP-validated (n=85)')

# Annotate 9 non-validated cliffs
ax_c.annotate('9 cliffs not\nMMP-validated',
              xy=(tan_sev_dp.max() * 0.97, 3),
              xytext=(3.2, 8),
              fontsize=6, color=PAL['gray_dark'], fontfamily='sans-serif',
              multialignment='center',
              arrowprops=dict(arrowstyle='->', lw=0.5, color=PAL['gray_dark']))

ax_c.set_xlabel('ΔpIC50', fontsize=7)
ax_c.set_ylabel('Count', fontsize=7)
ax_c.set_title('ΔpIC50: Tanimoto vs MMP severe cliffs', fontsize=7)
ax_c.legend(loc='upper right', fontsize=6)
plabel(ax_c, 'C')

# ── Panel D: Top 10 MMP cores by frequency ────────────────────────────────────
core_freq = (mmp.groupby('shared_core').size()
               .sort_values(ascending=False).head(10))
top10_vals = core_freq.values
top10_labels = [f'Core {i+1}' for i in range(len(top10_vals))]

ax_d.barh(range(len(top10_vals))[::-1], top10_vals,
          color=PAL['blue'], linewidth=0, height=0.65)
for i, (cnt, lbl) in enumerate(zip(top10_vals[::-1], top10_labels[::-1])):
    ax_d.text(cnt + 1.5, i, str(int(cnt)), va='center', fontsize=6,
              color=PAL['gray_dark'], fontfamily='sans-serif')

ax_d.set_yticks(range(len(top10_vals)))
ax_d.set_yticklabels(top10_labels[::-1], fontsize=6)
ax_d.set_xlabel('MMP pairs per core', fontsize=7)
ax_d.set_title('Top 10 shared MMP cores', fontsize=7)
ax_d.text(0.02, -0.14, 'Core SMARTS available in Supplementary Table S3.',
          transform=ax_d.transAxes, fontsize=5, color=PAL['gray_dark'],
          ha='left', va='top', style='italic', fontfamily='sans-serif')
plabel(ax_d, 'D')

save_fig(fig, 'fig8_mmp')
plt.close()
print("Fig 8 DONE")
