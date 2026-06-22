"""
fig06_mmp.py — Figure 6: MMP Analysis (4-panel, DOUBLE width)
Outputs: publication/figures/main/fig06_mmp.{png,pdf}
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, panel_label, save_fig, SEM, C, SINGLE, ONEHALF, DOUBLE

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

set_style()

HUB_IKS = {
    'A1': 'SMADULGDNOCLOP-GISFHXKWSA-N',
    'A2': 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
    'B1': 'UDCDEKJNAMHBFH-HSZRJFAPSA-N',
    'B2': 'DVCKJOQIVOGXEI-XMMPIXPASA-N',
}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'publication/figures/main'
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

fig, axes = plt.subplots(2, 2, figsize=(DOUBLE, 5.0), constrained_layout=True)
ax_a, ax_b, ax_c, ax_d = axes.flatten()

# ── Panel a: MMP pairs by cliff tier and mmp_type ─────────────────────────────
print("\n[Panel a] MMP pairs by tier/type ...")
tiers = ['severe', 'moderate', 'broad']
mmp_types = ['single_atom_change', 'small_substituent', 'medium_substituent']
type_colors = [C['navy'], C['teal'], C['orange']]
type_labels = ['Single atom', 'Small substituent', 'Medium substituent']

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
                       str(int(val)), ha='center', va='bottom', fontsize=5.5)

ax_a.set_xticks(x)
ax_a.set_xticklabels(['Severe\n(n=80)', 'Moderate\n(n=25)', 'Broad\n(n=2)'], fontsize=6)
ax_a.set_ylabel('MMP pairs')
ax_a.legend(fontsize=5.5, framealpha=0.7, edgecolor='none')
panel_label(ax_a, 'a')

# ── Panel b: MMP discontinuity score vs pIC50 ─────────────────────────────────
print("[Panel b] Discontinuity score vs pIC50 ...")
hub_iks_set = set(HUB_IKS.values())
non_hub = disc_joined[~disc_joined['inchi_key'].isin(hub_iks_set)]
hub_a = disc_joined[disc_joined['inchi_key'].isin({HUB_IKS['A1'], HUB_IKS['A2']})]
hub_b = disc_joined[disc_joined['inchi_key'].isin({HUB_IKS['B1'], HUB_IKS['B2']})]

ax_b.scatter(non_hub['pIC50'], non_hub['discontinuity_score'],
              s=18, c=C['grey'], alpha=0.5, zorder=2)
ax_b.scatter(hub_a['pIC50'], hub_a['discontinuity_score'],
              s=150, c=SEM['classA'], marker='*', zorder=5, edgecolors='white', lw=0.5)
ax_b.scatter(hub_b['pIC50'], hub_b['discontinuity_score'],
              s=110, c=SEM['classB'], marker='D', zorder=5, edgecolors='white', lw=0.5)

for lbl, ik in HUB_IKS.items():
    row = disc_joined[disc_joined['inchi_key'] == ik]
    if len(row) > 0:
        x_pt, y_pt = row['pIC50'].values[0], row['discontinuity_score'].values[0]
        color = SEM['classA'] if lbl.startswith('A') else SEM['classB']
        x_off = (+10, +10) if lbl != 'B2' else (+10, +10)
        ax_b.annotate(lbl, (x_pt, y_pt), xytext=(10, 10), textcoords='offset points',
                       fontsize=5.5, fontweight='bold', color=color)

ax_b.set_xlabel('pIC50')
ax_b.set_ylabel('Discontinuity score')
panel_label(ax_b, 'b')

# ── Panel c: DeltapIC50 all severe vs MMP-validated ─────────────────────────
print("[Panel c] DeltapIC50 distribution ...")
cliffs_all = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
severe_all_94 = cliffs_all[cliffs_all['cliff_tier'] == 'severe'].copy()

ecfp_ct = pd.read_csv(ROOT / 'outputs/stress_test/ECFP4_vs_ECFP6_all94_pairs.csv')

def _norm_pair(a, b):
    return tuple(sorted([a, b]))

ecfp_ct['pair'] = ecfp_ct.apply(lambda r: _norm_pair(r['ik_a'], r['ik_b']), axis=1)
severe_all_94['pair'] = severe_all_94.apply(
    lambda r: _norm_pair(r['inchi_key_a'], r['inchi_key_b']), axis=1)
mmp_validated_pairs = set(ecfp_ct.loc[ecfp_ct['mmp_validated'], 'pair'])
severe_all_94['mmp_validated_flag'] = severe_all_94['pair'].isin(mmp_validated_pairs)
severe_mmp = severe_all_94[severe_all_94['mmp_validated_flag']]
n_non = len(severe_all_94) - len(severe_mmp)
print(f"  All severe: {len(severe_all_94)} | MMP-validated: {len(severe_mmp)} | Non: {n_non}")

delta_all = severe_all_94['delta_pic50'].abs()
delta_mmp = severe_mmp['delta_pic50'].abs()

bins_c = np.linspace(delta_all.min() * 0.9, delta_all.max() * 1.05, 16)
ax_c.hist(delta_all, bins=bins_c, alpha=0.35, color=C['blue'],
           label=f'All severe (n={len(severe_all_94)})')
ax_c.hist(delta_mmp, bins=bins_c, alpha=0.80, color=C['teal'],
           label=f'MMP-validated (n={len(severe_mmp)})')

ax_c.set_xlabel('|ΔpIC50|')
ax_c.set_ylabel('Pairs')
ax_c.set_title(f'{n_non} of {len(severe_all_94)} cliffs not MMP-validated', fontsize=6, pad=3)
ax_c.legend(fontsize=6, framealpha=0.85, edgecolor='none', loc='upper right')
panel_label(ax_c, 'c')

# ── Panel d: Top 10 shared MMP cores ─────────────────────────────────────────
print("[Panel d] Top 10 shared MMP cores ...")
core_counts = mmp.groupby('shared_core').size().sort_values(ascending=False).head(10)
print(f"  Top cores:\n{core_counts}")

y_labels_d = [f'Core {i+1}' for i in range(len(core_counts))]
colors_d = [C['navy'], C['blue'], C['blue']] + [C['teal']] * 2 + [C['grey']] * 5
y_pos_d = range(len(core_counts))

ax_d.barh(list(y_pos_d), core_counts.values,
           color=colors_d[:len(core_counts)], height=0.6,
           edgecolor='white', lw=0.3)

x_max = core_counts.values[0]  # largest bar (Core 1)
for i, val in enumerate(core_counts.values):
    if i == 0 and val > 0.65 * x_max:
        # Place label inside the longest bar
        ax_d.text(val - 2, i, str(int(val)), va='center', ha='right', fontsize=5.5, color='white')
    else:
        ax_d.text(val + 0.2, i, str(int(val)), va='center', ha='left', fontsize=5.5)

ax_d.set_yticks(list(y_pos_d))
ax_d.set_yticklabels(y_labels_d, fontsize=6)
ax_d.text(0.98, 0.02, '† Cores ranked by pairs shared',
           transform=ax_d.transAxes, fontsize=5, ha='right', va='bottom',
           color=C['grey'], style='italic')
ax_d.set_xlabel('MMP pairs')
ax_d.invert_yaxis()
panel_label(ax_d, 'd')

# ── Save ──────────────────────────────────────────────────────────────────────
save_fig(fig, str(OUT / 'fig06_mmp'))
plt.close(fig)
print("\nFigure 6 complete.")
