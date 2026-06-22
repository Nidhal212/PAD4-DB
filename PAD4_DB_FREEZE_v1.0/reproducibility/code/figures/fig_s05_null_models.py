"""
fig_s05_null_models.py — Supplementary Figure S5.

Observed severe-cliff count and hub share vs three permutation nulls
(unrestricted, scaffold-constrained, assay-constrained). Data: null_model_comparison.csv.
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, panel_label, save_fig, SEM, C, DOUBLE
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

set_style()
ROOT = Path('/home/nidhal/PAD4-db_V2')
df = pd.read_csv(ROOT / 'outputs/audit/null_model_comparison.csv')
labels = ['Unrestricted', 'Scaffold-\nconstrained', 'Assay-\nconstrained']
x = np.arange(3)
cols = [C['grey'], C['navy'], C['teal']]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(DOUBLE * 0.78, 2.8), constrained_layout=True)

# Panel A: severe cliff count
axA.bar(x, df.null_cliffs_mean, yerr=df.null_cliffs_sd, color=cols, width=0.6,
        edgecolor='white', lw=0.5, error_kw=dict(lw=0.8, capsize=3))
axA.axhline(df.obs_cliffs.iloc[0], color=SEM['cliff'], lw=1.2, ls='--', zorder=5)
axA.text(2.4, df.obs_cliffs.iloc[0] + 40, f'observed = {int(df.obs_cliffs.iloc[0])}',
         color=SEM['cliff'], fontsize=6, ha='right', fontweight='bold')
for i, r in df.iterrows():
    p = 'n.s.' if r.cliff_p > 0.05 else f'{r.cliff_fold_depletion:.0f}×'
    axA.text(i, r.null_cliffs_mean + r.null_cliffs_sd + 30, p, ha='center', fontsize=6,
             color=(C['gray_dark'] if r.cliff_p > 0.05 else C['black']))
axA.set_xticks(x); axA.set_xticklabels(labels, fontsize=6)
axA.set_ylabel('Severe cliffs (null mean ± SD)')
axA.set_title('Cliff rarity', fontsize=7, pad=3)
panel_label(axA, 'a', x=-0.16, y=1.04)

# Panel B: hub share
axB.bar(x, df.null_hub_share_mean_pct, yerr=df.null_hub_share_sd_pct, color=cols, width=0.6,
        edgecolor='white', lw=0.5, error_kw=dict(lw=0.8, capsize=3))
axB.axhline(df.obs_hub_share_pct.iloc[0], color=SEM['cliff'], lw=1.2, ls='--', zorder=5)
axB.text(2.4, df.obs_hub_share_pct.iloc[0] + 1.5, f'observed = {df.obs_hub_share_pct.iloc[0]:.1f}%',
         color=SEM['cliff'], fontsize=6, ha='right', fontweight='bold')
for i, r in df.iterrows():
    pstr = 'P<0.001' if r.hub_p < 0.001 else f'P={r.hub_p:.3f}'
    axB.text(i, r.null_hub_share_mean_pct + r.null_hub_share_sd_pct + 1.5,
             f'{r.hub_fold_enrichment:.1f}×\n{pstr}', ha='center', fontsize=5.5, color=C['black'])
axB.set_xticks(x); axB.set_xticklabels(labels, fontsize=6)
axB.set_ylabel('Top-4 hub share of cliffs (%)')
axB.set_ylim(0, 62)
axB.set_title('Hub concentration', fontsize=7, pad=3)
panel_label(axB, 'b', x=-0.16, y=1.04)

save_fig(fig, str(ROOT / 'publication/figures/supplementary/fig_s05_null_models'))
plt.close(fig)
print("Figure S5 complete.")
