"""
fig03_potency.py — Figure 3: Potency Distribution (3-panel, DOUBLE width)

Panel A: Histogram + Black KDE. Dashed lines for Mean (6.55) and Median (6.84).
Panel B: Violin plots by source (PubChem, BindingDB, ChEMBL, Patent-only).
Panel C: Horizontal bar chart of Mechanism classes with n count labels.

Outputs: publication/figures/main/fig03_potency.{png,pdf}
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, panel_label, save_fig, SEM, C, DOUBLE

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path

set_style()

CANON = {'n_compounds': 3093, 'mean_pic50': 6.55, 'median_pic50': 6.84}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT  = ROOT / 'publication/figures/main'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FIGURE 3 — POTENCY DISTRIBUTION")
print("=" * 60)

df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
assert len(df) == CANON['n_compounds']

fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 2.7), constrained_layout=True)
ax_a, ax_b, ax_c = axes

# ── Panel A: Global histogram + Black KDE ─────────────────────────────────────
print("[Panel A] Histogram + KDE ...")
pic50 = df['pIC50'].dropna()
bins  = np.arange(2.0, 9.3, 0.30)   # ~24 bins for readability (was 29)
ax_a.hist(pic50, bins=bins, density=True,
          color=SEM['published'], alpha=0.55, edgecolor='white', lw=0.3)

x_kde = np.linspace(2.0, 9.0, 400)
kde_f = stats.gaussian_kde(pic50)
ax_a.plot(x_kde, kde_f(x_kde), color=C['black'], lw=1.2)

mean_v = pic50.mean()
med_v  = pic50.median()
sd_v   = pic50.std()

ax_a.axvline(mean_v, color=SEM['cliff'], lw=1.0, ls='--',
             label=f'Mean {mean_v:.2f}')
ax_a.axvline(med_v,  color=C['gray_dark'], lw=0.8, ls=':',
             label=f'Median {med_v:.2f}')

ax_a.set_xlabel('pIC50')
ax_a.set_ylabel('Density')
ax_a.set_xlim(2.0, 9.0)
ax_a.legend(fontsize=5.5, framealpha=0.85, edgecolor='none', loc='upper left')
ax_a.set_title(f'n = {len(pic50):,}  ·  SD = {sd_v:.2f}', fontsize=6, pad=3)
panel_label(ax_a, 'a')

# ── Panel B: Violin plots by source ──────────────────────────────────────────
print("[Panel B] Violin by source ...")

def src_vals(has_col):
    return df[df[has_col]]['pIC50'].dropna().values

df['has_pubchem'] = df['source_list'].str.contains('pubchem', na=False)
df['has_chembl']  = df['source_list'].str.contains('chembl',  na=False)
df['has_binding'] = df['source_list'].str.contains('bindingdb', na=False)

src_data   = [src_vals('has_pubchem'), src_vals('has_binding'), src_vals('has_chembl'),
              df[df['patent_flag']]['pIC50'].dropna().values]
src_labels = ['PubChem', 'BindingDB', 'ChEMBL', 'Patent\nonly']
src_colors = [C['blue'], C['teal'], C['navy'], C['orange']]
src_ns     = [len(d) for d in src_data]

vp = ax_b.violinplot(src_data, positions=range(len(src_data)),
                      showmedians=False, showextrema=False)
for pc, color in zip(vp['bodies'], src_colors):
    pc.set_facecolor(color)
    pc.set_alpha(0.60)
    pc.set_edgecolor(C['black'])
    pc.set_linewidth(0.4)

for i, data in enumerate(src_data):
    q1, med, q3 = np.percentile(data, [25, 50, 75])
    # IQR box + median tick so violins are not just colored blobs
    ax_b.vlines(i, q1, q3, color=C['black'], lw=3.2, zorder=3, alpha=0.55)
    ax_b.hlines(med, i - 0.28, i + 0.28, color=C['black'], lw=1.6, zorder=4)
    ax_b.hlines([q1, q3], i - 0.14, i + 0.14, color=C['black'], lw=0.7, zorder=4)

# Mann-Whitney: PubChem vs patent
mw_stat, mw_p = stats.mannwhitneyu(src_data[0], src_data[3], alternative='two-sided')
p_str = 'p < 0.001' if mw_p < 0.001 else f'p = {mw_p:.3f}'

ax_b.set_xticks(range(len(src_data)))
ax_b.set_xticklabels([f'{l}\nn={n:,}' for l, n in zip(src_labels, src_ns)], fontsize=5.5)
ax_b.set_ylabel('pIC50')
ax_b.set_title(f'PubChem vs Patent  {p_str}', fontsize=6, pad=3)
panel_label(ax_b, 'b')

# ── Panel C: Horizontal bar chart of Mechanism classes ────────────────────────
print("[Panel C] Mechanism class bar ...")
mech_order  = ['enzymatic', 'enzymatic_confirmed', 'fp_ic50', 'covalent']
mech_labels = ['Enzymatic\n(BAEE)', 'Enz. confirmed\n(RFMS)', 'FP binding', 'Covalent']
# Use blue for the dominant enzymatic class — grey (#BBBBBB) is invisible on white background
mech_colors = [C['blue'], SEM['enzymatic_confirmed'], SEM['fp_ic50'], SEM['covalent']]

mech_ns    = [int((df['mechanism_class'] == m).sum()) for m in mech_order]
mech_means = [df[df['mechanism_class'] == m]['pIC50'].mean() for m in mech_order]

y_pos = list(range(len(mech_order)))
bars_c = ax_c.barh(y_pos, mech_ns,
                    color=mech_colors, height=0.55,
                    edgecolor='white', lw=0.4, alpha=0.85)

# Annotate n and mean pIC50
for i, (n, mean_m) in enumerate(zip(mech_ns, mech_means)):
    ax_c.text(n + max(mech_ns) * 0.01, i,
              f'n={n:,}  (mean {mean_m:.2f})',
              va='center', ha='left', fontsize=5.5)

ax_c.set_yticks(y_pos)
ax_c.set_yticklabels(mech_labels, fontsize=6)
ax_c.set_xlabel('Compounds')
ax_c.set_xlim(0, max(mech_ns) * 1.65)
ax_c.invert_yaxis()
panel_label(ax_c, 'c')

# ── Save ──────────────────────────────────────────────────────────────────────
save_fig(fig, str(OUT / 'fig03_potency'))
plt.close(fig)
print("Figure 3 complete.")
