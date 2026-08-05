"""
fig03_potency.py — Figure 3: Potency Distribution (4-panel)
Outputs: outputs/figures/nature/fig03_potency.{png,pdf}
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
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
    'magenta':    '#EE3377',
    'red':        '#CC3311',
    'navy':       '#004488',
    'grey':       '#BBBBBB',
    'dark_grey':  '#555555',
    'light_grey': '#E8E8E8',
}

CANON = {'n_compounds': 3093, 'n_patent': 233, 'n_multi_06': 528, 'n_in_severe': 99}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/figures/nature_v2'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FIGURE 3 — POTENCY DISTRIBUTION")
print("=" * 60)

df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
assert len(df) == CANON['n_compounds']

# Load mmp_discontinuity_scores
disc = pd.read_csv(ROOT / 'outputs/mmp/mmp_discontinuity_scores.csv')
disc_joined = disc.merge(df[['inchi_key', 'pIC50']], on='inchi_key', how='left')
print(f"  disc scores: {len(disc)} | joined pIC50: {disc_joined['pIC50'].notna().sum()}")

fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
ax_a, ax_b, ax_c, ax_d = axes.flatten()
ax_c, ax_d = ax_d, ax_c  # FIX 4: scatter → bottom-left (c), violin → bottom-right (d)

# ── Panel a: Global histogram + KDE ──────────────────────────────────────────
print("\n[Panel a] Global histogram + KDE ...")
pic50 = df['pIC50'].dropna()
bins = np.arange(2.0, 9.25, 0.25)
n, bin_edges, patches = ax_a.hist(pic50, bins=bins, density=True,
                                    color=COLORS['blue'], alpha=0.6, edgecolor='white', lw=0.3)

x_kde = np.linspace(2.0, 9.0, 400)
kde = stats.gaussian_kde(pic50)
y_kde = kde(x_kde)
ax_a.plot(x_kde, y_kde, color=COLORS['navy'], lw=1.5)

med = pic50.median()
mean = pic50.mean()
sd = pic50.std()

ax_a.axvline(med, color=COLORS['red'], lw=1.2, ls='--', label=f'Median ({med:.2f})')
ax_a.axvline(mean, color=COLORS['dark_grey'], lw=1.0, ls=':', label=f'Mean ({mean:.2f})')

stats_txt = f'Median: {med:.2f}\nMean: {mean:.2f}\nSD: {sd:.2f}'
ax_a.text(0.97, 0.97, stats_txt, transform=ax_a.transAxes, fontsize=7.5,
          ha='right', va='top',
          bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=COLORS['grey'], lw=0.8))
ax_a.set_xlabel('pIC50', fontsize=9)
ax_a.set_ylabel('Density', fontsize=9)
ax_a.set_xlim(2.0, 9.0)
ax_a.legend(fontsize=7.5, framealpha=0.7, edgecolor='none')
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.text(-0.13, 1.04, 'a', transform=ax_a.transAxes, fontsize=11, fontweight='bold')

# ── Panel b: pIC50 by source coverage ────────────────────────────────────────
print("[Panel b] pIC50 by source coverage ...")
patent_vals = df[df['patent_flag']]['pIC50'].dropna()
multi_vals  = df[~df['patent_flag']]['pIC50'].dropna()
print(f"  Patent-only: {len(patent_vals)} | Multi-source: {len(multi_vals)}")

x_b = np.linspace(2.0, 9.0, 300)
kde_pat = stats.gaussian_kde(patent_vals)
kde_mul = stats.gaussian_kde(multi_vals)

ax_b.fill_between(x_b, kde_pat(x_b), alpha=0.3, color=COLORS['orange'])
ax_b.plot(x_b, kde_pat(x_b), lw=1.5, color=COLORS['orange'],
           label=f'Patent-exclusive (n={len(patent_vals)})')
ax_b.fill_between(x_b, kde_mul(x_b), alpha=0.3, color=COLORS['blue'])
ax_b.plot(x_b, kde_mul(x_b), lw=1.5, color=COLORS['blue'],
           label=f'Multi-source (n={len(multi_vals)})')

# Mann-Whitney U test
stat_mw, p_mw = stats.mannwhitneyu(patent_vals, multi_vals, alternative='two-sided')
p_str = f'p = {p_mw:.2e}' if p_mw >= 0.001 else f'p < 0.001'
ax_b.text(0.97, 0.97, f'Mann-Whitney U\n{p_str}',
          transform=ax_b.transAxes, fontsize=7.5, ha='right', va='top',
          bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=COLORS['grey'], lw=0.8))

ax_b.set_xlabel('pIC50', fontsize=9)
ax_b.set_ylabel('Density', fontsize=9)
ax_b.set_xlim(2.0, 9.0)
ax_b.legend(fontsize=7.5, framealpha=0.7, edgecolor='none')
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.text(-0.13, 1.04, 'b', transform=ax_b.transAxes, fontsize=11, fontweight='bold')

# ── Panel c: pIC50 by assay mechanism (violin) ────────────────────────────────
print("[Panel c] Violin by mechanism ...")
mech_order = ['enzymatic', 'enzymatic_confirmed', 'fp_ic50', 'covalent']
mech_colors = [COLORS['grey'], COLORS['teal'], COLORS['cyan'], COLORS['red']]
mech_labels = ['Enzymatic', 'Enzymatic\nconfirmed', 'FP-IC50', 'Covalent']
mech_data = [df[df['mechanism_class'] == m]['pIC50'].dropna().values for m in mech_order]
mech_ns = [len(d) for d in mech_data]
print(f"  Mechanism counts: {dict(zip(mech_order, mech_ns))}")

vp = ax_c.violinplot(mech_data, positions=range(len(mech_order)),
                      showmedians=False, showextrema=False)

for pc, color in zip(vp['bodies'], mech_colors):
    pc.set_facecolor(color)
    pc.set_alpha(0.7)
    pc.set_edgecolor(COLORS['dark_grey'])
    pc.set_linewidth(0.5)

for i, (data, color) in enumerate(zip(mech_data, mech_colors)):
    med = np.median(data)
    ax_c.hlines(med, i - 0.3, i + 0.3, color=COLORS['dark_grey'], lw=2, zorder=4)
    ax_c.text(i + 0.35, med, f'{med:.2f}', fontsize=7.5, va='center',
              color=COLORS['dark_grey'])

ax_c.set_xticks(range(len(mech_order)))
ax_c.set_xticklabels([f'{l}\n(n={n})' for l, n in zip(mech_labels, mech_ns)],
                      fontsize=7.5)
ax_c.set_ylabel('pIC50', fontsize=9)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
ax_c.text(-0.13, 1.04, 'd', transform=ax_c.transAxes, fontsize=11, fontweight='bold')

# ── Panel d: Cliff degree vs pIC50 ───────────────────────────────────────────
print("[Panel d] Cliff degree vs pIC50 ...")

HUB_IKS = {
    'A1': 'SMADULGDNOCLOP-GISFHXKWSA-N',
    'A2': 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
    'B1': 'UDCDEKJNAMHBFH-HSZRJFAPSA-N',
    'B2': 'DVCKJOQIVOGXEI-XMMPIXPASA-N',
}

hub_iks = set(HUB_IKS.values())
non_hub = disc_joined[~disc_joined['inchi_key'].isin(hub_iks)]
hub_a = disc_joined[disc_joined['inchi_key'].isin({HUB_IKS['A1'], HUB_IKS['A2']})]
hub_b = disc_joined[disc_joined['inchi_key'].isin({HUB_IKS['B1'], HUB_IKS['B2']})]

ax_d.scatter(non_hub['pIC50'], non_hub['discontinuity_score'],
              s=25, c=COLORS['blue'], alpha=0.5, zorder=2, label='Hub compounds')
ax_d.scatter(hub_a['pIC50'], hub_a['discontinuity_score'],
              s=200, c=COLORS['navy'], marker='*', zorder=5, edgecolors='white', lw=0.5)
ax_d.scatter(hub_b['pIC50'], hub_b['discontinuity_score'],
              s=180, c=COLORS['red'], marker='D', zorder=5, edgecolors='white', lw=0.5)

# Labels for hubs
for label, ik in HUB_IKS.items():
    row = disc_joined[disc_joined['inchi_key'] == ik]
    if len(row) > 0:
        x, y = row['pIC50'].values[0], row['discontinuity_score'].values[0]
        ax_d.annotate(label, (x, y), xytext=(5, 5), textcoords='offset points',
                       fontsize=8, fontweight='bold',
                       color=COLORS['navy'] if label.startswith('A') else COLORS['red'])

ax_d.set_xlabel('pIC50', fontsize=9)
ax_d.set_ylabel('Discontinuity score', fontsize=9)
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)

# Custom legend
h_nonhub = mpatches.Patch(color=COLORS['blue'], label='Non-hub')
h_a = plt.scatter([], [], s=100, c=COLORS['navy'], marker='*', label='Hub A')
h_b = plt.scatter([], [], s=80, c=COLORS['red'], marker='D', label='Hub B')
ax_d.legend(handles=[h_nonhub, h_a, h_b], fontsize=7.5, framealpha=0.7, edgecolor='none')
ax_d.text(-0.13, 1.04, 'c', transform=ax_d.transAxes, fontsize=11, fontweight='bold')

# ── Save ──────────────────────────────────────────────────────────────────────
for ext in ['png', 'pdf']:
    outpath = OUT / f'fig03_potency.{ext}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"  Saved: {outpath}")

plt.close(fig)
print("\nFigure 3 complete.")
