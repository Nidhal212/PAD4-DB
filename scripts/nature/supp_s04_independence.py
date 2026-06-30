"""
supp_s04_independence.py — S4: Source Independence (3-panel)
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
    'blue': '#0077BB', 'orange': '#EE7733', 'teal': '#009988',
    'navy': '#004488', 'grey': '#BBBBBB', 'dark_grey': '#555555',
    'light_grey': '#E8E8E8', 'red': '#CC3311',
}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/figures/nature_v2'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("SUPPLEMENTARY S4 — SOURCE INDEPENDENCE")
print("=" * 60)

df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
sis = df['source_independence_score'].dropna()
print(f"  SIS range: {sis.min():.3f} – {sis.max():.3f}")
print(f"  n>=0.6: {(sis>=0.6).sum()} | n>=0.7: {(sis>=0.7).sum()}")

fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(12, 4.5), constrained_layout=True)

# ── Panel a: Lollipop of SIS distribution ─────────────────────────────────────
print("[Panel a] SIS lollipop ...")
sis_rounded = sis.round(2)
sis_counts = sis_rounded.value_counts().sort_index()

ax_a.vlines(sis_counts.index, 0, sis_counts.values,
             color=COLORS['blue'], lw=1.0, alpha=0.6)
ax_a.scatter(sis_counts.index, sis_counts.values,
              s=20, c=COLORS['blue'], zorder=3, alpha=0.8)
ax_a.axvline(0.6, color=COLORS['orange'], lw=1.5, ls='--', label='Threshold 0.6')
ax_a.axvline(0.7, color=COLORS['red'], lw=1.5, ls=':', label='Threshold 0.7')
ax_a.set_xlabel('Source independence score', fontsize=9)
ax_a.set_ylabel('Number of compounds', fontsize=9)
ax_a.legend(fontsize=7.5, framealpha=0.7, edgecolor='none')
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.text(-0.15, 1.04, 'a', transform=ax_a.transAxes, fontsize=11, fontweight='bold')

# ── Panel b: Threshold comparison bar chart ───────────────────────────────────
print("[Panel b] Threshold bar chart ...")
n_06 = (sis >= 0.6).sum()
n_07 = (sis >= 0.7).sum()
labels_b = [f'Score ≥ 0.6\n(n={n_06})', f'Score ≥ 0.7\n(n={n_07})']
vals_b = [n_06, n_07]
colors_b = [COLORS['orange'], COLORS['red']]

bars_b = ax_b.bar([0, 1], vals_b, color=colors_b, alpha=0.8,
                   width=0.5, edgecolor='white', lw=0.3)
for bar, val in zip(bars_b, vals_b):
    ax_b.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5,
               str(val), ha='center', va='bottom', fontsize=9, fontweight='bold',
               color=COLORS['dark_grey'])

ax_b.set_xticks([0, 1])
ax_b.set_xticklabels(labels_b, fontsize=8.5)
ax_b.set_ylabel('Number of compounds', fontsize=9)
ax_b.set_ylim(0, max(vals_b) * 1.15)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.text(-0.18, 1.04, 'b', transform=ax_b.transAxes, fontsize=11, fontweight='bold')

# ── Panel c: KDE pIC50 multi vs single source ────────────────────────────────
print("[Panel c] pIC50 KDE by independence ...")
multi_vals = df[df['source_independence_score'] >= 0.6]['pIC50'].dropna()
single_vals = df[df['source_independence_score'] < 0.6]['pIC50'].dropna()
print(f"  Multi (>=0.6): {len(multi_vals)} | Single (<0.6): {len(single_vals)}")

x_c = np.linspace(2.0, 9.0, 300)
kde_multi = stats.gaussian_kde(multi_vals)
kde_single = stats.gaussian_kde(single_vals)

ax_c.fill_between(x_c, kde_multi(x_c), alpha=0.3, color=COLORS['blue'])
ax_c.plot(x_c, kde_multi(x_c), lw=1.5, color=COLORS['blue'],
           label=f'Score ≥ 0.6 (n={len(multi_vals)})')
ax_c.fill_between(x_c, kde_single(x_c), alpha=0.3, color=COLORS['grey'])
ax_c.plot(x_c, kde_single(x_c), lw=1.5, color=COLORS['dark_grey'],
           label=f'Score < 0.6 (n={len(single_vals)})')

# Mean lines
ax_c.axvline(multi_vals.mean(), color=COLORS['blue'], lw=1.0, ls='--', alpha=0.8)
ax_c.axvline(single_vals.mean(), color=COLORS['dark_grey'], lw=1.0, ls='--', alpha=0.8)

ax_c.set_xlabel('pIC50', fontsize=9)
ax_c.set_ylabel('Density', fontsize=9)
ax_c.set_xlim(2.0, 9.0)
ax_c.legend(fontsize=7.5, framealpha=0.7, edgecolor='none')
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
ax_c.text(-0.15, 1.04, 'c', transform=ax_c.transAxes, fontsize=11, fontweight='bold')

for ext in ['png', 'pdf']:
    outpath = OUT / f'supp_s04_independence.{ext}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"  Saved: {outpath}")

plt.close(fig)
print("S4 complete.")
