"""
supp_s04_independence.py — S4: Source Independence (3-panel)
Outputs: publication/figures/supplementary/fig_s04_independence.{png,pdf}
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, panel_label, save_fig, C, DOUBLE

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path

set_style()

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT  = ROOT / 'publication/figures/supplementary'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("SUPPLEMENTARY S4 — SOURCE INDEPENDENCE")
print("=" * 60)

df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
sis = df['source_independence_score'].dropna()
print(f"  SIS range: {sis.min():.3f} – {sis.max():.3f}")
print(f"  n>=0.6: {(sis>=0.6).sum()} | n>=0.7: {(sis>=0.7).sum()}")

fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(DOUBLE, 2.6),
                                         constrained_layout=True)

# ── Panel a: Lollipop of SIS distribution ────────────────────────────────────
print("[Panel a] SIS lollipop ...")
sis_counts = sis.round(2).value_counts().sort_index()

ax_a.vlines(sis_counts.index, 0, sis_counts.values,
            color=C['blue'], lw=1.2, alpha=0.7)
ax_a.scatter(sis_counts.index, sis_counts.values,
             s=22, c=C['blue'], zorder=3)
ax_a.axvline(0.6, color=C['orange'], lw=1.0, ls='--', label='Threshold 0.6')
ax_a.axvline(0.7, color=C['red'],    lw=1.0, ls=':',  label='Threshold 0.7')
ax_a.set_xlabel('Source independence score')
ax_a.set_ylabel('Compounds')
ax_a.legend(fontsize=6, framealpha=0.7, edgecolor='none')
panel_label(ax_a, 'a')

# ── Panel b: Threshold comparison bar chart ───────────────────────────────────
print("[Panel b] Threshold bar chart ...")
n_06 = int((sis >= 0.6).sum())
n_07 = int((sis >= 0.7).sum())

bars_b = ax_b.bar([0, 1], [n_06, n_07],
                   color=[C['orange'], C['red']],
                   alpha=0.85, width=0.45,
                   edgecolor='white', lw=0.3)
for bar, val in zip(bars_b, [n_06, n_07]):
    ax_b.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 4,
              str(val), ha='center', va='bottom', fontsize=6.5, fontweight='bold')

ax_b.set_xticks([0, 1])
ax_b.set_xticklabels([f'Score ≥ 0.6\n(n={n_06})', f'Score ≥ 0.7\n(n={n_07})'],
                      fontsize=6)
ax_b.set_ylabel('Compounds')
ax_b.set_ylim(0, max(n_06, n_07) * 1.18)
ax_b.set_title('Non-redundant at each threshold', fontsize=6, pad=3)
panel_label(ax_b, 'b')

# ── Panel c: KDE pIC50 by independence tier ───────────────────────────────────
print("[Panel c] pIC50 KDE ...")
multi_vals  = df[df['source_independence_score'] >= 0.6]['pIC50'].dropna()
single_vals = df[df['source_independence_score'] <  0.6]['pIC50'].dropna()

x_c = np.linspace(2.0, 9.0, 300)
ax_c.fill_between(x_c, stats.gaussian_kde(multi_vals)(x_c),
                  alpha=0.25, color=C['blue'])
ax_c.plot(x_c, stats.gaussian_kde(multi_vals)(x_c), lw=1.3, color=C['blue'],
          label=f'Non-redundant ≥0.6 (n={len(multi_vals)})')
ax_c.fill_between(x_c, stats.gaussian_kde(single_vals)(x_c),
                  alpha=0.25, color=C['grey'])
ax_c.plot(x_c, stats.gaussian_kde(single_vals)(x_c), lw=1.3, color=C['gray_dark'],
          label=f'Pipeline-redundant <0.6 (n={len(single_vals)})')

ax_c.axvline(multi_vals.mean(),  color=C['blue'],      lw=0.8, ls='--', alpha=0.7)
ax_c.axvline(single_vals.mean(), color=C['gray_dark'], lw=0.8, ls='--', alpha=0.7)

ax_c.set_xlabel('pIC50')
ax_c.set_ylabel('Density')
ax_c.set_xlim(2.0, 9.0)
ax_c.legend(fontsize=6, framealpha=0.7, edgecolor='none')
panel_label(ax_c, 'c')

save_fig(fig, str(OUT / 'fig_s04_independence'))
plt.close(fig)
print("S4 complete.")
