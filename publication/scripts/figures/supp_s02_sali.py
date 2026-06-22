"""
supp_s02_sali.py — S2: SALI Landscape (3-panel)
Outputs: publication/figures/supplementary/fig_s02_sali.{png,pdf}
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

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'publication/figures/supplementary'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("SUPPLEMENTARY S2 — SALI LANDSCAPE")
print("=" * 60)

pairs = pd.read_parquet(ROOT / 'data/processed/activity_pairs_with_sali.parquet')
print(f"  Pairs: {len(pairs)}")
print(f"  SALI range: {pairs['sali'].min():.2f} – {pairs['sali'].max():.2f}")

n_sali10 = (pairs['sali'] > 10).sum()
n_sali20 = (pairs['sali'] > 20).sum()
print(f"  SALI>10: {n_sali10} | SALI>20: {n_sali20}")

fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(DOUBLE, 3.0), constrained_layout=True)

# ── Panel a: SALI distribution histogram ────────────────────────────────────
print("[Panel a] SALI histogram ...")
sali_vals = pairs['sali'].clip(upper=pairs['sali'].quantile(0.999))
bins_a = np.linspace(0, sali_vals.max() * 1.05, 60)
ax_a.hist(sali_vals, bins=bins_a, color=SEM['enzymatic_confirmed'], alpha=0.8,
           edgecolor='white', lw=0.2)
ax_a.set_yscale('log')
ax_a.axvline(10, color=SEM['patent'], lw=1.2, ls='--', label=f'SALI>10 (n={n_sali10})')
ax_a.axvline(20, color=SEM['cliff'], lw=1.2, ls='--', label=f'SALI>20 (n={n_sali20})')
ax_a.set_xlabel('SALI')
ax_a.set_ylabel('Pairs (log scale)')
ax_a.legend(fontsize=5.5, framealpha=0.7, edgecolor='none')
panel_label(ax_a, 'a')

# ── Panel b: SALI vs delta_pic50, colored by tanimoto ────────────────────────
print("[Panel b] SALI vs delta_pic50 ...")
pairs_sample = pairs.sample(min(50000, len(pairs)), random_state=42)

sc = ax_b.scatter(pairs_sample['delta_pic50'].abs(), pairs_sample['sali'],
                    c=pairs_sample['tanimoto'], cmap='viridis',
                    s=0.5, alpha=0.25, rasterized=True, zorder=2)

high_sali = pairs[pairs['sali'] > 20]
ax_b.scatter(high_sali['delta_pic50'].abs(), high_sali['sali'],
              s=8, c=SEM['cliff'], alpha=0.9, zorder=5,
              label=f'SALI>20 (n={len(high_sali)})')

cbar = plt.colorbar(sc, ax=ax_b, shrink=0.8, pad=0.02)
cbar.set_label('Tanimoto', labelpad=8)
ax_b.set_xlabel('|ΔpIC50|')
ax_b.set_ylabel('SALI')
ax_b.legend(fontsize=5.5, framealpha=0.7, edgecolor='none')
panel_label(ax_b, 'b')

# ── Panel c: Top 20 SALI pairs ────────────────────────────────────────────────
print("[Panel c] Top 20 SALI pairs ...")
top20 = pairs.nlargest(20, 'sali').reset_index(drop=True)

tier_colors = {
    'severe':   SEM['cliff'],
    'moderate': SEM['patent'],
    'broad':    SEM['published'],
    None:       C['grey'],
}

bar_cols = [tier_colors.get(t, C['grey']) for t in top20['cliff_tier']]
y_pos_c = range(len(top20))

ax_c.barh(list(y_pos_c), top20['sali'].values,
           color=bar_cols, height=0.7, edgecolor='white', lw=0.3)

for i, val in enumerate(top20['sali']):
    ax_c.text(val + 0.2, i, f'{val:.1f}', va='center', ha='left', fontsize=5.5, clip_on=False)

# Use Pair N labels — InChIKey mapping available in deposited dataset
ax_c.set_yticks(list(y_pos_c))
ax_c.set_yticklabels([f'Pair {i+1}' for i in range(len(top20))], fontsize=6)
ax_c.set_xlabel('SALI')
ax_c.invert_yaxis()

for tier, color in tier_colors.items():
    if tier:
        ax_c.plot([], [], 's', color=color, markersize=5, label=tier.capitalize())
ax_c.legend(fontsize=5.5, framealpha=0.7, edgecolor='none', loc='lower right')
panel_label(ax_c, 'c', x=-0.3)

save_fig(fig, str(OUT / 'fig_s02_sali'))
plt.close(fig)
print("S2 complete.")
