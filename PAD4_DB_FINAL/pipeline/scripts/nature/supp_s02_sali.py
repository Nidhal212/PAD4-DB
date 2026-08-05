"""
supp_s02_sali.py — S2: SALI Landscape (3-panel)
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
    'blue': '#0077BB', 'orange': '#EE7733', 'teal': '#009988',
    'red': '#CC3311', 'navy': '#004488', 'grey': '#BBBBBB',
    'magenta': '#EE3377', 'dark_grey': '#555555', 'light_grey': '#E8E8E8',
}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/figures/nature_v2'
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

fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(10, 5.5), constrained_layout=True)

# ── Panel a: SALI distribution histogram ────────────────────────────────────
print("[Panel a] SALI histogram ...")
sali_vals = pairs['sali'].clip(upper=pairs['sali'].quantile(0.999))  # clip extreme outliers for display
bins_a = np.linspace(0, sali_vals.max() * 1.05, 60)
ax_a.hist(sali_vals, bins=bins_a, color=COLORS['teal'], alpha=0.8, edgecolor='white', lw=0.2)
ax_a.set_yscale('log')
ax_a.axvline(10, color=COLORS['orange'], lw=1.5, ls='--', label=f'SALI>10 (n={n_sali10})')
ax_a.axvline(20, color=COLORS['red'], lw=1.5, ls='--', label=f'SALI>20 (n={n_sali20})')
ax_a.set_xlabel('SALI', fontsize=9)
ax_a.set_ylabel('Number of pairs (log)', fontsize=9)
ax_a.legend(fontsize=7.5, framealpha=0.7, edgecolor='none')
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.text(-0.15, 1.04, 'a', transform=ax_a.transAxes, fontsize=11, fontweight='bold')

# ── Panel b: SALI vs delta_pic50, colored by tanimoto ────────────────────────
print("[Panel b] SALI vs delta_pic50 ...")
# Sample for performance
pairs_sample = pairs.sample(min(50000, len(pairs)), random_state=42)

sc = ax_b.scatter(pairs_sample['delta_pic50'].abs(), pairs_sample['sali'],
                    c=pairs_sample['tanimoto'], cmap='viridis',
                    s=1, alpha=0.3, rasterized=True, zorder=2)

# Highlight SALI>20 in red
high_sali = pairs[pairs['sali'] > 20]
ax_b.scatter(high_sali['delta_pic50'].abs(), high_sali['sali'],
              s=15, c=COLORS['red'], alpha=0.9, zorder=5,
              label=f'SALI>20 (n={len(high_sali)})')

plt.colorbar(sc, ax=ax_b, label='Tanimoto similarity', shrink=0.8, pad=0.02)
ax_b.set_xlabel('|ΔpIC50|', fontsize=9)
ax_b.set_ylabel('SALI', fontsize=9)
ax_b.legend(fontsize=7.5, framealpha=0.7, edgecolor='none')
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.text(-0.15, 1.04, 'b', transform=ax_b.transAxes, fontsize=11, fontweight='bold')

# ── Panel c: Top 20 SALI pairs ────────────────────────────────────────────────
print("[Panel c] Top 20 SALI pairs ...")
top20 = pairs.nlargest(20, 'sali').reset_index(drop=True)

tier_colors = {
    'severe': COLORS['red'],
    'moderate': COLORS['magenta'],
    'broad': COLORS['blue'],
    None: COLORS['grey'],
}

bar_cols = [tier_colors.get(t, COLORS['grey']) for t in top20['cliff_tier']]
y_labels_c = [f"{row['inchi_key_a'][:8]}.../{row['inchi_key_b'][:8]}..." for _, row in top20.iterrows()]
y_pos_c = range(len(top20))

ax_c.barh(list(y_pos_c), top20['sali'].values,
           color=bar_cols, height=0.7, edgecolor='white', lw=0.3)

for i, val in enumerate(top20['sali']):
    ax_c.text(val + 0.2, i, f'{val:.1f}', va='center', ha='left', fontsize=6.5)

ax_c.set_yticks(list(y_pos_c))
ax_c.set_yticklabels(y_labels_c, fontsize=5.5, family='monospace')
ax_c.set_xlabel('SALI', fontsize=9)
ax_c.invert_yaxis()
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)

# Legend for tiers
for tier, color in tier_colors.items():
    if tier:
        ax_c.plot([], [], 's', color=color, markersize=6, label=tier.capitalize())
ax_c.legend(fontsize=7.5, framealpha=0.7, edgecolor='none', loc='lower right')
ax_c.text(-0.25, 1.04, 'c', transform=ax_c.transAxes, fontsize=11, fontweight='bold')

for ext in ['png', 'pdf']:
    outpath = OUT / f'supp_s02_sali.{ext}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"  Saved: {outpath}")

plt.close(fig)
print("S2 complete.")
