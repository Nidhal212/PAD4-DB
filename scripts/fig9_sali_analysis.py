#!/usr/bin/env python
"""Figure 9 — SALI Analysis (supplementary) + Great Tables."""
import os, sys, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import scienceplots
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
plt.style.use(['science', 'nature', 'no-latex'])

from great_tables import GT, loc, style as gt_style
import great_tables
import numpy as np
import pandas as pd

os.makedirs('outputs/figures', exist_ok=True)
os.makedirs('outputs/tables', exist_ok=True)

print(f"SciencePlots: importable  |  Great Tables: {great_tables.__version__}")

# Load data
pairs = pd.read_parquet('data/processed/activity_pairs_with_sali.parquet')

# Locked number verification
n_sali_gt10  = int((pairs['sali'] > 10).sum())
n_sali_gt20  = int((pairs['sali'] > 20).sum())
sali_max_val = float(pairs['sali'].max())

LOCKED = [
    ("SALI > 10",           n_sali_gt10,              335),
    ("SALI > 20",           n_sali_gt20,              19),
    ("SALI max (round 2)",  round(sali_max_val, 2),   65.88),
]
all_pass = True
for label, actual, expected in LOCKED:
    ok = actual == expected
    print(f"  {label}: {actual}  {'PASS' if ok else f'FAIL (expected {expected})'}")
    if not ok:
        all_pass = False
if not all_pass:
    print("\nVERIFICATION FAILED — stopping.")
    sys.exit(1)
print()

# Subsets
valid   = pairs[pairs['sali'].notna()]
top20   = pairs.nlargest(20, 'sali').reset_index(drop=True)
extreme = pairs[pairs['sali'] > 20]    # 19 pairs

TIER_COLORS = {
    'severe':    '#E74C3C',
    'moderate':  '#F39C12',
    'broad':     '#4A90D9',
    'non_cliff': '#AAAAAA',
}

fig, axes = plt.subplots(1, 3, figsize=(14, 5))
ax_a, ax_b, ax_c = axes


def panel_label(ax, letter):
    ax.text(0.02, 0.96, letter, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top')


# Panel A: Top 20 pairs by SALI (horizontal bar)
top20['pair_label'] = (top20['inchi_key_a'].str[:10] + '.../' +
                       top20['inchi_key_b'].str[:10] + '...')
bar_colors = [TIER_COLORS.get(t, '#AAAAAA') for t in top20['cliff_tier']]
y_pos = np.arange(len(top20))

ax_a.barh(y_pos, top20['sali'], color=bar_colors, height=0.7,
          edgecolor='white', linewidth=0.3)
ax_a.set_yticks(y_pos)
ax_a.set_yticklabels(top20['pair_label'], fontsize=5.5)
ax_a.invert_yaxis()
ax_a.set_xlabel('SALI', fontsize=9)
ax_a.set_title('Top 20 Pairs by SALI', fontsize=9)

for i, (_, row) in enumerate(top20.iterrows()):
    ax_a.text(row['sali'] + 0.4, i, f"{row['sali']:.1f}",
              va='center', fontsize=5.5, color='#333333')
ax_a.set_xlim(0, sali_max_val + 10)

ax_a.axvline(sali_max_val, color='#E74C3C', linestyle='--',
             linewidth=0.8, alpha=0.6)
ax_a.text(sali_max_val - 1, 19.2, f'Max={sali_max_val:.2f}',
          ha='right', va='bottom', fontsize=6.5, color='#E74C3C')

legend_patches = [
    mpatches.Patch(color=TIER_COLORS['severe'],   label='Severe'),
    mpatches.Patch(color=TIER_COLORS['moderate'], label='Moderate'),
    mpatches.Patch(color=TIER_COLORS['broad'],    label='Broad'),
]
ax_a.legend(handles=legend_patches, loc='lower right', fontsize=7, framealpha=0.8)
panel_label(ax_a, 'A')

# Panel B: SALI vs delta_pic50 colored by Tanimoto
high_sali = valid[valid['sali'] > 1]
rng = np.random.default_rng(42)
n_samp = min(20000, len(high_sali))
idx = rng.choice(len(high_sali), size=n_samp, replace=False)
samp = high_sali.iloc[idx]

sc = ax_b.scatter(samp['delta_pic50'], samp['sali'],
                  c=samp['tanimoto'], cmap='viridis_r',
                  s=3, alpha=0.4, rasterized=True,
                  vmin=0.6, vmax=1.0)
ax_b.scatter(extreme['delta_pic50'], extreme['sali'],
             c='#E74C3C', s=30, zorder=5,
             edgecolors='white', linewidths=0.3,
             label=f'SALI>20 (n={len(extreme)})')

cb = plt.colorbar(sc, ax=ax_b, pad=0.02)
cb.set_label('Tanimoto Similarity', fontsize=8)
cb.ax.tick_params(labelsize=7)

ax_b.set_yscale('log')
ax_b.set_xlabel('ΔpIC50', fontsize=9)
ax_b.set_ylabel('SALI (log scale)', fontsize=9)
ax_b.set_title('SALI vs ΔpIC50', fontsize=9)
ax_b.legend(fontsize=7, framealpha=0.8, loc='upper left')
panel_label(ax_b, 'B')

# Panel C: Cumulative SALI distribution (ECDF)
all_sali = np.sort(valid['sali'].values)
cumfrac  = np.arange(1, len(all_sali) + 1) / len(all_sali)

mask_plot = all_sali <= 70
ax_c.plot(all_sali[mask_plot], cumfrac[mask_plot],
          color='#4A90D9', linewidth=1.5)

frac10 = float(np.interp(10.0, all_sali, cumfrac))
frac20 = float(np.interp(20.0, all_sali, cumfrac))

ax_c.axvline(10, color='#F39C12', linestyle='--', linewidth=1.0, label='SALI=10')
ax_c.axvline(20, color='#E74C3C', linestyle='--', linewidth=1.0, label='SALI=20')

ax_c.text(11.5, 0.88, f"335 pairs\n(SALI>10)",
          fontsize=7.5, color='#F39C12', va='top')
ax_c.text(21.5, 0.75, f"19 pairs\n(SALI>20)",
          fontsize=7.5, color='#E74C3C', va='top')

ax_c.set_xlim(0, 70)
ax_c.set_ylim(0, 1.02)
ax_c.set_xlabel('SALI', fontsize=9)
ax_c.set_ylabel('Cumulative Fraction', fontsize=9)
ax_c.set_title('Cumulative SALI Distribution', fontsize=9)
ax_c.legend(fontsize=7, framealpha=0.8)
panel_label(ax_c, 'C')

plt.tight_layout(pad=0.8)

PNG_PATH = 'outputs/figures/fig9_sali_analysis.png'
SVG_PATH = 'outputs/figures/fig9_sali_analysis.svg'
fig.savefig(PNG_PATH, dpi=300, bbox_inches='tight')
fig.savefig(SVG_PATH, bbox_inches='tight')
plt.close()
print(f"Saved: {PNG_PATH}  ({os.path.getsize(PNG_PATH)/1024:.1f} KB)")
print(f"Saved: {SVG_PATH}  ({os.path.getsize(SVG_PATH)/1024:.1f} KB)")
print()

# Great Tables: top 20 pairs
gt_df = top20[['inchi_key_a', 'inchi_key_b', 'tanimoto',
               'delta_pic50', 'sali', 'cliff_tier']].copy()
gt_df.insert(0, 'rank', range(1, 21))
gt_df['ik_a'] = gt_df['inchi_key_a'].str[:14]
gt_df['ik_b'] = gt_df['inchi_key_b'].str[:14]
gt_out = gt_df[['rank', 'ik_a', 'ik_b', 'tanimoto',
                'delta_pic50', 'sali', 'cliff_tier']]

gt = (
    GT(gt_out)
    .tab_header(
        title="PAD4-DB Top 20 Pairs by SALI",
        subtitle="SALI = |ΔpIC50| / (1 − Tanimoto); pairs with Tanimoto=1.0 excluded",
    )
    .cols_label(
        rank="Rank", ik_a="InChIKey A (14)", ik_b="InChIKey B (14)",
        tanimoto="Tanimoto", delta_pic50="ΔpIC50",
        sali="SALI", cliff_tier="Tier",
    )
    .fmt_number(columns=['tanimoto', 'delta_pic50', 'sali'], decimals=2)
    .tab_style(
        style=gt_style.fill(color="#FFEBEE"),
        locations=loc.body(rows=[0]),
    )
    .tab_source_note(
        "335 pairs with SALI>10; 19 pairs with SALI>20; "
        "max SALI=65.88 (GBABZCBYODXOKT vs nearest neighbor)."
    )
)

HTML_PATH = 'outputs/tables/fig9_sali_top_pairs.html'
with open(HTML_PATH, 'w') as f:
    f.write(gt.as_raw_html())
print(f"Great Tables: {HTML_PATH}  ({os.path.getsize(HTML_PATH)/1024:.1f} KB)")
print("TASK B: DONE")
