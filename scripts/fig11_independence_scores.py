#!/usr/bin/env python
"""Figure 11 — Source Independence Scores (supplementary) + Great Tables."""
import os, sys, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import scienceplots
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
plt.style.use(['science', 'nature', 'no-latex'])

from great_tables import GT, loc, style as gt_style
import great_tables
import numpy as np
import pandas as pd

os.makedirs('outputs/figures', exist_ok=True)
os.makedirs('outputs/tables', exist_ok=True)

print(f"SciencePlots: importable  |  Great Tables: {great_tables.__version__}")

# Load data
df = pd.read_parquet('data/processed/pad4_compounds.parquet')
assert len(df) == 3093

# Derive true_multi at threshold 0.6 (locked: 528 / 2565)
df['_true_multi'] = df['source_independence_score'] >= 0.6

n_true  = int(df['_true_multi'].sum())
n_false = int((~df['_true_multi']).sum())

print("=== Locked number verification ===")
LOCKED = [
    ("is_true_multi_source=True  (score≥0.6)", n_true,  528),
    ("is_true_multi_source=False (score<0.6)", n_false, 2565),
    ("Total compounds",                         len(df), 3093),
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

# Source combination labels (shortened)
LABEL_MAP = {
    'bindingdb|chembl|pubchem_confirmatory': 'BDB+ChEMBL+PC',
    'bindingdb|pubchem_confirmatory':        'BDB+PC',
    'pubchem_confirmatory':                  'PC only',
    'bindingdb|chembl':                      'BDB+ChEMBL',
    'bindingdb':                             'BDB only',
    'chembl|pubchem_confirmatory':           'ChEMBL+PC',
    'chembl':                                'ChEMBL only',
}
df['source_short'] = df['source_list'].map(LABEL_MAP).fillna(df['source_list'])

# Per-source-combination stats (score is deterministic per combination)
src_stats = (
    df.groupby('source_short')
    .agg(score=('source_independence_score', 'first'),
         n=('inchi_key', 'count'))
    .reset_index()
    .sort_values('score', ascending=False)
)
print("Source combination scores:")
print(src_stats.to_string(index=False))
print()

# Split groups for Panel C
true_multi  = df[df['_true_multi']]['pic50_consensus'].dropna().values
false_multi = df[~df['_true_multi']]['pic50_consensus'].dropna().values

# Figure
fig, axes = plt.subplots(1, 3, figsize=(13, 5))
ax_a, ax_b, ax_c = axes


def panel_label(ax, letter):
    ax.text(0.02, 0.96, letter, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top')


# Panel A: Score distribution histogram
scores = df['source_independence_score'].values
ax_a.hist(scores, bins=40, range=(0, 1),
          color='#4A90D9', alpha=0.7, edgecolor='white', linewidth=0.3)
ax_a.axvline(0.6, color='#E74C3C', linestyle='--', linewidth=1.2,
             label='Threshold (0.6)')
ax_a.axvspan(0.6, 1.0, alpha=0.08, color='#E74C3C', zorder=0)

ymax = ax_a.get_ylim()[1]
ax_a.text(0.80, ymax * 0.65,
          f"528 true\nmulti-source",
          ha='center', fontsize=8, color='#E74C3C',
          multialignment='center')
ax_a.text(0.32, ymax * 0.65,
          f"2,565 pipeline\nredundancy",
          ha='center', fontsize=8, color='#555555',
          multialignment='center')

ax_a.set_xlabel('Source Independence Score', fontsize=9)
ax_a.set_ylabel('Count', fontsize=9)
ax_a.set_title('Score Distribution', fontsize=9)
ax_a.legend(fontsize=7.5, framealpha=0.8)
panel_label(ax_a, 'A')

# Panel B: Score by source combination (bar chart — scores are deterministic)
bar_colors = ['#E74C3C' if s >= 0.6 else '#AAAAAA'
              for s in src_stats['score']]
y_pos = np.arange(len(src_stats))
bars = ax_b.barh(y_pos, src_stats['score'], color=bar_colors,
                 height=0.65, edgecolor='white', linewidth=0.3)
ax_b.set_yticks(y_pos)
ax_b.set_yticklabels(src_stats['source_short'], fontsize=8)
ax_b.invert_yaxis()
ax_b.axvline(0.6, color='#E74C3C', linestyle='--',
             linewidth=0.9, alpha=0.7)
ax_b.set_xlabel('Source Independence Score', fontsize=9)
ax_b.set_title('Score by Source Combination', fontsize=9)
ax_b.set_xlim(0, 1.15)

# Annotate n= counts on bars
for i, (_, row) in enumerate(src_stats.iterrows()):
    ax_b.text(row['score'] + 0.02, i,
              f"n={row['n']:,}", va='center', fontsize=7, color='#333333')
panel_label(ax_b, 'B')

# Panel C: pIC50 KDE — true multi-source vs pipeline redundancy
x_grid = np.linspace(1.5, 9.0, 400)
kde_true  = gaussian_kde(true_multi,  bw_method=0.25)
kde_false = gaussian_kde(false_multi, bw_method=0.25)
y_true  = kde_true(x_grid)
y_false = kde_false(x_grid)

ax_c.plot(x_grid, y_true,  color='#E74C3C', linewidth=1.5,
          label=f'True multi-source (n={n_true:,}, score≥0.6)')
ax_c.plot(x_grid, y_false, color='#4A90D9', linewidth=1.5,
          label=f'Pipeline redundancy (n={n_false:,}, score<0.6)')

mean_true  = float(np.mean(true_multi))
mean_false = float(np.mean(false_multi))
ax_c.axvline(mean_true,  color='#E74C3C', linestyle='--', linewidth=0.8, alpha=0.7)
ax_c.axvline(mean_false, color='#4A90D9', linestyle='--', linewidth=0.8, alpha=0.7)
ax_c.text(mean_true + 0.05,  max(y_true) * 0.85,
          f'{mean_true:.2f}',  fontsize=7, color='#E74C3C')
ax_c.text(mean_false - 0.08, max(y_false) * 0.85,
          f'{mean_false:.2f}', fontsize=7, color='#4A90D9', ha='right')

ax_c.set_xlabel('pIC50', fontsize=9)
ax_c.set_ylabel('Density', fontsize=9)
ax_c.set_title('pIC50: True Multi-source vs Pipeline Redundancy', fontsize=9)
ax_c.legend(fontsize=7, framealpha=0.8)
panel_label(ax_c, 'C')

plt.tight_layout(pad=0.8)

PNG_PATH = 'outputs/figures/fig11_independence_scores.png'
SVG_PATH = 'outputs/figures/fig11_independence_scores.svg'
fig.savefig(PNG_PATH, dpi=300, bbox_inches='tight')
fig.savefig(SVG_PATH, bbox_inches='tight')
plt.close()
print(f"Saved: {PNG_PATH}  ({os.path.getsize(PNG_PATH)/1024:.1f} KB)")
print(f"Saved: {SVG_PATH}  ({os.path.getsize(SVG_PATH)/1024:.1f} KB)")
print()

# Great Tables: independence summary
gt_rows = src_stats.rename(columns={
    'source_short': 'source_combination',
    'score':        'independence_score',
    'n':            'n_compounds',
}).copy()
gt_rows['true_multi'] = gt_rows['independence_score'] >= 0.6
gt_rows['pct_dataset'] = gt_rows['n_compounds'] / 3093 * 100

gt = (
    GT(gt_rows[['source_combination', 'independence_score',
                'n_compounds', 'pct_dataset', 'true_multi']])
    .tab_header(
        title="PAD4-DB Source Independence Score by Combination",
        subtitle=f"3,093 compounds; threshold ≥ 0.6 → 528 true multi-source",
    )
    .cols_label(
        source_combination="Source Combination",
        independence_score="Score",
        n_compounds="N",
        pct_dataset="% Dataset",
        true_multi="True Multi-source",
    )
    .fmt_number(columns=['independence_score', 'pct_dataset'], decimals=1)
    .fmt_integer(columns=['n_compounds'])
    .tab_style(
        style=gt_style.fill(color="#FFEBEE"),
        locations=loc.body(rows=[r for r, s in enumerate(gt_rows['independence_score'])
                                 if s >= 0.6]),
    )
    .tab_source_note(
        "Score reflects degree of independent experimental replication across sources. "
        "Score < 0.6 indicates pipeline overlap artifact (shared assay campaigns). "
        "Score ≥ 0.6 indicates genuine independent measurement."
    )
)

HTML_PATH = 'outputs/tables/fig11_independence_stats.html'
with open(HTML_PATH, 'w') as f:
    f.write(gt.as_raw_html())
print(f"Great Tables: {HTML_PATH}  ({os.path.getsize(HTML_PATH)/1024:.1f} KB)")

print()
print("=== Completion Report ===")
print(f"TASK D (Fig 11 independence): DONE")
print("TASK D: DONE")
