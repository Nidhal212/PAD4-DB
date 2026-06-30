"""
fig02_source_overlap.py — Figure 2: Source Overlap (2-panel)
Outputs: outputs/figures/nature/fig02_source_overlap.{png,pdf}
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

CANON = {'n_compounds': 3093, 'n_patent': 233, 'n_multi_06': 528}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/figures/nature_v2'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FIGURE 2 — SOURCE OVERLAP")
print("=" * 60)

df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
assert len(df) == CANON['n_compounds'], f"n: {len(df)}"
print(f"  {len(df)} compounds ✓")

# ── Derive source membership ───────────────────────────────────────────────────
def has_source(s, name):
    return s.str.contains(name, na=False)

df['has_pubchem'] = has_source(df['source_list'], 'pubchem')
df['has_chembl']  = has_source(df['source_list'], 'chembl')
df['has_binding'] = has_source(df['source_list'], 'bindingdb')

# Source coverage totals
n_pc = df['has_pubchem'].sum()
n_cb = df['has_chembl'].sum()
n_bd = df['has_binding'].sum()
print(f"  PubChem: {n_pc} | ChEMBL: {n_cb} | BindingDB: {n_bd}")

# Compute combinations
combos = []
for _, row in df.iterrows():
    c = []
    if row['has_pubchem']: c.append('PubChem')
    if row['has_chembl']:  c.append('ChEMBL')
    if row['has_binding']: c.append('BindingDB')
    combos.append('+'.join(sorted(c)))

df['combo'] = combos
combo_counts = df['combo'].value_counts().sort_values(ascending=False)
print(f"\n  Source combinations:\n{combo_counts}")
total = combo_counts.sum()
print(f"  Total = {total} (should be {CANON['n_compounds']})")
assert total == CANON['n_compounds'], f"Combo total mismatch: {total}"

# ── Build figure ──────────────────────────────────────────────────────────────
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10, 5.5),
                                   gridspec_kw={'width_ratios': [1.3, 1]})

# ── Panel a: UpSet-style manual plot ─────────────────────────────────────────
print("\n[Panel a] UpSet plot ...")

# Define order: largest to smallest
ordered_combos = combo_counts.index.tolist()
ordered_counts = combo_counts.values.tolist()

sources = ['PubChem', 'ChEMBL', 'BindingDB']
n_combos = len(ordered_combos)

# Top subplot: bars; bottom: matrix
gs_inner = ax_a.get_subplotspec().subgridspec(2, 1, height_ratios=[2.5, 1], hspace=0.05)
fig_gs = ax_a.get_gridspec()

# Remove original ax_a, replace with two sub-axes
ax_a.remove()
ax_bar = fig.add_subplot(gs_inner[0])
ax_mat = fig.add_subplot(gs_inner[1])

bar_colors = []
for combo in ordered_combos:
    if 'PubChem' in combo and 'ChEMBL' not in combo and 'BindingDB' not in combo:
        bar_colors.append(COLORS['orange'])
    else:
        bar_colors.append(COLORS['blue'])

bars = ax_bar.bar(range(n_combos), ordered_counts, color=bar_colors,
                   width=0.6, edgecolor='white', lw=0.5)

for i, (count, bar) in enumerate(zip(ordered_counts, bars)):
    ax_bar.text(i, count + max(ordered_counts) * 0.01, str(count),
                ha='center', va='bottom', fontsize=8, fontweight='bold',
                color=COLORS['dark_grey'])

ax_bar.set_xticks([])
ax_bar.set_ylabel('Number of compounds', fontsize=9)
ax_bar.set_xlim(-0.5, n_combos - 0.5)
ax_bar.spines['top'].set_visible(False)
ax_bar.spines['right'].set_visible(False)
ax_bar.text(-0.1, 1.04, 'a', transform=ax_bar.transAxes, fontsize=11, fontweight='bold')

# Matrix
ax_mat.set_xlim(-0.5, n_combos - 0.5)
ax_mat.set_ylim(-0.5, len(sources) - 0.5)
ax_mat.set_yticks(range(len(sources)))
ax_mat.set_yticklabels(sources[::-1], fontsize=8)
ax_mat.set_xticks([])
ax_mat.spines['top'].set_visible(False)
ax_mat.spines['right'].set_visible(False)
ax_mat.spines['bottom'].set_visible(False)

for i, combo in enumerate(ordered_combos):
    members = combo.split('+')
    for j, src in enumerate(sources[::-1]):
        filled = src in members
        ax_mat.scatter([i], [j],
                       s=80 if filled else 60,
                       c=COLORS['dark_grey'] if filled else COLORS['light_grey'],
                       zorder=3 if filled else 2,
                       edgecolors='none')
        if filled:
            # Draw connecting line between filled dots in same combo
            pass

# Draw vertical lines connecting filled dots
for i, combo in enumerate(ordered_combos):
    members = combo.split('+')
    filled_ys = [j for j, src in enumerate(sources[::-1]) if src in members]
    if len(filled_ys) > 1:
        ax_mat.plot([i, i], [min(filled_ys), max(filled_ys)],
                    color=COLORS['dark_grey'], lw=2, zorder=2)

# ── Panel b: Source coverage totals ──────────────────────────────────────────
print("[Panel b] Source coverage ...")

src_names = ['PubChem', 'BindingDB', 'ChEMBL']
src_counts = [n_pc, n_bd, n_cb]
src_colors = [COLORS['blue'], COLORS['orange'], COLORS['teal']]
total_n = CANON['n_compounds']

y_pos = range(len(src_names))
bars_b = ax_b.barh(list(y_pos), src_counts, color=src_colors,
                    height=0.5, edgecolor='white', lw=0.5)

for i, (count, bar) in enumerate(zip(src_counts, bars_b)):
    pct = count / total_n * 100
    ax_b.text(count + 30, i, f'{count:,} ({pct:.1f}%)',
              va='center', ha='left', fontsize=8, fontweight='bold',
              color=COLORS['dark_grey'])

ax_b.set_yticks(list(y_pos))
ax_b.set_yticklabels(src_names, fontsize=9)
ax_b.set_xlabel('Number of compounds', fontsize=9)
ax_b.set_xlim(0, max(src_counts) * 1.35)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

# Stats annotation
n_pc_only = combo_counts.get('PubChem', 0)
annotation = (f"233 patent-exclusive\ncompounds absent\nfrom ChEMBL and\nBindingDB")
ax_b.text(0.97, 0.25, annotation, transform=ax_b.transAxes,
          fontsize=7, ha='right', va='bottom',
          bbox=dict(boxstyle='round,pad=0.4', facecolor=COLORS['light_grey'],
                    edgecolor=COLORS['orange'], lw=1.5))
ax_b.text(-0.15, 1.04, 'b', transform=ax_b.transAxes, fontsize=11, fontweight='bold')

# ── Save ──────────────────────────────────────────────────────────────────────
for ext in ['png', 'pdf']:
    outpath = OUT / f'fig02_source_overlap.{ext}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"  Saved: {outpath}")

plt.close(fig)
print("\nFigure 2 complete.")
