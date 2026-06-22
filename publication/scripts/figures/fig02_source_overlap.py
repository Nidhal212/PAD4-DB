"""
fig02_source_overlap.py — Figure 2: Source Overlap + Independence Score (3-panel, DOUBLE width)

Panel a: UpSet plot
Panel b: Horizontal bars per source (with legend mapping colors to source names)
Panel c: Stacked bar — 528 non-redundant vs 2,565 redundant; 17.1% annotated

Outputs: publication/figures/main/fig02_source_overlap.{png,pdf}
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
from pathlib import Path

set_style()

CANON = {'n_compounds': 3093, 'n_nonredundant': 528, 'n_redundant': 2565}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT  = ROOT / 'publication/figures/main'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FIGURE 2 — SOURCE OVERLAP + INDEPENDENCE SCORE")
print("=" * 60)

df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
assert len(df) == CANON['n_compounds'], f"n={len(df)}"

# ── Derive source membership ───────────────────────────────────────────────────
df['has_pubchem'] = df['source_list'].str.contains('pubchem', na=False)
df['has_chembl']  = df['source_list'].str.contains('chembl',  na=False)
df['has_binding'] = df['source_list'].str.contains('bindingdb', na=False)

n_pc = int(df['has_pubchem'].sum())
n_cb = int(df['has_chembl'].sum())
n_bd = int(df['has_binding'].sum())

combos = []
for _, row in df.iterrows():
    c = []
    if row['has_pubchem']: c.append('PubChem')
    if row['has_chembl']:  c.append('ChEMBL')
    if row['has_binding']: c.append('BindingDB')
    combos.append('+'.join(sorted(c)))

df['combo'] = combos
combo_counts = df['combo'].value_counts().sort_values(ascending=False)
assert combo_counts.sum() == CANON['n_compounds']

# Independence score counts
n_nonred = int((df['source_independence_score'] >= 0.6).sum())
n_red    = CANON['n_compounds'] - n_nonred
assert n_nonred == CANON['n_nonredundant'], f"non-redundant: {n_nonred}"
assert n_red    == CANON['n_redundant'],    f"redundant: {n_red}"

# ── Layout ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(DOUBLE, 3.8), constrained_layout=True)
gs  = fig.add_gridspec(2, 3,
                        height_ratios=[2.4, 1.0],
                        width_ratios=[2.2, 1.0, 1.0])

ax_bar = fig.add_subplot(gs[0, 0])
ax_mat = fig.add_subplot(gs[1, 0])
ax_b   = fig.add_subplot(gs[:, 1])
ax_c   = fig.add_subplot(gs[:, 2])

# ── Panel a: UpSet-style plot ─────────────────────────────────────────────────
print("[Panel a] UpSet plot ...")
ordered_combos  = combo_counts.index.tolist()
ordered_counts  = combo_counts.values.tolist()
sources         = ['PubChem', 'ChEMBL', 'BindingDB']
n_combos        = len(ordered_combos)

def combo_color(combo):
    n = len(combo.split('+'))
    if n == 3: return C['teal']
    if n == 2: return C['blue']
    return C['orange']

bar_colors = [combo_color(c) for c in ordered_combos]
bars = ax_bar.bar(range(n_combos), ordered_counts,
                   color=bar_colors, width=0.6, edgecolor='white', lw=0.4)

for i, (count, bar) in enumerate(zip(ordered_counts, bars)):
    ax_bar.text(i, count + max(ordered_counts) * 0.01, str(count),
                ha='center', va='bottom', fontsize=5.5, fontweight='bold',
                color=C['black'])

ax_bar.set_xticks([])
ax_bar.set_ylabel('Compounds')
ax_bar.set_xlim(-0.5, n_combos - 0.5)
ax_bar.legend(handles=[
    mpatches.Patch(color=C['teal'],   label='3 sources'),
    mpatches.Patch(color=C['blue'],   label='2 sources'),
    mpatches.Patch(color=C['orange'], label='1 source'),
], fontsize=5.5, loc='upper right', framealpha=0.85, edgecolor='none')
panel_label(ax_bar, 'a', x=-0.10, y=1.04)

# Matrix
ax_mat.set_xlim(-0.5, n_combos - 0.5)
ax_mat.set_ylim(-0.5, len(sources) - 0.5)
ax_mat.set_yticks(range(len(sources)))
ax_mat.set_yticklabels(sources[::-1], fontsize=6)
ax_mat.set_xticks([])
ax_mat.spines['top'].set_visible(False)
ax_mat.spines['right'].set_visible(False)
ax_mat.spines['bottom'].set_visible(False)

for i, combo in enumerate(ordered_combos):
    members   = combo.split('+')
    col       = combo_color(combo)
    filled_ys = []
    for j, src in enumerate(sources[::-1]):
        filled = src in members
        ax_mat.scatter([i], [j],
                       s=35 if filled else 18,
                       c=col if filled else '#DDDDDD',
                       zorder=3 if filled else 2,
                       edgecolors='none')
        if filled:
            filled_ys.append(j)
    if len(filled_ys) > 1:
        ax_mat.plot([i, i], [min(filled_ys), max(filled_ys)],
                    color=col, lw=1.8, zorder=2)

# ── Panel b: Source coverage totals (with legend) ─────────────────────────────
print("[Panel b] Source coverage bars + legend ...")
src_names  = ['PubChem', 'BindingDB', 'ChEMBL']
src_counts = [n_pc, n_bd, n_cb]
# Colorblind-safe: PubChem=blue, BindingDB=teal, ChEMBL=navy
src_colors = [C['blue'], C['teal'], C['navy']]

y_pos  = list(range(len(src_names)))
bars_b = ax_b.barh(y_pos, src_counts,
                    color=src_colors, height=0.5,
                    edgecolor='white', lw=0.4)

for i, (count, color) in enumerate(zip(src_counts, src_colors)):
    pct = count / CANON['n_compounds'] * 100
    ax_b.text(count + 40, i, f'{count:,} ({pct:.1f}%)',
              va='center', ha='left', fontsize=5.5,
              color=C['black'])

ax_b.set_yticks(y_pos)
ax_b.set_yticklabels(src_names, fontsize=6.5)
ax_b.set_xlabel('Compounds')
ax_b.set_xlim(0, max(src_counts) * 1.70)   # wider to avoid label clip
# y-axis already labels each source — no redundant legend needed
panel_label(ax_b, 'b', x=-0.22, y=1.04)

# ── Panel c: Independence score stacked bar ────────────────────────────────────
print("[Panel c] Independence score stacked bar ...")
pct_nonred = n_nonred / CANON['n_compounds'] * 100   # 17.1%
pct_red    = n_red    / CANON['n_compounds'] * 100   # 82.9%

ax_c.barh([0], [n_nonred], color=C['blue'],  height=0.55, label=f'Non-redundant (score ≥ 0.6)')
ax_c.barh([0], [n_red],    left=[n_nonred],  color='#AAAAAA', height=0.55, label='Pipeline-redundant')

# Annotate percentages inside bars
ax_c.text(n_nonred / 2, 0, f'{pct_nonred:.1f}%\nn={n_nonred:,}',
          ha='center', va='center', fontsize=5.5, color='white', fontweight='bold')
ax_c.text(n_nonred + n_red / 2, 0, f'{pct_red:.1f}%\nn={n_red:,}',
          ha='center', va='center', fontsize=5.5, color='white', fontweight='bold')

ax_c.set_xlim(0, CANON['n_compounds'])
ax_c.set_yticks([])
ax_c.set_xlabel('Compounds')
ax_c.set_title('Source independence', fontsize=6, pad=3)
ax_c.legend(fontsize=5.5, loc='upper right',
            framealpha=0.85, edgecolor='none',
            bbox_to_anchor=(1.0, -0.12))
panel_label(ax_c, 'c', x=-0.22, y=1.04)

# ── Save ──────────────────────────────────────────────────────────────────────
save_fig(fig, str(OUT / 'fig02_source_overlap'))
plt.close(fig)
print("Figure 2 complete.")
