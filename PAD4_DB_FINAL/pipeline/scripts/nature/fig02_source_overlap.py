#!/usr/bin/env python
"""Nature Fig 2 — 3-Panel Source Overlap (Matplotlib-only)"""
import os, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
from pathlib import Path

NATURE_RC = {
    'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial'],
    'axes.linewidth': 0.75, 'xtick.major.width': 0.75, 'ytick.major.width': 0.75,
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'axes.labelsize': 7, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'legend.fontsize': 6, 'legend.frameon': False,
    'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    'pdf.fonttype': 42, 'ps.fonttype': 42
}
matplotlib.rcParams.update(NATURE_RC)

PAL = {
    'blue': '#0077BB', 'orange': '#EE7733', 'red': '#CC3311',
    'teal': '#009988', 'cyan': '#33BBEE', 'navy': '#1A237E',
    'grey': '#BBBBBB', 'dark_grey': '#555555', 'light_grey': '#E8E8E8'
}

OUT = Path('/home/nidhal/PAD4-db_V2/outputs/figures/nature')
OUT.mkdir(parents=True, exist_ok=True)

def save_fig(fig, name):
    for ext in ('png', 'svg', 'pdf'):
        p = OUT / f'{name}.{ext}'
        fig.savefig(p, dpi=600 if ext == 'png' else None, bbox_inches='tight')
    print(f"✅ Saved: {OUT / name}.png")

def plabel(ax, letter):
    ax.text(-0.15, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight='bold')

print("=" * 60)
print("FIGURE 2 — SOURCE OVERLAP (3-PANEL)")
print("=" * 60)

# ── Load Data ─────────────────────────────────────────────────────────────────
df = pd.read_parquet('data/processed/pad4_compounds.parquet')
print(f"Loaded {len(df)} compounds.")

# ── Helper: Source counts ────────────────────────────────────────────────────
def has_source(s, name):
    return s.str.contains(name, na=False)

df['pc'] = has_source(df['source_list'], 'pubchem')
df['cb'] = has_source(df['source_list'], 'chembl')
df['bd'] = has_source(df['source_list'], 'bindingdb')

# Derive combos for UpSet
def get_combo(row):
    c = []
    if row['pc']: c.append('PubChem')
    if row['cb']: c.append('ChEMBL')
    if row['bd']: c.append('BindingDB')
    return '+'.join(c)
df['combo'] = df.apply(get_combo, axis=1)

combo_counts = df['combo'].value_counts()
ordered_combos = combo_counts.index.tolist()
ordered_counts = combo_counts.values.tolist()
n_combos = len(ordered_combos)

# ── Build Figure (3 Panels) ──────────────────────────────────────────────────
fig = plt.figure(figsize=(10, 4.5))
gs = fig.add_gridspec(1, 3, width_ratios=[2.5, 1, 1], wspace=0.4)

# ════════════════════════════════════════════════════
# PANEL A: UpSet Plot (Manual)
# ════════════════════════════════════════════════════
ax_a = fig.add_subplot(gs[0])
ax_a.remove() # Remove blank frame
gs_inner = gs[0].subgridspec(2, 1, height_ratios=[2.5, 1], hspace=0.05)
ax_bar = fig.add_subplot(gs_inner[0])
ax_mat = fig.add_subplot(gs_inner[1])

sources = ['PubChem', 'ChEMBL', 'BindingDB']

# Top bar chart
bar_colors = []
for combo in ordered_combos:
    if combo == 'PubChem': 
        bar_colors.append(PAL['orange'])
    else: 
        bar_colors.append(PAL['blue'])

bars = ax_bar.bar(range(n_combos), ordered_counts, color=bar_colors, width=0.7, edgecolor='white', lw=0.5)
max_count = max(ordered_counts) if ordered_counts else 1
for i, (count, bar) in enumerate(zip(ordered_counts, bars)):
    ax_bar.text(i, count + max_count * 0.1, str(count),
                ha='center', va='bottom', fontsize=7, fontweight='bold')

ax_bar.set_ylabel('Compounds', fontsize=7)
ax_bar.set_xlim(-0.6, n_combos - 0.4)
ax_bar.set_xticks([])
ax_bar.spines['top'].set_visible(False)
ax_bar.spines['right'].set_visible(False)

# ── Legend for Panel A ──────────────────────────────
legend_handles = [
    mpatches.Patch(facecolor=PAL['blue'],   edgecolor='none', label='Multi-source combinations'),
    mpatches.Patch(facecolor=PAL['orange'], edgecolor='none', label='PubChem-only (n=233)'),
]
ax_bar.legend(handles=legend_handles, loc='upper right', fontsize=6, frameon=False,
              handlelength=1.2, borderpad=0.3)

# Bottom matrix (dots + connecting lines)
ax_mat.set_xlim(-0.6, n_combos - 0.4)
ax_mat.set_ylim(-0.5, len(sources) - 0.5)
ax_mat.set_yticks(range(len(sources)))
ax_mat.set_yticklabels(sources[::-1], fontsize=6)
ax_mat.set_xticks([])
ax_mat.spines[:].set_visible(False)

for i, combo in enumerate(ordered_combos):
    members = combo.split('+')
    for j, src in enumerate(sources[::-1]):
        filled = src in members
        ax_mat.scatter([i], [j],
                       s=80 if filled else 40,
                       c=PAL['dark_grey'] if filled else PAL['light_grey'],
                       zorder=3 if filled else 2, edgecolors='none')
    # Draw vertical lines connecting filled dots
    filled_ys = [j for j, src in enumerate(sources[::-1]) if src in members]
    if len(filled_ys) > 1:
        ax_mat.plot([i, i], [min(filled_ys), max(filled_ys)], 
                    color=PAL['dark_grey'], lw=1.5, zorder=2)

plabel(ax_bar, 'a')

# ════════════════════════════════════════════════════
# PANEL B: Per-source Coverage
# ════════════════════════════════════════════════════
ax_b = fig.add_subplot(gs[1])
src_names = ['PubChem', 'BindingDB', 'ChEMBL']
src_counts = [df['pc'].sum(), df['bd'].sum(), df['cb'].sum()]
src_colors = [PAL['blue'], PAL['teal'], PAL['navy']]

bars_b = ax_b.barh(range(len(src_names)), src_counts, color=src_colors, height=0.6, edgecolor='white', lw=0.5)
for i, (count, bar) in enumerate(zip(src_counts, bars_b)):
    ax_b.text(count + 20, i, f'{count:,}', va='center', ha='left', fontsize=7, fontweight='bold')

ax_b.set_yticks(range(len(src_names)))
ax_b.set_yticklabels(src_names, fontsize=7)
ax_b.set_xlim(0, max(src_counts) * 1.25)
ax_b.set_xlabel('Compounds', fontsize=7)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
plabel(ax_b, 'b')

# ════════════════════════════════════════════════════
# PANEL C: Source Independence
# ════════════════════════════════════════════════════
ax_c = fig.add_subplot(gs[2])
non_red = 528
pipe_red = 2565

ax_c.barh(0, non_red, height=0.6, color=PAL['blue'], edgecolor='white', lw=0.5)
ax_c.barh(1, pipe_red, height=0.6, color=PAL['grey'], edgecolor='white', lw=0.5)

ax_c.text(non_red + 10, 0, f'{non_red:,}', va='center', ha='left', fontsize=7, fontweight='bold')
ax_c.text(pipe_red + 10, 1, f'{pipe_red:,}', va='center', ha='left', fontsize=7, fontweight='bold')

ax_c.set_yticks([0, 1])
ax_c.set_yticklabels(['Non-redundant', 'Redundant'], fontsize=7)
ax_c.set_xlim(0, pipe_red * 1.15)
ax_c.set_xlabel('Compounds', fontsize=7)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
plabel(ax_c, 'c')

# ── Save Final Figure ─────────────────────────────────────────────────────────
save_fig(fig, 'fig02_source_overlap')
plt.close()
print("\n🎉 Figure 2 (3-panel) generated successfully!")