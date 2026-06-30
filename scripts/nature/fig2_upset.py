#!/usr/bin/env python
"""Nature Fig 2 — Source Overlap UpSet Plot."""
import os, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
from upsetplot import UpSet, from_memberships

NATURE_RC = {
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial','Helvetica','DejaVu Sans'],
    'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7, 'axes.linewidth': 0.75,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': False,
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'lines.linewidth': 0.75, 'patch.linewidth': 0.5,
    'legend.fontsize': 6, 'legend.frameon': False,
    'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    'figure.constrained_layout.use': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
}
matplotlib.rcParams.update(NATURE_RC)

PAL = {
    'blue': '#0077BB', 'orange': '#EE7733', 'red': '#CC3311',
    'teal': '#009988', 'navy': '#1A237E',
    'gray_light': '#BBBBBB', 'gray_dark': '#555555',
}
OUT = 'outputs/figures/nature'
os.makedirs(OUT, exist_ok=True)

def save_fig(fig, name):
    for ext in ('png', 'svg', 'pdf'):
        p = f'{OUT}/{name}.{ext}'
        fig.savefig(p, dpi=600 if ext == 'png' else None,
                    bbox_inches='tight', facecolor='white')
    sz = os.path.getsize(f'{OUT}/{name}.png') / 1024
    print(f"Saved {name}: {sz:.0f} KB")

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_parquet('data/processed/pad4_compounds.parquet')

# Build membership list (pipe-delimited → list of sets)
source_map = {
    'pubchem_confirmatory': 'PubChem',
    'chembl':               'ChEMBL',
    'bindingdb':            'BindingDB',
}
memberships = []
for sl in df['source_list']:
    parts = [source_map.get(s.strip(), s.strip()) for s in sl.split('|')]
    memberships.append(parts)

data = from_memberships(memberships, data=pd.Series(np.ones(len(df), dtype=int)))

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(7.2, 4))
upset = UpSet(data, subset_size='count', show_counts=True,
              totals_plot_elements=3, sort_by='cardinality',
              element_size=30)
d = upset.plot(fig)

# Color bars: PubChem-only (count=233) → orange; all others → blue
ax_inter = d['intersections']
for patch in ax_inter.patches:
    h = patch.get_height()
    if abs(h - 233) < 1:
        patch.set_facecolor(PAL['orange'])
    else:
        patch.set_facecolor(PAL['blue'])
    patch.set_edgecolor('none')

# Remove internal grid
ax_inter.yaxis.grid(False)
ax_inter.set_ylim(0, 1500)
ax_inter.set_ylabel('Compound count', fontsize=7)

# Label the totals (left) bar chart axis
if 'totals' in d:
    d['totals'].set_xlabel('Compounds\nper source', fontsize=6)

# Style the totals bars
if 'totals' in d:
    ax_tot = d['totals']
    for patch in ax_tot.patches:
        patch.set_facecolor(PAL['gray_light'])
        patch.set_edgecolor('none')

# Legend and footnote
handles = [
    mpatches.Patch(facecolor=PAL['blue'],   edgecolor='none', label='Multi-source combinations'),
    mpatches.Patch(facecolor=PAL['orange'], edgecolor='none', label='PubChem-only (n=233)'),
]
ax_inter.legend(handles=handles, loc='upper right', fontsize=6, frameon=False,
                handlelength=1.2, borderpad=0.3)

fig.text(0.02, 0.01,
         "N=3,093 curated PAD4 inhibitors. Intersection sizes above bars. "
         "Source labels: PubChem=95 confirmatory AIDs; ChEMBL=CHEMBL6111; BindingDB=Q9UM07.",
         fontsize=5.5, color=PAL['gray_dark'], ha='left', va='bottom',
         fontfamily='sans-serif')

save_fig(fig, 'fig2_upset')
plt.close()
print("Fig 2 DONE")
