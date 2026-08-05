#!/usr/bin/env python
"""Nature Fig 2 — Source Overlap (UpSet + Coverage + Independence)"""
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

# ── RC Params ─────────────────────────────────────────────────────────────────
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

def plabel(ax, letter, x=-0.15, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='bottom', ha='right')

# ── Load Data ─────────────────────────────────────────────────────────────────
df = pd.read_parquet('data/processed/pad4_compounds.parquet')

# ── Prepare data for UpSet (Panel A) ─────────────────────────────────────────
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

# ── Prepare data for Panel B (Coverage) ──────────────────────────────────────
def has_source(s, name):
    return s.str.contains(name, na=False)
n_pc = has_source(df['source_list'], 'pubchem').sum()
n_cb = has_source(df['source_list'], 'chembl').sum()
n_bd = has_source(df['source_list'], 'bindingdb').sum()

# ── Prepare data for Panel C (Source Independence) ──────────────────────────
# Hardcoded based on the paper's exact numbers
n_non_redundant = 528
n_pipeline_redundant = 2565

# ── Build Figure with GridSpec (3 columns: Wide, Narrow, Narrow) ─────────────
fig = plt.figure(figsize=(10, 4.5))
gs = fig.add_gridspec(1, 3, width_ratios=[3.5, 1.1, 1.1], wspace=0.4)

# ── Panel A: UpSet Plot ──────────────────────────────────────────────────────
ax_a = fig.add_subplot(gs[0])
upset = UpSet(data, subset_size='count', show_counts=True,
              totals_plot_elements=3, sort_by='cardinality',
              element_size=25)
d = upset.plot(fig)

# Transfer the UpSet axes to the specific subplot slot (tricky with upsetplot).
# To avoid breaking the layout, let's manually reposition the UpSet axes.
# The `upset.plot(fig)` will create its own axes. We need to rely on `d['intersections']` etc.
# Instead of fighting upsetplot layout, we can just draw Panel A on the left of the figure.
# Re-creating the UpSet plot directly onto the gs[0] can be done, but to keep it bulletproof:
# We generate the UpSet in a standalone figure, then extract it. 
# But for simplicity, let's just let upsetplot handle Panel A, and we create Panels B & C separately.

# **CRITICAL CHANGE TO MAKE THIS WORK WITHOUT WRECKING LAYOUT:**
# We will manually craft Panel A using `upsetplot` but *without* the totals plot to keep it clean.
# Wait, standard `upsetplot` doesn't easily let you put into a custom subplot.
# Let me suggest an easier fix: Keep the code simple. Combine Panels B and C into the manual script.
# Actually, your second script uses `plt.figure()`. Let's just use `subplots`.
# Let's write a clean, non-breaking 3-panel approach.

# ----- NEW APPROACH FOR PANEL A (Manual UpSet using the library's `plot` method) -----
# We create the UpSet plot, get the figure, and crop it. 
# But the BEST fix for you: Use standard Matplotlib for Panels B and C, and use the Upset library for Panel A.
# Let's generate the Upset in its own figure, then save it separately if needed.
# Actually, writing a custom UpSet in matplotlib is risky. Let's just provide the code that adds Panels B & C.

# Wait, I see what the user wants. I'll rewrite the `fig2_upset.py` script to include Panel B and C.
# Since `upsetplot` hijacks the figure, we can put Panel B and C in a *separate* figure, or use subplots.
# Let's use a simpler, robust approach: generate the Upset figure, and manually add the side plots using `gridspec`.

fig = plt.figure(figsize=(10, 4.5))
gs = fig.add_gridspec(1, 3, width_ratios=[3, 1, 1.2], wspace=0.3)

# Panel A (UpSet)
ax_a = fig.add_subplot(gs[0])
upset = UpSet(data, subset_size='count', show_counts=True,
              totals_plot_elements=3, sort_by='cardinality',
              element_size=25)
d = upset.plot(fig)

# The Upset plot will overlay on the whole figure by default. 
# We need to force it onto the left subplot by creating a new figure for it, OR
# Actually, the standard `upset.plot` doesn't support `ax=` directly.
# **The simplest, bulletproof solution for you is to keep Panel A as its own plot, and Panels B & C as a separate figure.**
# OR: Use the manual matplotlib script for all 3 panels (Panel A manually drawn) to guarantee 1 figure.

# Since I must give you a working solution, let's go with the Matplotlib-only version (Script 1) but fix it completely and add Panel C.

# Let's fix `fig02_source_overlap.py` instead.
# I will provide a completely rewritten, robust version of `fig02_source_overlap.py` that includes all 3 panels, fixes the text overlap, and uses clean matplotlib.

# Let's craft a custom combined Matplotlib-only script for Figure 2 (3 panels).