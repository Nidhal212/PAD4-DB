#!/usr/bin/env python
"""Nature Fig 7 — Cliff Network + Degree Bar (manual axes layout)."""
import os, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

NATURE_RC = {
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial','Helvetica','DejaVu Sans'],
    'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7, 'axes.linewidth': 0.75,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': False,
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'lines.linewidth': 0.75, 'lines.markersize': 4, 'patch.linewidth': 0.5,
    'legend.fontsize': 6, 'legend.frameon': False,
    'legend.handlelength': 1.5, 'legend.handletextpad': 0.5,
    'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    'figure.constrained_layout.use': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
}
matplotlib.rcParams.update(NATURE_RC)

PAL = {
    'blue': '#0077BB', 'orange': '#EE7733', 'red': '#CC3311',
    'teal': '#009988', 'cyan': '#33BBEE', 'navy': '#1A237E',
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

def plabel(ax, letter, x=-0.08, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=8, fontweight='bold', va='bottom', ha='right',
            fontfamily='sans-serif')

# ── Load network data ─────────────────────────────────────────────────────────
try:
    import networkx as nx
    has_nx = True
except ImportError:
    has_nx = False

nodes_df = pd.read_csv('outputs/figures/fig7_nodes.csv')
edges_df = pd.read_csv('outputs/figures/fig7_edges.csv')

HUB_A = {'SMADULGDNOCLOP-GISFHXKWSA-N', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N'}
HUB_B = {'UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'DVCKJOQIVOGXEI-XMMPIXPASA-N'}
hub_labels = {
    'SMADULGDNOCLOP-GISFHXKWSA-N': 'A1',
    'RAVBZQAQTVGKIV-XBPDSQQVSA-N': 'A2',
    'UDCDEKJNAMHBFH-HSZRJFAPSA-N': 'B1',
    'DVCKJOQIVOGXEI-XMMPIXPASA-N': 'B2',
}
mech_colors = {
    'enzymatic':           PAL['gray_light'],
    'enzymatic_confirmed': PAL['teal'],
    'fp_ic50':             PAL['cyan'],
    'covalent':            PAL['red'],
}

# ── Build graph and layout ─────────────────────────────────────────────────────
if has_nx:
    G = nx.Graph()
    for _, row in nodes_df.iterrows():
        G.add_node(row['inchi_key'], **row.to_dict())
    for _, row in edges_df.iterrows():
        G.add_edge(row['inchi_key_a'], row['inchi_key_b'],
                   same_mechanism=row.get('same_mechanism', True),
                   hub_class_involved=row.get('hub_class_involved', 'none'))

    np.random.seed(42)
    try:
        pos = nx.spring_layout(G, k=1.8, iterations=100, seed=42)
    except Exception:
        pos = nx.spectral_layout(G)
    node_order = list(G.nodes())
else:
    # Fallback: random positions
    rng = np.random.default_rng(42)
    iks = nodes_df['inchi_key'].tolist()
    pos = {ik: rng.uniform(-1, 1, 2) for ik in iks}
    node_order = iks

# ── Figure: manual axes (network left 60%, bar right 38%) ─────────────────────
fig = plt.figure(figsize=(7.2, 3.8))
ax_net = fig.add_axes([0.01, 0.08, 0.57, 0.88])
ax_bar = fig.add_axes([0.63, 0.08, 0.35, 0.88])

ax_net.set_axis_off()

# Draw edges
for _, row in edges_df.iterrows():
    ia, ib = row['inchi_key_a'], row['inchi_key_b']
    if ia not in pos or ib not in pos:
        continue
    xa, ya = pos[ia]
    xb, yb = pos[ib]
    is_cross = not bool(row.get('same_mechanism', True))
    hc = row.get('hub_class_involved', 'none')
    if is_cross:
        clr, lw, ls = PAL['gray_dark'], 0.5, (0, (2, 2))  # tight dash
    elif hc == 'A':
        clr, lw, ls = PAL['navy'], 0.65, 'solid'
    elif hc == 'B':
        clr, lw, ls = PAL['red'], 0.65, 'solid'
    else:
        clr, lw, ls = PAL['gray_light'], 0.4, 'solid'
    ax_net.plot([xa, xb], [ya, yb], color=clr, linewidth=lw, linestyle=ls,
                alpha=0.55, zorder=1, solid_capstyle='round')

# Draw non-hub nodes first
for ik in node_order:
    if ik in HUB_A or ik in HUB_B:
        continue
    xi, yi = pos[ik]
    row = nodes_df[nodes_df['inchi_key'] == ik].iloc[0]
    mech = row.get('mechanism_class', 'enzymatic')
    clr = mech_colors.get(mech, PAL['gray_light'])
    ax_net.scatter(xi, yi, s=7, c=clr, linewidths=0.3, edgecolors='white',
                   alpha=0.85, zorder=2)

# Draw hub nodes on top with labels
for ik, lbl in hub_labels.items():
    if ik not in pos:
        continue
    xi, yi = pos[ik]
    clr = PAL['navy'] if ik in HUB_A else PAL['red']
    ax_net.scatter(xi, yi, s=55, c=clr, linewidths=0.5, edgecolors='white',
                   alpha=1.0, zorder=4)
    ax_net.text(xi, yi, lbl, ha='center', va='center', fontsize=5.5,
                fontweight='bold', color='white', fontfamily='sans-serif', zorder=5)

# Legend
leg_elements = [
    mpatches.Patch(facecolor=PAL['gray_light'], label='Enzymatic (non-confirmed)'),
    mpatches.Patch(facecolor=PAL['teal'],       label='Enzymatic confirmed'),
    mpatches.Patch(facecolor=PAL['cyan'],        label='FP-based IC50'),
    mpatches.Patch(facecolor=PAL['navy'],        label='Hub Class A'),
    mpatches.Patch(facecolor=PAL['red'],         label='Hub Class B'),
    mlines.Line2D([], [], color=PAL['navy'],      lw=0.65, label='Hub A edge'),
    mlines.Line2D([], [], color=PAL['red'],       lw=0.65, label='Hub B edge'),
    mlines.Line2D([], [], color=PAL['gray_light'],lw=0.4,  label='Non-hub edge'),
    mlines.Line2D([], [], color=PAL['gray_dark'], lw=0.5,
                  linestyle=(0,(2,2)), label='Cross-mechanism (n=4)'),
]
# Trim to 6 legend entries
leg_elements_6 = [
    mpatches.Patch(facecolor=PAL['navy'],        label='Hub Class A (n=2, 27 pairs)'),
    mpatches.Patch(facecolor=PAL['red'],         label='Hub Class B (n=2, 23 pairs)'),
    mpatches.Patch(facecolor=PAL['teal'],        label='Enzymatic confirmed'),
    mpatches.Patch(facecolor=PAL['gray_light'],  label='Enzymatic / FP-IC50'),
    mlines.Line2D([], [], color=PAL['gray_light'], lw=0.4, label='Non-hub edge'),
    mlines.Line2D([], [], color=PAL['gray_dark'],  lw=0.5,
                  linestyle=(0,(2,2)), label='Cross-mechanism (n=4)'),
]
ax_net.legend(handles=leg_elements_6, loc='lower left', fontsize=5.0,
              frameon=False, handlelength=0.9, labelspacing=0.18,
              borderpad=0.2, ncol=1, handletextpad=0.4)
# No panel title — removed per reviewer

# Panel label A
ax_net.text(-0.03, 1.01, 'A', transform=ax_net.transAxes,
            fontsize=8, fontweight='bold', va='bottom', ha='left')

# ── Panel B: Degree distribution bar ─────────────────────────────────────────
top_n_deg = 20
from collections import Counter
all_iks = edges_df['inchi_key_a'].tolist() + edges_df['inchi_key_b'].tolist()
deg_count = Counter(all_iks)

deg_df_sorted = pd.DataFrame({'ik': list(deg_count.keys()),
                               'deg': list(deg_count.values())})
deg_df_sorted = deg_df_sorted.sort_values('deg', ascending=False).head(top_n_deg)

bar_colors = []
for ik in deg_df_sorted['ik']:
    if ik in HUB_A:
        bar_colors.append(PAL['navy'])
    elif ik in HUB_B:
        bar_colors.append(PAL['red'])
    else:
        bar_colors.append(PAL['gray_light'])

bar_labels = []
for ik in deg_df_sorted['ik']:
    if ik in hub_labels:
        bar_labels.append(hub_labels[ik])
    else:
        bar_labels.append(ik[:8] + '…')

y_pos = range(len(deg_df_sorted))
ax_bar.barh(list(y_pos), deg_df_sorted['deg'].values,
            color=bar_colors, linewidth=0, height=0.7)
ax_bar.set_yticks(list(y_pos))
ax_bar.set_yticklabels(bar_labels, fontsize=5.2)
ax_bar.invert_yaxis()
ax_bar.set_xlabel('Number of severe cliff partners', fontsize=7)
# No panel title — reviewer request
ax_bar.spines['left'].set_visible(True)
ax_bar.spines['bottom'].set_visible(True)

# Color y-tick labels for hubs
for ytick, ik in zip(ax_bar.get_yticklabels(), deg_df_sorted['ik']):
    if ik in HUB_A:
        ytick.set_color(PAL['navy'])
        ytick.set_fontweight('bold')
    elif ik in HUB_B:
        ytick.set_color(PAL['red'])
        ytick.set_fontweight('bold')

ax_bar.text(-0.18, 1.01, 'B', transform=ax_bar.transAxes,
            fontsize=8, fontweight='bold', va='bottom', ha='left')

save_fig(fig, 'fig7_network')
plt.close()
print("Fig 7 DONE")
