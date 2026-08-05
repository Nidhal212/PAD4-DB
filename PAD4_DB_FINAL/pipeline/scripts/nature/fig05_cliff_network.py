"""
fig05_cliff_network.py — Figure 5: Cliff Network (2-panel)
Outputs: outputs/figures/nature/fig05_cliff_network.{png,pdf}
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from pathlib import Path

try:
    import networkx as nx
    HAS_NX = True
except ImportError:
    HAS_NX = False
    print("WARNING: networkx not available")

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
    'navy':       '#004488',
    'red':        '#CC3311',
    'grey':       '#BBBBBB',
}

HUB_IKS = {
    'A1': 'SMADULGDNOCLOP-GISFHXKWSA-N',
    'A2': 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
    'B1': 'UDCDEKJNAMHBFH-HSZRJFAPSA-N',
    'B2': 'DVCKJOQIVOGXEI-XMMPIXPASA-N',
}

# Manual per-hub label offsets (data units). B1/B2 sit close together in the
# layout (Tanimoto 0.975 near-duplicates), so they get pulled apart hardest.
# Tune these if a future re-run of spring_layout shifts node positions.
HUB_LABEL_OFFSET = {
    'A1': (0.35, 0.35),
    'A2': (0.35, -0.35),
    'B1': (-0.55, 0.45),
    'B2': (0.20, -0.55),
}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/figures/nature_v2'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FIGURE 5 — CLIFF NETWORK (Fixed Legend/Label Layout)")
print("=" * 60)

def get_col(df, candidates):
    for col in candidates:
        if col in df.columns:
            return col
    raise KeyError(f"None of {candidates} found in DataFrame columns: {df.columns.tolist()}")

df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
cliffs = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
disc = pd.read_csv(ROOT / 'outputs/mmp/mmp_discontinuity_scores.csv')

mech_col = get_col(df, ['mechanism_class', 'assay_mechanism_classes'])
tier_col = get_col(cliffs, ['cliff_tier', 'Cliff_Tier'])
severe = cliffs[cliffs[tier_col].str.lower() == 'severe'].copy()
print(f"  Severe pairs: {len(severe)}")

if HAS_NX:
    G = nx.Graph()
    for _, row in severe.iterrows():
        G.add_edge(row['inchi_key_a'], row['inchi_key_b'],
                   delta=abs(row['delta_pic50']))

    pic50_map = dict(zip(df['inchi_key'], df['pIC50']))
    nx.set_node_attributes(G, pic50_map, 'pIC50')

    # Sanity check: warn if any node's pIC50 falls outside the fixed
    # colour-scale range (2-8) used below, so points don't get silently
    # clipped to the scale's endpoint colour.
    pic50_vals = [G.nodes[n]['pIC50'] for n in G.nodes()]
    if min(pic50_vals) < 2 or max(pic50_vals) > 8:
        print(f"  WARNING: pIC50 range [{min(pic50_vals):.2f}, {max(pic50_vals):.2f}] "
              f"exceeds fixed colour scale [2, 8] — those nodes will clip.")

    fig, axes = plt.subplots(1, 2, figsize=(11, 5),
                              gridspec_kw={'width_ratios': [1.6, 1]},
                              constrained_layout=True)
    ax_a, ax_b = axes

    # ── Panel a: Force-directed network ────────────────────────────────────────
    print("\n[Panel a] Drawing force-directed network...")
    pos = nx.spring_layout(G, seed=17, k=1.5, scale=4.5, iterations=200)

    # Edge color mapping
    edge_deltas = [d['delta'] for u, v, d in G.edges(data=True)]
    min_delta, max_delta = min(edge_deltas), max(edge_deltas)
    cmap_edges = mcolors.LinearSegmentedColormap.from_list('edge_cmap', ['#BBBBBB', '#CC3311'])

    # Draw Edges
    for u, v, data in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        norm_delta = (data['delta'] - min_delta) / (max_delta - min_delta + 1e-9)
        color = cmap_edges(norm_delta)
        ax_a.plot([x0, x1], [y0, y1], color=color, lw=0.8, alpha=0.6, zorder=1)

    # Draw Non-hub Nodes
    hub_set = set(HUB_IKS.values())
    for node in G.nodes():
        if node not in hub_set:
            x, y = pos[node]
            pic50 = G.nodes[node]['pIC50']
            size = G.degree(node) * 20 + 10
            ax_a.scatter(x, y, s=size, c=pic50, cmap='viridis', vmin=2, vmax=8,
                          edgecolors='white', lw=0.3, alpha=0.9, zorder=2)

    # Draw Hub Nodes
    hub_a_nodes = [HUB_IKS['A1'], HUB_IKS['A2']]
    hub_b_nodes = [HUB_IKS['B1'], HUB_IKS['B2']]

    for node in hub_a_nodes:
        if node in G:
            x, y = pos[node]
            pic50 = G.nodes[node]['pIC50']
            ax_a.scatter(x, y, s=450, marker='*', c=pic50, cmap='viridis', vmin=2, vmax=8,
                          edgecolors='white', lw=1.5, zorder=3)

    for node in hub_b_nodes:
        if node in G:
            x, y = pos[node]
            pic50 = G.nodes[node]['pIC50']
            ax_a.scatter(x, y, s=400, marker='D', c=pic50, cmap='viridis', vmin=2, vmax=8,
                          edgecolors='white', lw=1.5, zorder=3)

    # Hub Labels — annotate with a leader line so labels stay legible even
    # when two hubs sit close together in the layout (e.g. B1/B2).
    for lbl, ik in HUB_IKS.items():
        if ik in pos:
            x, y = pos[ik]
            dx, dy = HUB_LABEL_OFFSET.get(lbl, (0.3, 0.3))
            ax_a.annotate(
                lbl, xy=(x, y), xytext=(x + dx, y + dy),
                fontsize=10, fontweight='bold', ha='center', va='center', color='#333333',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#999999', lw=0.5, alpha=0.95),
                arrowprops=dict(arrowstyle='-', color='#777777', lw=0.6, shrinkA=1, shrinkB=8),
                zorder=6,
            )

    # Colorbar (placed between panel A and B)
    norm = mcolors.Normalize(vmin=2, vmax=8)
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    sm._A = []
    cbar = fig.colorbar(sm, ax=ax_a, shrink=0.6, pad=0.02, aspect=20)
    cbar.set_label('pIC50', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax_a.set_title(f'{G.number_of_nodes()} nodes · {G.number_of_edges()} severe cliff edges', fontsize=9, pad=10)

    # ── Legend: reserved dead-space band below all nodes ─────────────────────
    # Spring layout is organic — there is no corner that's reliably empty of
    # nodes across reruns/data updates. Instead of guessing a "clear" corner,
    # carve out an explicit band strictly below every node's y-coordinate
    # (edges are straight lines between two nodes, so they never dip below
    # the lower node's y either — the band is guaranteed empty by construction).
    x_vals = [p[0] for p in pos.values()]
    y_vals = [p[1] for p in pos.values()]
    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)
    padding_x = (x_max - x_min) * 0.15
    padding_y = (y_max - y_min) * 0.15
    legend_band = (y_max - y_min) * 0.34  # extra room reserved for legend

    ax_a.set_xlim(x_min - padding_x, x_max + padding_x)
    ax_a.set_ylim(y_min - padding_y - legend_band, y_max + padding_y)

    # Single combined legend (node encoding + edge encoding as one block).
    # A two-box stacked legend (title/handles, title/handles) is fragile:
    # the first box's true rendered height depends on font metrics and isn't
    # known until draw time, so a fixed axes-fraction gap between the two
    # anchors will eventually clip one into the other (it did, on the first
    # pass here — the "Non-hub" row was cut by the edge-legend box on top of
    # it). One legend with descriptive labels avoids the stacking problem
    # entirely. Edge ΔpIC50 endpoints are read from the actual plotted
    # min/max rather than hardcoded, so the legend can't drift out of sync
    # with the data on a re-run.
    legend_handles = [
        Line2D([0], [0], marker='*', color='w', label='Hub A — series floor (n=2)',
               markeredgecolor=COLORS['navy'], markerfacecolor=COLORS['navy'], markersize=12),
        Line2D([0], [0], marker='D', color='w', label='Hub B — singleton attractor (n=2)',
               markeredgecolor=COLORS['red'], markerfacecolor=COLORS['red'], markersize=10),
        Line2D([0], [0], marker='o', color='w', label='Non-hub (size scales w/ degree, colour = pIC50)',
               markeredgecolor='#888888', markerfacecolor='#888888', markersize=8),
        mpatches.Patch(color='#BBBBBB', label=f'Edge ΔpIC50, low (~{min_delta:.2f})'),
        mpatches.Patch(color='#CC3311', label=f'Edge ΔpIC50, high (~{max_delta:.2f})'),
    ]
    leg = ax_a.legend(handles=legend_handles,
                       loc='upper left', bbox_to_anchor=(0.0, 0.0), bbox_transform=ax_a.transAxes,
                       fontsize=7.5, framealpha=0.95, edgecolor='#CCCCCC', borderaxespad=0.3)
    leg.set_zorder(10)
    # ──────────────────────────────────────────────────────────────────────────

    ax_a.set_xticks([])
    ax_a.set_yticks([])
    for spine in ax_a.spines.values():
        spine.set_visible(False)
    ax_a.text(-0.05, 1.05, 'a', transform=ax_a.transAxes, fontsize=11, fontweight='bold')

    # ── Panel b: Degree bar chart ──────────────────────────────────────────────
    print("\n[Panel b] Degree bar chart ...")
    disc_top = disc.sort_values('severe_cliff_degree', ascending=False).head(12)

    y_labels_b = []
    bar_colors_b = []
    nonhub_counter = 0
    id_map = []  # (short_id, full_inchi_key) for non-hub compounds — write out alongside the figure

    for ik in disc_top['inchi_key']:
        if ik == HUB_IKS['A1']:
            y_labels_b.append('Hub A1')
            bar_colors_b.append(COLORS['navy'])
        elif ik == HUB_IKS['A2']:
            y_labels_b.append('Hub A2')
            bar_colors_b.append(COLORS['navy'])
        elif ik == HUB_IKS['B1']:
            y_labels_b.append('Hub B1')
            bar_colors_b.append(COLORS['red'])
        elif ik == HUB_IKS['B2']:
            y_labels_b.append('Hub B2')
            bar_colors_b.append(COLORS['red'])
        else:
            nonhub_counter += 1
            short_id = f'C{nonhub_counter}'
            y_labels_b.append(short_id)
            bar_colors_b.append(COLORS['grey'])
            id_map.append((short_id, ik))

    y_pos = range(len(disc_top))
    ax_b.barh(list(y_pos), disc_top['severe_cliff_degree'].values,
               color=bar_colors_b, height=0.6, edgecolor='white', lw=0.3)

    for i, (val, ik) in enumerate(zip(disc_top['severe_cliff_degree'], disc_top['inchi_key'])):
        ax_b.text(val + 0.1, i, str(int(val)), va='center', ha='left', fontsize=7.5)

    ax_b.set_yticks(list(y_pos))
    ax_b.set_yticklabels(y_labels_b, fontsize=8, family='sans-serif')
    ax_b.set_xlabel('Severe cliff partners (network degree)', fontsize=9)
    ax_b.invert_yaxis()
    ax_b.spines['top'].set_visible(False)
    ax_b.spines['right'].set_visible(False)
    ax_b.text(-0.15, 1.02, 'b', transform=ax_b.transAxes, fontsize=11, fontweight='bold')

    # ── Save ────────────────────────────────────────────────────────────────────
    for ext in ['png', 'pdf']:
        outpath = OUT / f'fig05_cliff_network.{ext}'
        fig.savefig(outpath, dpi=300, bbox_inches='tight')
        print(f"  Saved: {outpath}")

    plt.close(fig)

    # Write the C1..C8 -> full InChIKey mapping so panel b labels are
    # traceable. Previously these were truncated 14-char keys with no
    # lookup anywhere in the manuscript (Table 4 only covers the 4 hubs).
    if id_map:
        map_path = OUT / 'fig05_panel_b_nonhub_id_map.csv'
        pd.DataFrame(id_map, columns=['short_id', 'inchi_key']).to_csv(map_path, index=False)
        print(f"  Saved: {map_path}  (panel b non-hub short-ID lookup — cite as Supplementary Table)")

    print(f"\nFigure 5 complete.")
else:
    print("ERROR: networkx is required to generate this figure.")