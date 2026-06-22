"""
fig05_cliff_network.py — Figure 5: Cliff Network (multi-encoding, DOUBLE width)

Encodings:
  Node size     = severe cliff degree
  Node color    = compound consensus pIC50 (viridis colormap)
  Node border   = thick navy for Hub A, thick red for Hub B, thin grey for others
  Edge color    = ΔpIC50 (light grey → #AA0044 dark magenta)
  Inset panel b = top-12 degree bar chart

Outputs: publication/figures/main/fig05_cliff_network.{png,pdf}
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
import matplotlib.lines as mlines
from matplotlib.colors import Normalize, to_rgba
from matplotlib import cm
from pathlib import Path

set_style()

try:
    import networkx as nx
    HAS_NX = True
except ImportError:
    HAS_NX = False
    print("WARNING: networkx not available — bar chart fallback")

CANON  = {'n_compounds': 3093, 'n_severe': 94, 'n_in_severe': 99}

HUB_IKS = {
    'A1': 'SMADULGDNOCLOP-GISFHXKWSA-N',
    'A2': 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
    'B1': 'UDCDEKJNAMHBFH-HSZRJFAPSA-N',
    'B2': 'DVCKJOQIVOGXEI-XMMPIXPASA-N',
}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT  = ROOT / 'publication/figures/main'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FIGURE 5 — CLIFF NETWORK (multi-encoding)")
print("=" * 60)

df    = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
cliffs = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
disc  = pd.read_csv(ROOT / 'outputs/mmp/mmp_discontinuity_scores.csv')

sev   = cliffs[cliffs['cliff_tier'] == 'severe'].copy()
assert len(sev) == CANON['n_severe']
sev_iks = set(sev['inchi_key_a'].tolist() + sev['inchi_key_b'].tolist())
assert len(sev_iks) == CANON['n_in_severe']

pic50_map = dict(zip(df['inchi_key'], df['pIC50']))
hub_a_set = {HUB_IKS['A1'], HUB_IKS['A2']}
hub_b_set = {HUB_IKS['B1'], HUB_IKS['B2']}

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(DOUBLE, 3.6),
                                   gridspec_kw={'width_ratios': [1.7, 1.0]},
                                   constrained_layout=True)

# ── Panel a: Network ──────────────────────────────────────────────────────────
if HAS_NX:
    G = nx.Graph()
    for _, row in sev.iterrows():
        G.add_edge(row['inchi_key_a'], row['inchi_key_b'],
                   delta=abs(row['delta_pic50']))

    nodes = list(G.nodes())
    n     = len(nodes)
    print(f"  Graph: {n} nodes, {G.number_of_edges()} edges")

    # Spring layout with k=0.25 (per spec)
    pos = nx.spring_layout(G, seed=42, k=0.25)

    # Manual offset for near-identical node pairs (Hub B Tan=0.975, Hub A Tan=0.761)
    for pair in [('B1', 'B2'), ('A1', 'A2')]:
        ik0, ik1 = HUB_IKS[pair[0]], HUB_IKS[pair[1]]
        if ik0 in pos and ik1 in pos:
            p0, p1 = np.array(pos[ik0]), np.array(pos[ik1])
            if np.linalg.norm(p0 - p1) < 0.10:
                mid = (p0 + p1) / 2
                off = np.array([0.06, 0.06])
                pos[ik0] = list(mid + off)
                pos[ik1] = list(mid - off)

    # ── Edges coloured by ΔpIC50 ─────────────────────────────────────────────
    edge_cmap = cm.get_cmap('Reds')
    delta_vals = np.array([d['delta'] for _, _, d in G.edges(data=True)])
    delta_min, delta_max = delta_vals.min(), delta_vals.max()
    edge_norm  = Normalize(vmin=delta_min, vmax=delta_max)

    for (u, v, data) in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        t = edge_norm(data['delta'])
        # Blend: grey (#888888) → dark magenta (#AA0044)
        # Low ΔpIC50 → medium grey; high ΔpIC50 → dark magenta, fully opaque
        edge_rgba = (0.53 + t * (0.67 - 0.53),
                     0.53 + t * (0.00 - 0.53),
                     0.53 + t * (0.27 - 0.53),
                     0.40 + t * 0.40)   # alpha 0.40 (low ΔpIC50) -> 0.80 (high)
        ax_a.plot([x0, x1], [y0, y1],
                  color=edge_rgba, lw=0.7, zorder=1, rasterized=True)

    # ── Nodes coloured by pIC50 (viridis) ────────────────────────────────────
    node_cmap = cm.viridis
    pic50_min, pic50_max = 2.0, 8.52
    node_norm  = Normalize(vmin=pic50_min, vmax=pic50_max)

    # Node sizes proportional to degree
    degrees = dict(G.degree())
    min_deg, max_deg = min(degrees.values()), max(degrees.values())
    def node_size(ik):
        d = degrees.get(ik, 1)
        # scale: min=10, max=140 — keeps hub markers legible without dominating the plot
        return 10 + (d - min_deg) / max(1, max_deg - min_deg) * 130

    # Draw non-hub nodes
    hub_all = hub_a_set | hub_b_set
    other_nodes = [nd for nd in nodes if nd not in hub_all]
    if other_nodes:
        ox  = np.array([pos[nd][0] for nd in other_nodes])
        oy  = np.array([pos[nd][1] for nd in other_nodes])
        oz  = np.array([pic50_map.get(nd, 6.0) for nd in other_nodes])
        osz = np.array([node_size(nd) for nd in other_nodes])
        sc_other = ax_a.scatter(ox, oy, s=osz, c=oz,
                                 cmap=node_cmap, norm=node_norm,
                                 edgecolors='#888888', linewidths=0.4,
                                 zorder=3, rasterized=True)

    # Hub A (navy border, thick) — fixed size so star doesn't overwhelm
    HUB_A_SZ = 180   # star marker; visually matches ~circle s≈80
    for ik in hub_a_set:
        if ik in pos:
            ax_a.scatter([pos[ik][0]], [pos[ik][1]],
                          s=HUB_A_SZ, marker='*',
                          c=[[pic50_map.get(ik, 5.4)]],
                          cmap=node_cmap, norm=node_norm,
                          edgecolors=SEM['classA'], linewidths=2.0, zorder=6)

    # Hub B (red border, thick) — diamond at fixed size; smaller than star visually
    HUB_B_SZ = 52    # diamond is visually larger per point than circle → keep small
    for ik in hub_b_set:
        if ik in pos:
            ax_a.scatter([pos[ik][0]], [pos[ik][1]],
                          s=HUB_B_SZ, marker='D',
                          c=[[pic50_map.get(ik, 4.3)]],
                          cmap=node_cmap, norm=node_norm,
                          edgecolors=SEM['classB'], linewidths=1.5, zorder=6)

    # pIC50 colorbar
    sm_nodes = cm.ScalarMappable(cmap=node_cmap, norm=node_norm)
    sm_nodes.set_array([])
    cbar = fig.colorbar(sm_nodes, ax=ax_a, shrink=0.40, aspect=14,
                        pad=0.02, orientation='vertical')
    cbar.set_label('pIC50', fontsize=6)
    cbar.ax.tick_params(labelsize=5.5)

    # Hub labels
    for lbl, ik in HUB_IKS.items():
        if ik in pos:
            x, y = pos[ik]
            color = SEM['classA'] if lbl.startswith('A') else SEM['classB']
            offsets = {'A1': (0.12, 0.10), 'A2': (-0.14, -0.10),
                       'B1': (0.12, 0.10), 'B2': (-0.14, -0.10)}
            xo, yo = offsets[lbl]
            ax_a.annotate(lbl, xy=(x, y), xytext=(x + xo, y + yo),
                          fontsize=6, fontweight='bold', color=color,
                          arrowprops=dict(arrowstyle='-', lw=0.4, color=color))

    ax_a.set_title(f'{n} nodes · {G.number_of_edges()} severe cliff edges', fontsize=6, pad=4)
    n_nodes = n
    n_edges = G.number_of_edges()

else:
    ax_a.text(0.5, 0.5, 'networkx not available', ha='center', va='center',
              transform=ax_a.transAxes)
    n_nodes = n_edges = 'N/A'

for sp in ['top', 'right', 'left', 'bottom']:
    ax_a.spines[sp].set_visible(False)
ax_a.set_xticks([])
ax_a.set_yticks([])

# Split legend: NODE encoding (markers/colour) vs EDGE encoding (ΔpIC50)
node_handles = [
    mlines.Line2D([], [], color=SEM['classA'], marker='*', linestyle='None',
                  markersize=8, label='Hub A — series floor (n=2)'),
    mlines.Line2D([], [], color=SEM['classB'], marker='D', linestyle='None',
                  markersize=6, label='Hub B — singleton attractor (n=2)'),
    mlines.Line2D([], [], color='#888888', marker='o', linestyle='None',
                  markersize=5, alpha=0.7, label='Non-hub (size scales w/ degree)'),
]
leg_node = ax_a.legend(handles=node_handles, fontsize=5.5, loc='lower left',
                       title='Nodes (colour = pIC50)', framealpha=0.90,
                       edgecolor='none', handletextpad=0.4, labelspacing=0.35)
leg_node.get_title().set_fontsize(5.5)
leg_node.get_title().set_fontweight('bold')
ax_a.add_artist(leg_node)

# Edge-encoding legend (grey → dark magenta with ΔpIC50)
edge_handles = [
    mlines.Line2D([], [], color='#888888', lw=1.2, label=f'low (≈{delta_min:.1f})'),
    mlines.Line2D([], [], color=SEM['cliff'], lw=1.2, label=f'high (≈{delta_max:.1f})'),
]
leg_edge = ax_a.legend(handles=edge_handles, fontsize=5.5, loc='lower right',
                       title='Edge ΔpIC50', framealpha=0.90,
                       edgecolor='none', handletextpad=0.4, labelspacing=0.35)
leg_edge.get_title().set_fontsize(5.5)
leg_edge.get_title().set_fontweight('bold')
panel_label(ax_a, 'a', x=-0.02, y=1.04)

# ── Panel b: Top-12 degree bar chart ─────────────────────────────────────────
print("[Panel b] Degree bar chart ...")
disc_top = disc.sort_values('severe_cliff_degree', ascending=False).head(12)

hub_label_map = {
    HUB_IKS['A1']: 'Hub A1', HUB_IKS['A2']: 'Hub A2',
    HUB_IKS['B1']: 'Hub B1', HUB_IKS['B2']: 'Hub B2',
}
bar_colors_b = []
y_labels_b   = []
for i, (_, row) in enumerate(disc_top.iterrows()):
    ik = row['inchi_key']
    if ik in hub_a_set:
        bar_colors_b.append(SEM['classA'])
    elif ik in hub_b_set:
        bar_colors_b.append(SEM['classB'])
    else:
        bar_colors_b.append(C['grey'])
    # Hubs get named labels; non-hubs show their InChIKey skeleton (first block) for traceability
    y_labels_b.append(hub_label_map.get(ik, ik.split('-')[0]))

y_pos = list(range(len(disc_top)))
ax_b.barh(y_pos, disc_top['severe_cliff_degree'].values,
           color=bar_colors_b, height=0.58, edgecolor='white', lw=0.3)

for i, val in enumerate(disc_top['severe_cliff_degree']):
    ax_b.text(val + 0.15, i, str(int(val)), va='center', ha='left', fontsize=6)

ax_b.set_yticks(y_pos)
ax_b.set_yticklabels(y_labels_b, fontsize=6)
ax_b.set_xlabel('Severe cliff partners (network degree)')
ax_b.invert_yaxis()
panel_label(ax_b, 'b', x=-0.22, y=1.04)

# ── Save ──────────────────────────────────────────────────────────────────────
save_fig(fig, str(OUT / 'fig05_cliff_network'))
plt.close(fig)
print(f"Figure 5 complete. Nodes={n_nodes}, Edges={n_edges}")
