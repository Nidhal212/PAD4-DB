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
    'purple':     '#AA4499',
}

CANON = {
    'n_compounds': 3093,
    'n_severe': 94,
    'n_in_severe': 99,
}

HUB_IKS = {
    'A1': 'SMADULGDNOCLOP-GISFHXKWSA-N',
    'A2': 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
    'B1': 'UDCDEKJNAMHBFH-HSZRJFAPSA-N',
    'B2': 'DVCKJOQIVOGXEI-XMMPIXPASA-N',
}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/figures/nature_v2'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FIGURE 5 — CLIFF NETWORK")
print("=" * 60)

df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
cliffs = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
disc = pd.read_csv(ROOT / 'outputs/mmp/mmp_discontinuity_scores.csv')

# Get severe cliffs
severe = cliffs[cliffs['cliff_tier'] == 'severe'].copy()
print(f"  Severe pairs: {len(severe)} (canonical: {CANON['n_severe']})")
assert len(severe) == CANON['n_severe'], f"n_severe mismatch: {len(severe)}"

severe_iks = set(severe['inchi_key_a'].tolist() + severe['inchi_key_b'].tolist())
print(f"  Severe nodes: {len(severe_iks)} (canonical: {CANON['n_in_severe']})")
assert len(severe_iks) == CANON['n_in_severe'], f"n_in_severe: {len(severe_iks)}"

# Map mechanism to nodes
mech_map = dict(zip(df['inchi_key'], df['mechanism_class']))

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11, 6),
                                   gridspec_kw={'width_ratios': [1.5, 1]},
                                   constrained_layout=True)

# ── Panel a: Network graph ─────────────────────────────────────────────────────
print("\n[Panel a] Network graph ...")

if HAS_NX:
    G = nx.Graph()
    for _, row in severe.iterrows():
        G.add_edge(row['inchi_key_a'], row['inchi_key_b'],
                   delta=abs(row['delta_pic50']),
                   mech_a=mech_map.get(row['inchi_key_a'], ''),
                   mech_b=mech_map.get(row['inchi_key_b'], ''))

    nodes = list(G.nodes())
    print(f"  Graph: {len(nodes)} nodes, {G.number_of_edges()} edges")

    pos = nx.spring_layout(G, seed=42, k=2.5)

    hub_a_nodes = [HUB_IKS['A1'], HUB_IKS['A2']]
    hub_b_nodes = [HUB_IKS['B1'], HUB_IKS['B2']]
    hub_all = set(hub_a_nodes + hub_b_nodes)

    # Classify nodes
    conf_nodes = [n for n in nodes if mech_map.get(n, '') == 'enzymatic_confirmed' and n not in hub_all]
    other_nodes = [n for n in nodes if n not in hub_all and n not in conf_nodes]
    hub_a_in = [n for n in hub_a_nodes if n in G]
    hub_b_in = [n for n in hub_b_nodes if n in G]

    # Draw edges
    for u, v, data in G.edges(data=True):
        mech_u = mech_map.get(u, '')
        mech_v = mech_map.get(v, '')
        x0, y0 = pos[u]
        x1, y1 = pos[v]

        is_hub_a = u in hub_a_nodes or v in hub_a_nodes
        is_hub_b = u in hub_b_nodes or v in hub_b_nodes
        cross_mech = (mech_u != mech_v and mech_u and mech_v)

        if is_hub_a and not is_hub_b:
            color, lw, alpha, ls = COLORS['navy'], 0.8, 0.4, 'solid'
        elif is_hub_b and not is_hub_a:
            color, lw, alpha, ls = COLORS['red'], 0.8, 0.4, 'solid'
        elif cross_mech:
            color, lw, alpha, ls = COLORS['purple'], 1.0, 0.6, 'dashed'
        else:
            color, lw, alpha, ls = COLORS['grey'], 0.3, 0.15, 'solid'

        ax_a.plot([x0, x1], [y0, y1], color=color, lw=lw, alpha=alpha, ls=ls, zorder=1)

    # Draw nodes
    if other_nodes:
        ox = [pos[n][0] for n in other_nodes]
        oy = [pos[n][1] for n in other_nodes]
        ax_a.scatter(ox, oy, s=15, c=COLORS['grey'], alpha=0.5, zorder=2)

    if conf_nodes:
        cx = [pos[n][0] for n in conf_nodes]
        cy = [pos[n][1] for n in conf_nodes]
        ax_a.scatter(cx, cy, s=20, c=COLORS['teal'], alpha=0.7, zorder=3)

    # Hub A
    for ik in hub_a_in:
        ax_a.scatter([pos[ik][0]], [pos[ik][1]], s=400, c=COLORS['navy'],
                      marker='*', edgecolors='white', lw=1.5, zorder=5)

    # Hub B
    for ik in hub_b_in:
        ax_a.scatter([pos[ik][0]], [pos[ik][1]], s=300, c=COLORS['red'],
                      marker='D', edgecolors='white', lw=1.5, zorder=5)

    # Labels
    labels = {'A1': HUB_IKS['A1'], 'A2': HUB_IKS['A2'],
               'B1': HUB_IKS['B1'], 'B2': HUB_IKS['B2']}
    for lbl, ik in labels.items():
        if ik in pos:
            x, y = pos[ik]
            color = COLORS['navy'] if lbl.startswith('A') else COLORS['red']
            ax_a.text(x + 0.05, y + 0.05, lbl, fontsize=10, fontweight='bold', color=color)

    # Stats box
    # Get degrees from disc scores
    deg_map = {}
    for _, row in disc.iterrows():
        deg_map[row['inchi_key']] = row['severe_cliff_degree']

    a1_deg = deg_map.get(HUB_IKS['A1'], '?')
    a2_deg = deg_map.get(HUB_IKS['A2'], '?')
    b1_deg = deg_map.get(HUB_IKS['B1'], '?')
    b2_deg = deg_map.get(HUB_IKS['B2'], '?')

    stats_txt = (f'Nodes: {len(nodes)}\nEdges: {G.number_of_edges()}\n'
                 f'Hub A1 degree: {a1_deg}\nHub A2 degree: {a2_deg}\n'
                 f'Hub B1 degree: {b1_deg}\nHub B2 degree: {b2_deg}')
    ax_a.text(0.02, 0.98, stats_txt, transform=ax_a.transAxes,
              fontsize=7.5, ha='left', va='top', family='monospace',
              bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                        edgecolor=COLORS['grey'], lw=0.8))

else:
    ax_a.text(0.5, 0.5, 'networkx not available\n(placeholder)',
              transform=ax_a.transAxes, ha='center', va='center', fontsize=10)

ax_a.set_xticks([])
ax_a.set_yticks([])
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.spines['left'].set_visible(False)
ax_a.spines['bottom'].set_visible(False)

# Legend
leg_handles = [
    plt.scatter([], [], s=80, c=COLORS['navy'], marker='*', label='Hub A'),
    plt.scatter([], [], s=60, c=COLORS['red'], marker='D', label='Hub B'),
    mpatches.Patch(color=COLORS['teal'], label='Enzymatic confirmed', alpha=0.7),
    mpatches.Patch(color=COLORS['grey'], label='Other', alpha=0.5),
]
ax_a.legend(handles=leg_handles, fontsize=7, loc='lower right', framealpha=0.7, edgecolor='none')
ax_a.text(-0.05, 1.02, 'a', transform=ax_a.transAxes, fontsize=11, fontweight='bold')

# ── Panel b: Degree bar chart (top 12) ────────────────────────────────────────
print("\n[Panel b] Degree bar chart ...")
disc_top = disc.sort_values('severe_cliff_degree', ascending=False).head(12)

bar_colors_b = []
for ik in disc_top['inchi_key']:
    if ik in {HUB_IKS['A1'], HUB_IKS['A2']}:
        bar_colors_b.append(COLORS['navy'])
    elif ik in {HUB_IKS['B1'], HUB_IKS['B2']}:
        bar_colors_b.append(COLORS['red'])
    else:
        bar_colors_b.append(COLORS['grey'])

y_labels = [ik[:14] + '...' for ik in disc_top['inchi_key']]
y_pos = range(len(disc_top))

ax_b.barh(list(y_pos), disc_top['severe_cliff_degree'].values,
           color=bar_colors_b, height=0.6, edgecolor='white', lw=0.3)

for i, (val, ik) in enumerate(zip(disc_top['severe_cliff_degree'], disc_top['inchi_key'])):
    ax_b.text(val + 0.1, i, str(int(val)), va='center', ha='left', fontsize=7.5)

ax_b.set_yticks(list(y_pos))
ax_b.set_yticklabels(y_labels, fontsize=6.5, family='monospace')
ax_b.set_xlabel('Number of severe cliff partners', fontsize=9)
ax_b.invert_yaxis()
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.text(-0.2, 1.02, 'b', transform=ax_b.transAxes, fontsize=11, fontweight='bold')

# ── Save ──────────────────────────────────────────────────────────────────────
for ext in ['png', 'pdf']:
    outpath = OUT / f'fig05_cliff_network.{ext}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"  Saved: {outpath}")

plt.close(fig)
print(f"\nFigure 5 complete. Nodes={len(nodes) if HAS_NX else 'N/A'}, Edges={G.number_of_edges() if HAS_NX else 'N/A'}")
