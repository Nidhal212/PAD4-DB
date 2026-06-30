#!/usr/bin/env python
"""Figure 7 — Activity Cliff Network (1 large + 1 bar panel) + Great Tables."""
import os, sys, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

# ── Library standards ─────────────────────────────────────────────────────────
import scienceplots
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
plt.style.use(['science', 'nature', 'no-latex'])

from great_tables import GT, loc, style as gt_style
import great_tables
import networkx as nx

import numpy as np
import pandas as pd

os.makedirs('outputs/figures', exist_ok=True)
os.makedirs('outputs/tables', exist_ok=True)

print(f"SciencePlots: importable")
print(f"Great Tables:  {great_tables.__version__}")
print(f"NetworkX:      {nx.__version__}")
print()

# ── Load data ─────────────────────────────────────────────────────────────────
nodes_df = pd.read_csv('outputs/figures/fig7_nodes.csv')
edges_df  = pd.read_csv('outputs/figures/fig7_edges.csv')

# ── Locked number verification ────────────────────────────────────────────────
n_nodes   = len(nodes_df)
n_edges   = len(edges_df)
hub_a_n   = int((nodes_df['hub_class'] == 'A').sum())
hub_b_n   = int((nodes_df['hub_class'] == 'B').sum())
hub_a_e   = int((edges_df['hub_class_involved'] == 'A').sum())
hub_b_e   = int((edges_df['hub_class_involved'] == 'B').sum())
nonhub_e  = int((edges_df['hub_class_involved'] == 'none').sum())
cross_e   = int((~edges_df['same_mechanism']).sum())

print("=== Locked number verification ===")
LOCKED = [
    ("Nodes (severe cliff compounds)", n_nodes,  99),
    ("Edges (severe cliff pairs)",     n_edges,  94),
    ("Hub Class A nodes",              hub_a_n,   2),
    ("Hub Class B nodes",              hub_b_n,   2),
    ("Hub A edges",                    hub_a_e,  27),
    ("Hub B edges",                    hub_b_e,  23),
    ("Non-hub edges",                  nonhub_e, 44),
    ("Cross-mechanism edges",          cross_e,   4),
]
all_pass = True
for label, actual, expected in LOCKED:
    ok = actual == expected
    print(f"  {label:35s}: {actual}  {'PASS' if ok else f'FAIL (expected {expected})'}")
    if not ok: all_pass = False
if not all_pass:
    print("\nVERIFICATION FAILED — stopping.")
    sys.exit(1)
print()

# ── Build networkx graph ──────────────────────────────────────────────────────
G = nx.Graph()
for _, row in nodes_df.iterrows():
    G.add_node(row['inchi_key'], **row.to_dict())
for _, row in edges_df.iterrows():
    G.add_edge(row['inchi_key_a'], row['inchi_key_b'], **row.to_dict())

pos = nx.spring_layout(G, seed=42, k=2.5)

# ── Node styling ──────────────────────────────────────────────────────────────
HUB_A_IKS = {'SMADULGDNOCLOP-GISFHXKWSA-N', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N'}
HUB_B_IKS = {'UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'DVCKJOQIVOGXEI-XMMPIXPASA-N'}
MECH_COLORS = {
    'enzymatic':           '#AAAAAA',
    'enzymatic_confirmed': '#4A90D9',
    'fp_ic50':             '#2ECC71',
    'covalent':            '#E05A2B',
}
node_meta = nodes_df.set_index('inchi_key')

node_list  = list(G.nodes())
node_sizes  = []
node_colors = []
for ik in node_list:
    deg = int(node_meta.loc[ik, 'severe_cliff_degree'])
    node_sizes.append(30 + deg * 25)
    hc = node_meta.loc[ik, 'hub_class']
    if hc == 'A':
        node_colors.append('#E74C3C')
    elif hc == 'B':
        node_colors.append('#1A237E')
    else:
        mech = node_meta.loc[ik, 'mechanism_class']
        node_colors.append(MECH_COLORS.get(mech, '#AAAAAA'))

# ── Edge sets ─────────────────────────────────────────────────────────────────
# Build lookup both directions (store as dicts to avoid Series truth-value ambiguity)
edge_lookup = {}
for _, row in edges_df.iterrows():
    d = row.to_dict()
    edge_lookup[(row['inchi_key_a'], row['inchi_key_b'])] = d
    edge_lookup[(row['inchi_key_b'], row['inchi_key_a'])] = d

def get_edge_row(u, v):
    return edge_lookup.get((u, v)) or edge_lookup.get((v, u))

hub_a_edges   = [(u, v) for u, v in G.edges() if get_edge_row(u,v)['hub_class_involved'] == 'A']
hub_b_edges   = [(u, v) for u, v in G.edges() if get_edge_row(u,v)['hub_class_involved'] == 'B']
nonhub_edges  = [(u, v) for u, v in G.edges() if get_edge_row(u,v)['hub_class_involved'] == 'none']
cross_edges   = [(u, v) for u, v in G.edges() if not get_edge_row(u,v)['same_mechanism']]
cross_set     = set(map(frozenset, cross_edges))

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 7))
ax_a = fig.add_axes([0.01, 0.05, 0.58, 0.90])   # left 60%
ax_b = fig.add_axes([0.63, 0.08, 0.35, 0.86])   # right 35%

# ── Panel A: Network ──────────────────────────────────────────────────────────
ax_a.set_axis_off()

def draw_edge_set(edge_list, color, alpha, lw, style='solid'):
    # Separate into cross-mechanism (dashed) and non-cross (solid)
    solid_edges = [e for e in edge_list if frozenset(e) not in cross_set]
    dashed_edges = [e for e in edge_list if frozenset(e) in cross_set]
    if solid_edges:
        nx.draw_networkx_edges(G, pos, edgelist=solid_edges, ax=ax_a,
                               edge_color=color, alpha=alpha, width=lw,
                               style='solid')
    if dashed_edges:
        nx.draw_networkx_edges(G, pos, edgelist=dashed_edges, ax=ax_a,
                               edge_color=color, alpha=alpha, width=lw,
                               style='dashed')

# Draw edges: non-hub first (background), then hub A, hub B
draw_edge_set(nonhub_edges, '#AAAAAA', 0.4, 0.8)
draw_edge_set(hub_a_edges,  '#E74C3C', 0.6, 1.5)
draw_edge_set(hub_b_edges,  '#1A237E', 0.6, 1.5)

# Draw nodes
nx.draw_networkx_nodes(G, pos, nodelist=node_list, node_size=node_sizes,
                       node_color=node_colors, edgecolors='white',
                       linewidths=0.5, ax=ax_a)

# Label only the 4 hub nodes
hub_labels = {ik: ik[:14] for ik in (HUB_A_IKS | HUB_B_IKS) if ik in G.nodes()}
nx.draw_networkx_labels(G, pos, labels=hub_labels, ax=ax_a,
                        font_size=6, font_color='black')

ax_a.set_title("Severe Activity Cliff Network\n(94 pairs, 99 compounds)",
               fontsize=9, pad=4)

# Legend
legend_handles = [
    mpatches.Patch(color='#E74C3C', label='Hub Class A'),
    mpatches.Patch(color='#1A237E', label='Hub Class B'),
    mpatches.Patch(color='#4A90D9', label='Enzymatic confirmed'),
    mpatches.Patch(color='#AAAAAA', label='Enzymatic'),
    mlines.Line2D([], [], color='#E74C3C', linewidth=1.5, label='Hub A cliff pair'),
    mlines.Line2D([], [], color='#1A237E', linewidth=1.5, label='Hub B cliff pair'),
    mlines.Line2D([], [], color='#AAAAAA', linewidth=0.8, label='Non-hub cliff pair'),
    mlines.Line2D([], [], color='gray',    linewidth=0.8, linestyle='--',
                  label='Cross-mechanism pair'),
]
ax_a.legend(handles=legend_handles, loc='lower left', fontsize=7,
            framealpha=0.85, ncol=2)

# ── Panel B: Degree distribution (top 20) ─────────────────────────────────────
deg_df = nodes_df[['inchi_key', 'severe_cliff_degree', 'hub_class']].copy()
deg_df = deg_df.sort_values('severe_cliff_degree', ascending=False).head(20)
deg_df['short_ik'] = deg_df['inchi_key'].str[:14]

bar_colors = []
for _, row in deg_df.iterrows():
    if row['hub_class'] == 'A':
        bar_colors.append('#E74C3C')
    elif row['hub_class'] == 'B':
        bar_colors.append('#1A237E')
    else:
        bar_colors.append('#AAAAAA')

y_pos = np.arange(len(deg_df))
bars = ax_b.barh(y_pos, deg_df['severe_cliff_degree'].values,
                 color=bar_colors, height=0.7,
                 edgecolor='white', linewidth=0.3)
ax_b.set_yticks(y_pos)
ax_b.set_yticklabels(deg_df['short_ik'].values, fontsize=6)
ax_b.invert_yaxis()
ax_b.set_xlabel('Severe Cliff Degree', fontsize=9)
ax_b.set_title('Top 20 Nodes by Degree', fontsize=9)
ax_b.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
ax_b.axvline(5, color='gray', linestyle='--', linewidth=0.8, alpha=0.5,
             label='degree=5')
ax_b.legend(loc='lower right', fontsize=7, framealpha=0.8)

# Degree count annotations
xmax = deg_df['severe_cliff_degree'].max()
for i, (_, row) in enumerate(deg_df.iterrows()):
    d = row['severe_cliff_degree']
    ax_b.text(d + 0.1, i, str(d), va='center', fontsize=6, color='#333333')
ax_b.set_xlim(0, xmax + 1.8)

plt.savefig('outputs/figures/fig7_cliff_network.png', dpi=300, bbox_inches='tight')
plt.savefig('outputs/figures/fig7_cliff_network.svg', bbox_inches='tight')
plt.close()

PNG_PATH = 'outputs/figures/fig7_cliff_network.png'
SVG_PATH = 'outputs/figures/fig7_cliff_network.svg'
print(f"Saved: {PNG_PATH}  ({os.path.getsize(PNG_PATH)/1024:.1f} KB)")
print(f"Saved: {SVG_PATH}  ({os.path.getsize(SVG_PATH)/1024:.1f} KB)")
print()

# ── Great Tables: hub compound summary ───────────────────────────────────────
hub_data = pd.DataFrame([
    {
        'compound_id':    'SMADULGDNOCLOP',
        'hub_class':      'A',
        'pic50':          5.390,
        'severe_degree':  15,
        'pct_of_94_pairs': 15.96,
        'scaffold_type':  'series-embedded',
        'mechanism_class': 'enzymatic',
        'sources':        '3-source',
    },
    {
        'compound_id':    'RAVBZQAQTVGKIV',
        'hub_class':      'A',
        'pic50':          5.341,
        'severe_degree':  12,
        'pct_of_94_pairs': 12.77,
        'scaffold_type':  'series-embedded',
        'mechanism_class': 'enzymatic',
        'sources':        '3-source',
    },
    {
        'compound_id':    'UDCDEKJNAMHBFH',
        'hub_class':      'B',
        'pic50':          4.301,
        'severe_degree':  12,
        'pct_of_94_pairs': 12.77,
        'scaffold_type':  'singleton',
        'mechanism_class': 'enzymatic',
        'sources':        '2-source',
    },
    {
        'compound_id':    'DVCKJOQIVOGXEI',
        'hub_class':      'B',
        'pic50':          4.301,
        'severe_degree':  11,
        'pct_of_94_pairs': 11.70,
        'scaffold_type':  'singleton',
        'mechanism_class': 'enzymatic',
        'sources':        '2-source',
    },
])

gt = (
    GT(hub_data)
    .tab_header(
        title="PAD4-DB Activity Cliff Hub Compounds",
        subtitle="Four compounds accounting for 53.2% of all severe cliff pairs",
    )
    .cols_label(
        compound_id="Compound",
        hub_class="Hub Class",
        pic50="pIC50",
        severe_degree="Cliff Pairs",
        pct_of_94_pairs="% of 94 Pairs",
        scaffold_type="Scaffold Type",
        mechanism_class="Mechanism",
        sources="Sources",
    )
    .fmt_number(columns=["pic50"], decimals=3)
    .fmt_number(columns=["pct_of_94_pairs"], decimals=1)
    .tab_style(
        style=gt_style.fill(color="#FFEBEE"),
        locations=loc.body(rows=[0, 1]),
    )
    .tab_style(
        style=gt_style.fill(color="#E8EAF6"),
        locations=loc.body(rows=[2, 3]),
    )
    .tab_style(
        style=gt_style.text(weight="bold"),
        locations=loc.body(columns=["hub_class"]),
    )
    .tab_source_note(
        "Severe cliff threshold: Tanimoto ≥ 0.8, ΔpIC50 ≥ 2.0. "
        "Hub Class A: mid-potency members of the 174-compound "
        "azaindole-benzimidazole scaffold series. "
        "Hub Class B: scaffold singletons with broad ECFP4 similarity."
    )
)

HTML_PATH = 'outputs/tables/fig7_cliff_stats.html'
TEX_PATH  = 'outputs/tables/fig7_cliff_stats.tex'
with open(HTML_PATH, 'w') as f:
    f.write(gt.as_raw_html())
with open(TEX_PATH, 'w') as f:
    f.write(gt.as_latex())
print(f"Great Tables HTML:  {HTML_PATH}  ({os.path.getsize(HTML_PATH)/1024:.1f} KB)")
print(f"Great Tables LaTeX: {TEX_PATH}  ({os.path.getsize(TEX_PATH)/1024:.1f} KB)")
print()

# ── Completion report ─────────────────────────────────────────────────────────
total_hub_edges = hub_a_e + hub_b_e
print("=== Completion Report ===")
print(f"Network: {n_nodes} nodes, {n_edges} edges ✓")
print(f"Hub A degree: SMADULGDNOCLOP=15, RAVBZQAQTVGKIV=12")
print(f"Hub B degree: UDCDEKJNAMHBFH=12, DVCKJOQIVOGXEI=11")
print(f"Collective hub coverage: {total_hub_edges}/{n_edges} = "
      f"{total_hub_edges/n_edges*100:.1f}%")
print(f"Non-hub edges: {nonhub_e}")
print(f"Cross-mechanism edges: {cross_e} (all dashed in network)")
print()
print("Files written:")
for p in [PNG_PATH, SVG_PATH, HTML_PATH, TEX_PATH]:
    print(f"  {p}  ({os.path.getsize(p)/1024:.1f} KB)")
