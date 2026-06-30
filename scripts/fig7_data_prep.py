#!/usr/bin/env python
import os
import sys
os.chdir('/home/nidhal/PAD4-db_V2')

import pandas as pd
import numpy as np

os.makedirs('outputs/figures', exist_ok=True)

HUB_CLASS_A = {'SMADULGDNOCLOP-GISFHXKWSA-N', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N'}
HUB_CLASS_B = {'UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'DVCKJOQIVOGXEI-XMMPIXPASA-N'}

# ── Load data ───────────────────────────────────────────────────────────────
cliffs = pd.read_parquet('data/processed/activity_cliffs.parquet')
comps = pd.read_parquet('data/processed/pad4_compounds.parquet')

# ── Severe cliff pairs (edges) ──────────────────────────────────────────────
severe = cliffs[cliffs['cliff_tier'] == 'severe'].copy().reset_index(drop=True)

# Enrich edges with mechanism info from pad4_compounds
mech_map = comps.set_index('inchi_key')['mechanism_class'].to_dict()
severe['mech_a'] = severe['inchi_key_a'].map(mech_map)
severe['mech_b'] = severe['inchi_key_b'].map(mech_map)
severe['same_mechanism'] = severe['mech_a'] == severe['mech_b']

def hub_class_of(ik):
    if ik in HUB_CLASS_A:
        return 'A'
    if ik in HUB_CLASS_B:
        return 'B'
    return 'none'

severe['hub_class_a'] = severe['inchi_key_a'].map(hub_class_of)
severe['hub_class_b'] = severe['inchi_key_b'].map(hub_class_of)

def combined_hub_class(row):
    ca, cb = row['hub_class_a'], row['hub_class_b']
    if ca == 'none' and cb == 'none':
        return 'none'
    if ca != 'none' and cb != 'none':
        if ca == cb:
            return ca
        return 'both'
    return ca if ca != 'none' else cb

severe['hub_class_involved'] = severe.apply(combined_hub_class, axis=1)
severe['hub_involved'] = severe['hub_class_involved'] != 'none'

edges_out = severe[[
    'inchi_key_a', 'inchi_key_b', 'tanimoto', 'delta_pic50',
    'mech_a', 'mech_b', 'same_mechanism', 'hub_involved', 'hub_class_involved',
]].copy()

# ── Severe cliff nodes ──────────────────────────────────────────────────────
all_iks = set(severe['inchi_key_a']) | set(severe['inchi_key_b'])

# Severe cliff degree per compound
deg_a = severe['inchi_key_a'].value_counts()
deg_b = severe['inchi_key_b'].value_counts()
degree = (deg_a.add(deg_b, fill_value=0)).astype(int)

node_data = comps[comps['inchi_key'].isin(all_iks)].copy()
node_data = node_data.set_index('inchi_key')
node_data['severe_cliff_degree'] = degree
node_data['hub_class'] = node_data.index.map(hub_class_of)

nodes_out = node_data[[
    'smiles_std', 'pic50_consensus', 'mechanism_class',
    'is_covalent', 'fragment_flag', 'severe_cliff_degree', 'hub_class',
]].copy().reset_index().rename(columns={'index': 'inchi_key'})
nodes_out = nodes_out.sort_values('severe_cliff_degree', ascending=False).reset_index(drop=True)

# ── Verification ────────────────────────────────────────────────────────────
n_nodes = len(nodes_out)
n_edges = len(edges_out)
hub_a_nodes = nodes_out['hub_class'].eq('A').sum()
hub_b_nodes = nodes_out['hub_class'].eq('B').sum()
hub_a_edges = severe['hub_class_involved'].eq('A').sum()
hub_b_edges = severe['hub_class_involved'].eq('B').sum()
non_hub_edges = severe['hub_class_involved'].eq('none').sum()
cross_mech_edges = (~severe['same_mechanism']).sum()

hub_and_cross = int((severe['hub_involved'] & ~severe['same_mechanism']).sum())

# Sanity checks (data-derived invariants)
assert n_edges == 94, f'Expected 94 severe pairs, got {n_edges}'
assert int(hub_a_nodes) == 2, f'Expected 2 Hub A nodes, got {hub_a_nodes}'
assert int(hub_b_nodes) == 2, f'Expected 2 Hub B nodes, got {hub_b_nodes}'
assert int(hub_a_edges) + int(hub_b_edges) + int(non_hub_edges) == 94, 'Edge counts do not sum to 94'

# ── Write outputs ────────────────────────────────────────────────────────────
nodes_out.to_csv('outputs/figures/fig7_nodes.csv', index=False)
edges_out.to_csv('outputs/figures/fig7_edges.csv', index=False)
print(f'Written: outputs/figures/fig7_nodes.csv  ({n_nodes} rows)')
print(f'Written: outputs/figures/fig7_edges.csv  ({n_edges} rows)')

# ── Print summary ────────────────────────────────────────────────────────────
print()
print('=== Final Summary ===')
print(f'  Nodes:                              {n_nodes}')
print(f'  Edges:                              {n_edges}')
print(f'  Hub A edges:                        {int(hub_a_edges)}')
print(f'  Hub B edges:                        {int(hub_b_edges)}')
print(f'  Non-hub edges:                      {int(non_hub_edges)}  (expected 44)')
print(f'  Cross-mechanism edges:              {int(cross_mech_edges)}')
print(f'  Hub-involved AND cross-mechanism:   {hub_and_cross}')
