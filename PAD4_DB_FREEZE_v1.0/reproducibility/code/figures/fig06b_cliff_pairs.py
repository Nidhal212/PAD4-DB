"""
fig06b_cliff_pairs.py — Figure 6 (merged 6a + 6b), DOUBLE width.

Layout:
  Top 10%  — horizontal bar chart of MMP change types (was separate Fig 6a)
  Bottom 90% — 2×2 grid of 4 representative cliff pairs

Diff-atom coloring:
  Gain-of-function atoms (molecule with HIGHER pIC50) → Green (#00CC33)
  Loss-of-function atoms (molecule with LOWER  pIC50) → Red  (#CC0000)

Selection: 4 MMP-confirmed severe pairs (warhead-free, non-ecfp4-only)
  - Panels (i)+(ii): 2 highest-ΔpIC50 single_atom_change pairs
  - Panels (iii)+(iv): 2 highest-ΔpIC50 small_substituent pairs
  Each pair: ΔpIC50 = X.XX, Tan = 0.XX under structures.
  Connectivity pre-filter: both diff atom sets must be connected subgraphs.

Outputs: publication/figures/main/fig06b_cliff_pairs.{png,pdf}
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, save_fig, C, SEM, DOUBLE

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdFMCS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
import io
from PIL import Image

set_style()

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT  = ROOT / 'publication/figures/main'
OUT.mkdir(parents=True, exist_ok=True)

CLASS_A = {'SMADULGDNOCLOP-GISFHXKWSA-N', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N'}
CLASS_B = {'UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'DVCKJOQIVOGXEI-XMMPIXPASA-N'}


def hubclass(ik):
    return 'A' if ik in CLASS_A else ('B' if ik in CLASS_B else 'none')


def has_warhead(v):
    if v is None: return False
    if isinstance(v, float) and math.isnan(v): return False
    return str(v).lower() not in ('none', 'nan', '')


def ukey(a, b):
    return tuple(sorted([a, b]))


def _is_connected(mol, idxs):
    if len(idxs) <= 1: return True
    idxs = set(idxs)
    start = next(iter(idxs))
    visited = {start}
    queue = [start]
    while queue:
        cur = queue.pop()
        for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
            ni = nbr.GetIdx()
            if ni in idxs and ni not in visited:
                visited.add(ni)
                queue.append(ni)
    return visited == idxs


def compute_diff(mA, mB, shared_core_smarts):
    res = rdFMCS.FindMCS(
        [mA, mB], timeout=20,
        seedSmarts=shared_core_smarts,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        matchChiralTag=False,
    )
    if not res.smartsString: return None
    core = Chem.MolFromSmarts(res.smartsString)
    if core is None: return None
    mA_m = set(mA.GetSubstructMatch(core, useChirality=False))
    mB_m = set(mB.GetSubstructMatch(core, useChirality=False))
    if not mA_m or not mB_m: return None
    hlA = [i for i in range(mA.GetNumAtoms()) if i not in mA_m]
    hlB = [i for i in range(mB.GetNumAtoms()) if i not in mB_m]
    if not hlA or not hlB: return None
    if not _is_connected(mA, hlA) or not _is_connected(mB, hlB): return None
    if len(hlA) > 12 or len(hlB) > 12: return None
    return hlA, hlB, core


print("=" * 60)
print("FIG 6 — MMP BAR + REPRESENTATIVE CLIFF PAIRS")
print("=" * 60)

# ── Load data ────────────────────────────────────────────────────────────────
ac  = pd.read_parquet(ROOT / 'publication/data/activity_cliffs.parquet')
mmp = pd.read_csv(ROOT / 'outputs/mmp/mmp_pairs_cliff99.csv')
pad = pd.read_parquet(ROOT / 'publication/data/pad4_compounds.parquet')
wh  = dict(zip(pad.inchi_key, pad.warhead_class))
smi = dict(zip(pad.inchi_key, pad.smiles_std))
pic = dict(zip(pad.inchi_key, pad.pic50_consensus))

mmp_map = {
    ukey(r.inchi_key_a, r.inchi_key_b): {'mmp_type': r.mmp_type, 'shared_core': r.shared_core}
    for _, r in mmp.iterrows()
}

# ── MMP type counts for bar panel ────────────────────────────────────────────
# Count only MMP-confirmed SEVERE cliff pairs (not all MMP pairs among cliff compounds)
# Join: severe cliff pairs × mmp_pairs_cliff99 on canonical pair key
sev_ac = ac[ac['cliff_tier'] == 'severe'].copy()
sev_ac['pair_key'] = sev_ac.apply(
    lambda r: ukey(r['inchi_key_a'], r['inchi_key_b']), axis=1)
mmp_keys_in_sev = {ukey(r.inchi_key_a, r.inchi_key_b): r.mmp_type
                   for _, r in mmp.iterrows()}
sev_mmp_types = sev_ac['pair_key'].map(mmp_keys_in_sev).dropna()
mmp_type_counts = sev_mmp_types[sev_mmp_types.isin(
    ['single_atom_change', 'small_substituent', 'medium_substituent']
)].value_counts()
print(f"  MMP-confirmed severe cliff pairs by type: {mmp_type_counts.to_dict()}")
print(f"  Total MMP-confirmed: {mmp_type_counts.sum()}")

# ── Candidate pool ───────────────────────────────────────────────────────────
sev = ac[(ac.cliff_tier == 'severe') & (~ac.ecfp4_only_cliff)].copy()
rows = []
for _, r in sev.iterrows():
    k = ukey(r.inchi_key_a, r.inchi_key_b)
    if k not in mmp_map: continue
    a, b = r.inchi_key_a, r.inchi_key_b
    if has_warhead(wh.get(a)) or has_warhead(wh.get(b)): continue
    rows.append({'a': a, 'b': b, 'delta': abs(r.delta_pic50), 'tan': float(r.tanimoto),
                 'mmp_type': mmp_map[k]['mmp_type'],
                 'shared_core': mmp_map[k]['shared_core'],
                 'hub_a': hubclass(a), 'hub_b': hubclass(b)})

pool = pd.DataFrame(rows).sort_values('delta', ascending=False).reset_index(drop=True)
print(f"  Pool: {len(pool)} candidate pairs")

# Pre-check connectivity
viable = []
for _, r in pool.iterrows():
    a, b = r['a'], r['b']
    if a not in smi or b not in smi: continue
    mA = Chem.MolFromSmiles(smi[a])
    mB = Chem.MolFromSmiles(smi[b])
    if mA is None or mB is None: continue
    AllChem.Compute2DCoords(mA)
    result = compute_diff(mA, mB, r['shared_core'])
    if result is None: continue
    hlA, hlB, mcs_core = result
    viable.append({**r.to_dict(), 'hlA': hlA, 'hlB': hlB, 'mcs_core': mcs_core})

viable_df = pd.DataFrame(viable)
print(f"  Viable (connected diffs): {len(viable_df)} pairs")

chosen, used_cpds = [], set()

def add(row, why):
    a, b = row['a'], row['b']
    if a in used_cpds or b in used_cpds: return False
    used_cpds.add(a); used_cpds.add(b)
    chosen.append({**row, 'rule': why})
    return True

for _, row in viable_df[viable_df.mmp_type == 'single_atom_change'].iterrows():
    if add(row.to_dict(), f'{len(chosen)+1}_single_atom'):
        if sum(1 for c in chosen if c['mmp_type'] == 'single_atom_change') == 2: break

for _, row in viable_df[viable_df.mmp_type == 'small_substituent'].iterrows():
    if add(row.to_dict(), f'{len(chosen)+1}_small_subst'):
        if sum(1 for c in chosen if c['mmp_type'] == 'small_substituent') == 2: break

assert len(chosen) == 4, f"expected 4 pairs, got {len(chosen)}"

# Persist selection
sel = []
for c in chosen:
    sel.append({k: v for k, v in c.items() if k not in ('hlA', 'hlB', 'mcs_core')})
    sel[-1].update({
        'pic50_a': round(float(pic[c['a']]), 3),
        'pic50_b': round(float(pic[c['b']]), 3),
        'smiles_a': smi[c['a']], 'smiles_b': smi[c['b']],
    })
(ROOT / 'outputs/audit').mkdir(parents=True, exist_ok=True)
json.dump(sel, open(ROOT / 'outputs/audit/E4_cliff_pairs.json', 'w'), indent=2)


def render_pair(mA, mB, hlA, hlB, picA, picB):
    """Draw aligned pair. Higher-pIC50 diff atoms = green (gain), lower = red (loss)."""
    d = rdMolDraw2D.MolDraw2DCairo(1800, 780, 900, 780)
    opts = d.drawOptions()
    opts.legendFontSize  = 34          # larger pIC50 legends to survive journal scaling
    opts.padding         = 0.12
    opts.bondLineWidth   = 2.4
    opts.highlightRadius = 0.50

    # Colourblind-safe blue/orange (replaces red/green)
    blue   = (0.00, 0.467, 0.733)  # #0077BB — gain (higher pIC50)
    orange = (0.933, 0.467, 0.20)  # #EE7733 — loss (lower pIC50)

    # Higher pIC50 = gain (blue), lower = loss (orange)
    if picA >= picB:
        colA, colB = blue, orange   # mA is the more potent → its diff atoms are "gain context"
    else:
        colA, colB = orange, blue

    hlcolA = {i: colA for i in hlA}
    hlcolB = {i: colB for i in hlB}
    hlbondsA = [b.GetIdx() for b in mA.GetBonds()
                if b.GetBeginAtomIdx() in set(hlA) and b.GetEndAtomIdx() in set(hlA)]
    hlbondsB = [b.GetIdx() for b in mB.GetBonds()
                if b.GetBeginAtomIdx() in set(hlB) and b.GetEndAtomIdx() in set(hlB)]
    hlbondcolA = {i: colA for i in hlbondsA}
    hlbondcolB = {i: colB for i in hlbondsB}

    d.DrawMolecules([mA, mB],
                    highlightAtoms=[hlA, hlB],
                    highlightAtomColors=[hlcolA, hlcolB],
                    highlightBonds=[hlbondsA, hlbondsB],
                    highlightBondColors=[hlbondcolA, hlbondcolB],
                    legends=[f'pIC50 = {picA:.2f}', f'pIC50 = {picB:.2f}'])
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText()))


panel_imgs = []
panel_meta = []
for c in chosen:
    mA = Chem.MolFromSmiles(smi[c['a']])
    mB = Chem.MolFromSmiles(smi[c['b']])
    AllChem.Compute2DCoords(mA)
    mcs_core = c['mcs_core']
    if mcs_core is not None and mA.HasSubstructMatch(mcs_core) and mB.HasSubstructMatch(mcs_core):
        try:
            AllChem.GenerateDepictionMatching2DStructure(mB, mA, refPatt=mcs_core)
        except Exception:
            AllChem.Compute2DCoords(mB)
    else:
        AllChem.Compute2DCoords(mB)

    picA = pic[c['a']]
    picB = pic[c['b']]
    panel_imgs.append(render_pair(mA, mB, c['hlA'], c['hlB'], picA, picB))
    hubs = [f'Hub {h}' for h in [c['hub_a'], c['hub_b']] if h != 'none']
    hub_str = ' + '.join(hubs) if hubs else 'non-hub'
    panel_meta.append({
        'label': (f"ΔpIC50 = {abs(c['delta']):.2f}  Tan = {c['tan']:.3f}\n"
                  f"{c['mmp_type'].replace('_', ' ')} · {hub_str}"),
        'delta': abs(c['delta']), 'tan': c['tan'], 'mmp_type': c['mmp_type'],
    })

# ── Build figure: 10% top bar + 90% 2×2 grid ─────────────────────────────────
fig = plt.figure(figsize=(DOUBLE, 6.3), constrained_layout=True)
gs  = fig.add_gridspec(2, 1, height_ratios=[0.15, 0.85], hspace=0.04)
ax_bar = fig.add_subplot(gs[0])
gs_grid = gs[1].subgridspec(2, 2, hspace=0.10, wspace=0.04)

# MMP change-type bar (horizontal)
mmp_types   = ['single_atom_change', 'small_substituent', 'medium_substituent']
mmp_ns      = [int(mmp_type_counts.get(t, 0)) for t in mmp_types]
mmp_total   = sum(mmp_ns)
mmp_labels  = ['Single atom change', 'Small substituent', 'Medium substituent']
mmp_colors  = [SEM['single_atom'], SEM['small_subst'], SEM['medium']]

bars = ax_bar.barh([0, 1, 2], mmp_ns,
                    color=mmp_colors, height=0.62,
                    edgecolor='white', lw=0.4, alpha=0.90)
# count + percent of the 80 MMP-confirmed severe pairs (no redundant n= in y-label)
for i, n in enumerate(mmp_ns):
    ax_bar.text(n + 1.0, i, f'{n}  ({n/mmp_total*100:.0f}%)',
                va='center', ha='left', fontsize=6)

ax_bar.set_yticks([0, 1, 2])
ax_bar.set_yticklabels(mmp_labels, fontsize=6.5)
ax_bar.set_xlabel(f'MMP-confirmed severe cliff pairs (n={mmp_total})', fontsize=6.5)
ax_bar.set_xlim(0, max(mmp_ns) * 1.32)
ax_bar.invert_yaxis()
ax_bar.text(-0.085, 1.18, 'a', transform=ax_bar.transAxes,
            fontsize=8, fontweight='bold', va='bottom', ha='right')

# 2×2 grid
letters = ['i', 'ii', 'iii', 'iv']
for idx, (img, meta, letter) in enumerate(zip(panel_imgs, panel_meta, letters)):
    r, col = divmod(idx, 2)
    ax = fig.add_subplot(gs_grid[r, col])
    ax.imshow(img, aspect='auto')
    ax.axis('off')
    ax.text(-0.02, 1.03, f'({letter})', transform=ax.transAxes,
            fontsize=8, fontweight='bold', va='bottom', ha='right')
    ax.text(0.02, 0.99, meta['label'], transform=ax.transAxes,
            fontsize=6.5, va='top', ha='left', fontfamily='sans-serif')

# Legend for diff-atom colors
legend_handles = [
    mpatches.Patch(color='#0077BB', label='Higher-pIC50 molecule diff atoms (gain)'),
    mpatches.Patch(color='#EE7733', label='Lower-pIC50 molecule diff atoms (loss)'),
]
fig.legend(handles=legend_handles, fontsize=7, loc='lower center', ncol=2,
           framealpha=0.88, edgecolor='none', bbox_to_anchor=(0.5, -0.01))

# ── Save ──────────────────────────────────────────────────────────────────────
save_fig(fig, str(OUT / 'fig06b_cliff_pairs'))
plt.close(fig)

print("\n=== SELECTED PAIRS ===")
for s in sel:
    print(f"  [{s['rule']}] {s['a']} (hub {s['hub_a']}) <-> "
          f"{s['b']} (hub {s['hub_b']})  "
          f"d={s['delta']:.3f}  Tan={s['tan']:.3f}  {s['mmp_type']}")
print("DONE")
