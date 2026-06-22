"""
supp_s07_physicochemical.py — S7: Physicochemical Property Landscape (8-panel grid)
2x4 grid with constrained_layout, panel labels a-h.
Outputs: publication/figures/supplementary/fig_s07_physicochemical.{png,pdf}
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, panel_label, save_fig, SEM, C, DOUBLE

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

set_style()
np.random.seed(42)

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'publication/figures/supplementary'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("SUPPLEMENTARY S7 — PHYSICOCHEMICAL PROPERTY LANDSCAPE")
print("=" * 60)

df = pd.read_parquet(ROOT / 'publication/data/pad4_compounds.parquet')
print(f"  Compounds: {len(df)}")

props = {'MolWt': [], 'MolLogP': [], 'TPSA': [], 'NumHAcceptors': [],
         'NumHDonors': [], 'NumRotatableBonds': [], 'FractionCSP3': [],
         'NumAromaticRings': []}
fail = 0
for smi in df['smiles_std']:
    m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    if m is None:
        fail += 1
        continue
    props['MolWt'].append(Descriptors.MolWt(m))
    props['MolLogP'].append(Crippen.MolLogP(m))
    props['TPSA'].append(rdMolDescriptors.CalcTPSA(m))
    props['NumHAcceptors'].append(rdMolDescriptors.CalcNumHBA(m))
    props['NumHDonors'].append(rdMolDescriptors.CalcNumHBD(m))
    props['NumRotatableBonds'].append(rdMolDescriptors.CalcNumRotatableBonds(m))
    props['FractionCSP3'].append(rdMolDescriptors.CalcFractionCSP3(m))
    props['NumAromaticRings'].append(rdMolDescriptors.CalcNumAromaticRings(m))

print(f"  Parse failures: {fail}")
P = {k: np.array(v) for k, v in props.items()}
n = len(P['MolWt'])

panels = [
    ('MolWt',            'Molecular weight (Da)', 'a', False),
    ('MolLogP',          'cLogP (Crippen)',        'b', False),
    ('TPSA',             'TPSA (Å²)',              'c', False),
    ('NumHAcceptors',    'H-bond acceptors',       'd', True),
    ('NumHDonors',       'H-bond donors',          'e', True),
    ('NumRotatableBonds','Rotatable bonds',         'f', True),
    ('FractionCSP3',     'Fraction Csp3',          'g', False),
    ('NumAromaticRings', 'Aromatic rings',          'h', True),
]

fig, axes = plt.subplots(2, 4, figsize=(DOUBLE, 3.4), constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=0.08, h_pad=0.06)
axes = axes.ravel()

for ax, (key, label, letter, integer) in zip(axes, panels):
    v = P[key]
    med = np.median(v)
    if integer:
        lo, hi = int(v.min()), int(v.max())
        bins = np.arange(lo - 0.5, hi + 1.5, 1)
    else:
        bins = 40
    ax.hist(v, bins=bins, color=SEM['published'], edgecolor='white', linewidth=0.2, alpha=0.8)
    ax.axvline(med, color=SEM['cliff'], linewidth=1.0, linestyle='--')

    # Place annotation in the empty half of the plot (avoid bars)
    # Determine which side has more space
    v_min, v_max = v.min(), v.max()
    midpoint = (v_min + v_max) / 2
    if med > midpoint:
        # More space on the left
        ann_x = 0.04
        ann_ha = 'left'
    else:
        ann_x = 0.96
        ann_ha = 'right'
    lbl_str = f'median {med:.2f}' if not integer else f'median {med:.0f}'
    ax.text(ann_x, 0.97, lbl_str,
            transform=ax.transAxes, ha=ann_ha, va='top',
            fontsize=5.5, color=SEM['cliff'])
    ax.set_xlabel(label)
    ax.set_ylabel('Compounds')
    if key == 'FractionCSP3':
        ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.tick_params(axis='x', labelsize=5.5)
    panel_label(ax, letter, x=-0.18, y=1.06)

save_fig(fig, str(OUT / 'fig_s07_physicochemical'))
plt.close(fig)

# Compute Lipinski Ro5 and Veber compliance
from rdkit.Chem import Descriptors as _Desc
ro5_violations = []
veber_pass = []
for smi in df['smiles_std']:
    m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    if m is None:
        continue
    mw   = _Desc.MolWt(m)
    logp = Crippen.MolLogP(m)
    hbd  = rdMolDescriptors.CalcNumHBD(m)
    hba  = rdMolDescriptors.CalcNumHBA(m)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(m)
    tpsa = rdMolDescriptors.CalcTPSA(m)
    ro5_violations.append(sum([mw > 500, logp > 5, hbd > 5, hba > 10]))
    veber_pass.append(rotb <= 10 and tpsa <= 140)

ro5_arr = np.array(ro5_violations)
n_strict  = int((ro5_arr == 0).sum())
n_classic = int((ro5_arr <= 1).sum())
n_veber   = int(sum(veber_pass))
total     = len(ro5_arr)

print(f"  Lipinski strict (0 violations): {n_strict}/{total} = {n_strict/total*100:.1f}%")
print(f"  Lipinski classic (<=1 violation): {n_classic}/{total} = {n_classic/total*100:.1f}%")
print(f"  Veber compliant: {n_veber}/{total} = {n_veber/total*100:.1f}%")

jpath = ROOT / 'outputs/audit/E2_physchem.json'
jpath.parent.mkdir(parents=True, exist_ok=True)
meta = {
    "n_compounds": total,
    "parse_failures": fail,
    "descriptors": {
        k: {
            "median": float(round(np.median(P[k]), 3)),
            "mean":   float(round(np.mean(P[k]), 3)),
            "p5":     float(round(np.percentile(P[k], 5), 3)),
            "p95":    float(round(np.percentile(P[k], 95), 3)),
        }
        for k in P
    },
    "lipinski_ro5": {
        "n_strict": n_strict, "n_classic": n_classic,
        "n_total": total,
        "fraction_strict": round(n_strict / total, 4),
        "fraction": round(n_strict / total, 4),
        "fraction_classic": round(n_classic / total, 4),
    },
    "veber": {"n_compliant": n_veber, "n_total": total,
              "fraction": round(n_veber / total, 4)},
}
with open(jpath, 'w') as f:
    json.dump(meta, f, indent=2)
print(f"  Wrote: {jpath}")
print("DONE")
