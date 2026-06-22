"""
supp_hub_properties.py — Supplementary Analysis: Hub Compound Physicochemical Properties

Computes RDKit descriptors for:
  - 4 Hub compounds (Class A + Class B)
  - 95 Non-hub severe cliff compounds

Performs Wilcoxon rank-sum test (Mann-Whitney U) for each property.

Outputs:
  outputs/tables/supp_hub_properties.csv  — Property × Hub vs Non-hub table
  (Markdown table printed to stdout for direct copy into manuscript)
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path('/home/nidhal/PAD4-db_V2')

HUB_IKS = {
    'SMADULGDNOCLOP-GISFHXKWSA-N': 'A',
    'RAVBZQAQTVGKIV-XBPDSQQVSA-N': 'A',
    'UDCDEKJNAMHBFH-HSZRJFAPSA-N': 'B',
    'DVCKJOQIVOGXEI-XMMPIXPASA-N': 'B',
}

print("=" * 60)
print("SUPP: Hub Compound Physicochemical Characterization")
print("=" * 60)

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("WARNING: RDKit not available — using pre-computed descriptors from parquet")

# ── Load data ────────────────────────────────────────────────────────────────
pad   = pd.read_parquet(ROOT / 'publication/data/pad4_compounds.parquet')
ac    = pd.read_parquet(ROOT / 'publication/data/activity_cliffs.parquet')
sev   = ac[ac['cliff_tier'] == 'severe']
sev_iks = set(sev['inchi_key_a'].tolist() + sev['inchi_key_b'].tolist())

hub_iks    = set(HUB_IKS.keys())
nonhub_iks = sev_iks - hub_iks

pad_hub    = pad[pad['inchi_key'].isin(hub_iks)].copy()
pad_nonhub = pad[pad['inchi_key'].isin(nonhub_iks)].copy()

print(f"  Hub compounds:        {len(pad_hub)}")
print(f"  Non-hub cliff cmpds:  {len(pad_nonhub)}")

# ── Compute RDKit descriptors ────────────────────────────────────────────────
DESCRIPTOR_NAMES = [
    ('MW',      'Molecular Weight'),
    ('cLogP',   'Crippen cLogP'),
    ('TPSA',    'TPSA (Å²)'),
    ('HBD',     'H-bond donors'),
    ('HBA',     'H-bond acceptors'),
    ('RotBonds','Rotatable bonds'),
    ('ArRings', 'Aromatic rings'),
    ('Fsp3',    'Fraction Csp3'),
    ('pIC50',   'pIC50'),
]

def compute_descriptors(df):
    rows = []
    for _, row in df.iterrows():
        d = {'inchi_key': row['inchi_key'], 'pIC50': row.get('pic50_consensus', np.nan)}
        if HAS_RDKIT:
            mol = Chem.MolFromSmiles(row['smiles_std'])
            if mol:
                d['MW']       = Descriptors.MolWt(mol)
                d['cLogP']    = Descriptors.MolLogP(mol)
                d['TPSA']     = Descriptors.TPSA(mol)
                d['HBD']      = rdMolDescriptors.CalcNumHBD(mol)
                d['HBA']      = rdMolDescriptors.CalcNumHBA(mol)
                d['RotBonds'] = rdMolDescriptors.CalcNumRotatableBonds(mol)
                d['ArRings']  = rdMolDescriptors.CalcNumAromaticRings(mol)
                d['Fsp3']     = rdMolDescriptors.CalcFractionCSP3(mol)
            else:
                for k in ('MW','cLogP','TPSA','HBD','HBA','RotBonds','ArRings','Fsp3'):
                    d[k] = np.nan
        else:
            # Fallback to pre-computed columns if available
            d['MW']    = row.get('mol_weight', np.nan)
            for k in ('cLogP','TPSA','HBD','HBA','RotBonds','ArRings','Fsp3'):
                d[k] = np.nan
        rows.append(d)
    return pd.DataFrame(rows)

print("\n  Computing hub descriptors ...")
desc_hub    = compute_descriptors(pad_hub)
print("  Computing non-hub cliff descriptors ...")
desc_nonhub = compute_descriptors(pad_nonhub)

# ── Statistical comparison ────────────────────────────────────────────────────
prop_keys = ['MW','cLogP','TPSA','HBD','HBA','RotBonds','ArRings','Fsp3','pIC50']
prop_names = dict(DESCRIPTOR_NAMES)

results = []
for pk in prop_keys:
    h_vals  = desc_hub[pk].dropna().values
    nh_vals = desc_nonhub[pk].dropna().values
    if len(h_vals) == 0 or len(nh_vals) == 0:
        results.append({'Property': prop_names[pk],
                        'Hub Mean (n=4)': 'N/A', 'Non-hub Mean (n=95)': 'N/A',
                        'p-value': 'N/A'})
        continue
    _, p = stats.mannwhitneyu(h_vals, nh_vals, alternative='two-sided')
    p_str = '< 0.001' if p < 0.001 else f'{p:.3f}'
    results.append({
        'Property':             prop_names[pk],
        'Hub Mean (n=4)':       f'{np.mean(h_vals):.2f}',
        'Non-hub Mean (n=95)':  f'{np.mean(nh_vals):.2f}',
        'p-value':              p_str,
    })

result_df = pd.DataFrame(results)
out_dir = ROOT / 'outputs/tables'
out_dir.mkdir(parents=True, exist_ok=True)
result_df.to_csv(out_dir / 'supp_hub_properties.csv', index=False)

# ── Print Markdown table ──────────────────────────────────────────────────────
print("\n=== Supplementary Table SX: Hub vs Non-hub Physicochemical Properties ===\n")
print("| Property | Hub Average (n=4) | Non-hub Average (n=95) | p-value (MWU) |")
print("|----------|-------------------|------------------------|---------------|")
for _, r in result_df.iterrows():
    print(f"| {r['Property']} | {r['Hub Mean (n=4)']} | {r['Non-hub Mean (n=95)']} | {r['p-value']} |")
print()
print("*Wilcoxon rank-sum (Mann-Whitney U) test, two-sided.*")
print(f"\nOutput: {out_dir / 'supp_hub_properties.csv'}")
print("DONE")
