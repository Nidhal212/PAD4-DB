"""
audit_export_package.py — Phase 7.4

Assembles the publication-ready supplementary resource package:
  activity_cliffs.csv, cliff_hubs.csv, scaffold_ruggedness.csv,
  mmp_transformations.csv, source_independence_scores.csv,
  standardized_structures.sdf

Output dir: outputs/supplementary_package/
"""
from pathlib import Path
import pandas as pd
from rdkit import Chem

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/supplementary_package'
OUT.mkdir(parents=True, exist_ok=True)

assets = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
ac = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')

print("=" * 64); print("PHASE 7.4 — SUPPLEMENTARY RESOURCE PACKAGE"); print("=" * 64)

# 1. activity cliffs (all tiers)
ac.to_csv(OUT / 'activity_cliffs.csv', index=False)

# 2. cliff hubs
HUBS = {'SMADULGDNOCLOP-GISFHXKWSA-N': 'A', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N': 'A',
        'UDCDEKJNAMHBFH-HSZRJFAPSA-N': 'B', 'DVCKJOQIVOGXEI-XMMPIXPASA-N': 'B'}
hub_df = assets[assets.inchi_key.isin(HUBS)][
    ['inchi_key', 'smiles_std', 'pIC50', 'murcko_smiles', 'mol_weight', 'source_list']].copy()
hub_df['hub_class'] = hub_df.inchi_key.map(HUBS)
hub_df.to_csv(OUT / 'cliff_hubs.csv', index=False)

# 3. scaffold ruggedness (from Phase 2)
rug = ROOT / 'outputs/tables/scaffold_ruggedness_table.csv'
if rug.exists():
    pd.read_csv(rug).to_csv(OUT / 'scaffold_ruggedness.csv', index=False)

# 4. mmp transformations (from Phase 3)
for src, dst in [('outputs/audit/top25_dangerous_transformations.csv', 'mmp_transformations.csv'),
                 ('outputs/mmp/mmp_pairs_cliff99.csv', 'mmp_pairs_all.csv')]:
    p = ROOT / src
    if p.exists():
        pd.read_csv(p).to_csv(OUT / dst, index=False)

# 5. source independence scores
assets[['inchi_key', 'source_list', 'n_sources', 'source_independence_score',
        'is_true_multi_source', 'multi_source', 'concordant']].to_csv(
    OUT / 'source_independence_scores.csv', index=False)

# 6. standardized structures SDF
w = Chem.SDWriter(str(OUT / 'standardized_structures.sdf'))
n = 0
for _, r in assets.iterrows():
    m = Chem.MolFromSmiles(r.smiles_std)
    if m is None:
        continue
    m.SetProp('_Name', r.inchi_key)
    m.SetProp('pIC50', f"{r.pIC50:.3f}")
    m.SetProp('murcko_scaffold', str(r.murcko_smiles))
    m.SetProp('source_independence_score', f"{r.source_independence_score}")
    if 'hub_class' in assets.columns:
        m.SetProp('hub_class', str(r.hub_class))
    w.write(m); n += 1
w.close()

print(f"  exported {n} structures to SDF")
for f in sorted(OUT.glob('*')):
    print(f"  {f.name:42s} {f.stat().st_size/1024:8.1f} KB")
print("DONE")
