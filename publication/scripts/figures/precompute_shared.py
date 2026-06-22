"""
precompute_shared.py — Build shared_assets.parquet for all figure scripts.

Outputs:
  data/interim/shared_assets.parquet  — enriched compound table
  data/interim/scaffold_family_map.csv — scaffold → family mapping
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

# ── Canonical numbers ──────────────────────────────────────────────────────────
CANON = {
    'n_compounds': 3093,
    'pic50_min': 2.00,
    'pic50_max': 8.52,
    'n_pairs_sim_ge06': 358416,
    'n_severe': 94,
    'n_moderate': 193,
    'n_broad': 580,
    'n_in_severe': 99,
    'max_delta_severe': 3.045,
    'max_delta_all': 3.228,
    'n_scaffolds': 1244,
    'n_series': 375,
    'n_singletons': 869,
    'largest_series': 174,
    'n_in_series': 2224,
    'n_patent': 233,
    'n_patent_scaffolds': 103,
    'n_patent_cliffs': 1,
    'n_multi_06': 528,
    'n_multi_07': 361,
    'hub_a1_degree': 15,
    'hub_a2_degree': 12,
    'hub_b1_degree': 12,
    'hub_b2_degree': 11,
    'n_mmp_validated': 85,
    'n_single_rgroup': 49,
    'n_hts': 327336,
    'n_hts_confirmed': 1453,
    'n_unique_ik': 328976,
}

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
    'olive':      '#999933',
    'purple':     '#AA4499',
}

fam_to_color = {
    'azaindole-benzimidazole biaryl amide derivatives': COLORS['navy'],
    'indazole-N-alkylindole biaryl amide derivatives': COLORS['blue'],
    'indazole-azaindole biaryl amide derivatives': COLORS['blue'],
    'indole-benzimidazole biaryl amide derivatives': COLORS['blue'],
    'chalcone-oxindole derivatives': COLORS['orange'],
    'chalcone-bicyclic lactam derivatives': COLORS['orange'],
    'benzimidazolyl-dihydroisoquinolinone derivatives': COLORS['teal'],
    'bis-benzimidazolyl biaryl diamide derivatives': COLORS['cyan'],
    'Other': COLORS['grey'],
}

ROOT = Path('/home/nidhal/PAD4-db_V2')

print("=" * 60)
print("PRECOMPUTE SHARED ASSETS")
print("=" * 60)

# ── Load main compound table ───────────────────────────────────────────────────
print("\n[1] Loading pad4_compounds.parquet ...")
df = pd.read_parquet(ROOT / 'data/processed/pad4_compounds.parquet')
assert len(df) == CANON['n_compounds'], f"n_compounds mismatch: {len(df)} != {CANON['n_compounds']}"
print(f"    Loaded {len(df)} compounds ✓")

# ── Column aliases ─────────────────────────────────────────────────────────────
df['pIC50'] = df['pic50_consensus']
df['assay_mechanism'] = df['mechanism_class']

# ── Patent flag ────────────────────────────────────────────────────────────────
df['patent_flag'] = ~df['source_list'].str.contains('chembl|bindingdb', na=False)
n_patent = df['patent_flag'].sum()
assert n_patent == CANON['n_patent'], f"n_patent mismatch: {n_patent} != {CANON['n_patent']}"
print(f"    Patent flag: {n_patent} compounds ✓")

# ── Validate pIC50 range ───────────────────────────────────────────────────────
pic50_min = round(df['pIC50'].min(), 2)
pic50_max = round(df['pIC50'].max(), 2)
assert pic50_min == CANON['pic50_min'], f"pic50_min: {pic50_min}"
assert pic50_max == CANON['pic50_max'], f"pic50_max: {pic50_max}"
print(f"    pIC50 range: {pic50_min}–{pic50_max} ✓")

# ── Validate multi-source ──────────────────────────────────────────────────────
n_multi_06 = (df['source_independence_score'] >= 0.6).sum()
n_multi_07 = (df['source_independence_score'] >= 0.7).sum()
assert n_multi_06 == CANON['n_multi_06'], f"n_multi_06: {n_multi_06}"
assert n_multi_07 == CANON['n_multi_07'], f"n_multi_07: {n_multi_07}"
print(f"    Multi-source ≥0.6: {n_multi_06} ✓ | ≥0.7: {n_multi_07} ✓")

# ── Compute Murcko scaffold SMILES ─────────────────────────────────────────────
print("\n[2] Computing Murcko scaffold SMILES ...")

def get_murcko_smiles(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        core = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(core)
    except Exception:
        return None

df['murcko_smiles'] = df['smiles_std'].apply(get_murcko_smiles)
n_ok = df['murcko_smiles'].notna().sum()
print(f"    Computed Murcko SMILES for {n_ok}/{len(df)} compounds")

# ── Load scaffold summary ──────────────────────────────────────────────────────
print("\n[3] Loading scaffold summary ...")
scaffold_sum = pd.read_csv(ROOT / 'outputs/tables/05_scaffold_summary.csv')
assert len(scaffold_sum) == CANON['n_scaffolds'], f"n_scaffolds: {len(scaffold_sum)}"
scaffold_sum = scaffold_sum.reset_index().rename(columns={'index': 'scaffold_rank'})
scaffold_sum['scaffold_rank'] = scaffold_sum['scaffold_rank'] + 1  # 1-indexed

# Validate
n_series = (scaffold_sum['n_compounds'] >= 2).sum()
n_singletons = (scaffold_sum['n_compounds'] == 1).sum()
largest = scaffold_sum['n_compounds'].max()
print(f"    Scaffolds: {len(scaffold_sum)} | Series: {n_series} | Singletons: {n_singletons} | Largest: {largest}")
assert n_series == CANON['n_series'], f"n_series: {n_series} != {CANON['n_series']}"
assert n_singletons == CANON['n_singletons'], f"n_singletons: {n_singletons}"
assert largest == CANON['largest_series'], f"largest_series: {largest}"
print(f"    ✓")

# ── Load top-50 scaffold review ────────────────────────────────────────────────
print("\n[4] Loading scaffold_top50_review.csv ...")
top50 = pd.read_csv(ROOT / 'outputs/tables/scaffold_top50_review.csv')
print(f"    Top-50 columns: {list(top50.columns)}")

# Normalize scaffold_smiles for matching
def norm_smiles(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            return Chem.MolToSmiles(mol)
    except Exception:
        pass
    return smi

print("    Normalizing scaffold SMILES for matching ...")
scaffold_sum['scaffold_smiles_canon'] = scaffold_sum['scaffold_smiles'].apply(norm_smiles)
top50['scaffold_smiles_canon'] = top50['scaffold_smiles'].apply(norm_smiles)

# Merge scaffold_rank and family onto scaffold_sum
top50_merge = top50[['scaffold_smiles_canon', 'scaffold_rank', 'scaffold_family']].copy()
top50_merge.columns = ['scaffold_smiles_canon', 'top50_rank', 'scaffold_family']

# Use top50 ranks directly
scaffold_sum_merged = scaffold_sum.merge(
    top50_merge, on='scaffold_smiles_canon', how='left'
)
# Fill family for non-top-50
scaffold_sum_merged['scaffold_family'] = scaffold_sum_merged['scaffold_family'].fillna('Other')

# For top50 matches, use top50 rank; otherwise keep scaffold_rank
scaffold_sum_merged['scaffold_rank_final'] = scaffold_sum_merged['top50_rank'].fillna(
    scaffold_sum_merged['scaffold_rank']
).astype(int)

print(f"    Matched {scaffold_sum_merged['top50_rank'].notna().sum()} top-50 scaffolds")
print(f"    Family distribution:\n{scaffold_sum_merged['scaffold_family'].value_counts().head(10)}")

# ── Map family back to compounds via Murcko scaffold ─────────────────────────
print("\n[5] Joining scaffold family to compounds ...")

# Build compound → scaffold mapping
scaffold_to_family = dict(zip(
    scaffold_sum_merged['scaffold_smiles_canon'],
    scaffold_sum_merged['scaffold_family']
))
scaffold_to_rank = dict(zip(
    scaffold_sum_merged['scaffold_smiles_canon'],
    scaffold_sum_merged['scaffold_rank_final']
))
scaffold_to_size = dict(zip(
    scaffold_sum_merged['scaffold_smiles_canon'],
    scaffold_sum_merged['n_compounds']
))

def norm_smiles_safe(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            return Chem.MolToSmiles(mol)
    except Exception:
        pass
    return smi

df['murcko_smiles_canon'] = df['murcko_smiles'].apply(
    lambda x: norm_smiles_safe(x) if pd.notna(x) else None
)

df['scaffold_family'] = df['murcko_smiles_canon'].map(scaffold_to_family).fillna('Other')
df['scaffold_rank'] = df['murcko_smiles_canon'].map(scaffold_to_rank)
df['scaffold_series_size'] = df['murcko_smiles_canon'].map(scaffold_to_size)

# Map family to simplified group for visualization
FAM_GROUPS = {
    'azaindole-benzimidazole biaryl amide derivatives': 'azaindole-benzimidazole biaryl amide derivatives',
    'indazole-N-alkylindole biaryl amide derivatives': 'indazole-N-alkylindole biaryl amide derivatives',
    'indazole-azaindole biaryl amide derivatives': 'indazole-azaindole biaryl amide derivatives',
    'indole-benzimidazole biaryl amide derivatives': 'indole-benzimidazole biaryl amide derivatives',
    'chalcone-oxindole derivatives': 'chalcone-oxindole derivatives',
    'chalcone-bicyclic lactam derivatives': 'chalcone-bicyclic lactam derivatives',
    'benzimidazolyl-dihydroisoquinolinone derivatives': 'benzimidazolyl-dihydroisoquinolinone derivatives',
    'bis-benzimidazolyl biaryl diamide derivatives': 'bis-benzimidazolyl biaryl diamide derivatives',
}
df['scaffold_family_group'] = df['scaffold_family'].apply(
    lambda x: x if x in FAM_GROUPS else 'Other'
)

print(f"    Family group distribution:\n{df['scaffold_family_group'].value_counts()}")

# ── Validate n_in_series ──────────────────────────────────────────────────────
n_in_series = (df['scaffold_series_size'] >= 2).sum()
print(f"\n    n_in_series = {n_in_series} (canonical: {CANON['n_in_series']})")
# Allow slight discrepancy since spec notes 2221 vs 2224
if abs(n_in_series - CANON['n_in_series']) > 10:
    print(f"    WARNING: n_in_series differs by {abs(n_in_series - CANON['n_in_series'])}")
else:
    print(f"    ✓ within tolerance")

# ── Save scaffold family map ──────────────────────────────────────────────────
print("\n[6] Saving scaffold_family_map.csv ...")
fam_map = scaffold_sum_merged[['scaffold_smiles', 'scaffold_smiles_canon', 'scaffold_family',
                                 'scaffold_rank_final', 'n_compounds', 'mean_pic50']].copy()
fam_map.to_csv(ROOT / 'data/interim/scaffold_family_map.csv', index=False)
print(f"    Saved {len(fam_map)} rows")

# ── Save shared_assets.parquet ────────────────────────────────────────────────
print("\n[7] Saving shared_assets.parquet ...")
out_cols = list(df.columns)
df.to_parquet(ROOT / 'data/interim/shared_assets.parquet', index=False)
print(f"    Saved {len(df)} rows × {len(df.columns)} cols")

# ── Final summary ─────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"  n_compounds:    {len(df)} (canonical: {CANON['n_compounds']}) {'✓' if len(df)==CANON['n_compounds'] else '✗'}")
print(f"  n_patent:       {df['patent_flag'].sum()} (canonical: {CANON['n_patent']}) {'✓' if df['patent_flag'].sum()==CANON['n_patent'] else '✗'}")
print(f"  n_multi_06:     {n_multi_06} (canonical: {CANON['n_multi_06']}) {'✓' if n_multi_06==CANON['n_multi_06'] else '✗'}")
print(f"  n_scaffolds:    {len(scaffold_sum)} (canonical: {CANON['n_scaffolds']}) {'✓' if len(scaffold_sum)==CANON['n_scaffolds'] else '✗'}")
print(f"  n_series:       {n_series} (canonical: {CANON['n_series']}) {'✓' if n_series==CANON['n_series'] else '✗'}")
print(f"  n_singletons:   {n_singletons} (canonical: {CANON['n_singletons']}) {'✓' if n_singletons==CANON['n_singletons'] else '✗'}")
print(f"  largest_series: {largest} (canonical: {CANON['largest_series']}) {'✓' if largest==CANON['largest_series'] else '✗'}")
print(f"  n_in_series:    {n_in_series} (canonical: {CANON['n_in_series']}) {'~✓' if abs(n_in_series-CANON['n_in_series'])<=10 else '✗'}")
print("\nDone. shared_assets.parquet ready.")
