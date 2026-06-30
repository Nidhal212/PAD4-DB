#!/usr/bin/env python3
"""
PAD4-DB v2 — Add derived columns to master parquet
Columns added: is_covalent, warhead_class, mechanism_class, fragment_flag
All values verified by diagnostic_mechanism_v2.py (9/9 PASS)
Run from project root with pad4bench activated.
"""

import pandas as pd
from rdkit import Chem
import shutil, datetime, os

# ── Backup first ──────────────────────────────────────────────────────────────
src = 'data/processed/pad4_compounds.parquet'
bak = f'data/processed/pad4_compounds_pre_columns_{datetime.date.today()}.parquet'
shutil.copy2(src, bak)
print(f"Backup written: {bak}")

# ── Load ──────────────────────────────────────────────────────────────────────
compounds = pd.read_parquet(src)
print(f"Loaded: {len(compounds):,} compounds, {len(compounds.columns)} columns")
assert len(compounds) == 3093, f"Row count wrong: {len(compounds)}"

# ── SMARTS panel v2 (all validated) ──────────────────────────────────────────
WARHEADS = {
    'chloroacetamidine':  'NC(=N)CCl',
    'fluoroacetamidine':  'NC(=N)CF',
    'haloacetyl':         '[F,Cl,Br]CC(=O)N',
    'vinyl_sulfone':      'C=CS(=O)(=O)',
    'amidine_any_halide': 'NC(=N)C[F,Cl,Br]',
    'alpha_bromoketone':  'BrCC(=O)',
    'enaminone':          'NC=CC(=O)',
}
compiled = {k: Chem.MolFromSmarts(v) for k, v in WARHEADS.items()}
assert all(v is not None for v in compiled.values()), "SMARTS compile failure"
print(f"SMARTS panel: {len(compiled)} warheads compiled")

def get_warhead(smi):
    if not isinstance(smi, str):
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    for name, pat in compiled.items():
        if mol.HasSubstructMatch(pat):
            return name
    return None

def classify_mechanism(x):
    x = str(x)
    if 'covalent_irreversible' in x:
        return 'covalent'
    elif 'rfms_enzymatic' in x:
        return 'enzymatic_confirmed'
    elif 'baee_colorimetric' in x:
        return 'enzymatic'
    elif 'fp_binding' in x:
        return 'fp_ic50'
    return 'unknown'

# ── Compute columns ───────────────────────────────────────────────────────────
print("Computing warhead_class...")
compounds['warhead_class']   = compounds['smiles_std'].apply(get_warhead)

print("Computing is_covalent...")
compounds['is_covalent']     = compounds['warhead_class'].notna()

print("Computing mechanism_class...")
compounds['mechanism_class'] = compounds['assay_mechanism_classes'].apply(
    classify_mechanism)

print("Computing fragment_flag...")
compounds['fragment_flag']   = (
    (compounds['mol_weight'] < 200) &
    (compounds['pic50_consensus'] < 4.0)
)

# ── Verification ──────────────────────────────────────────────────────────────
print("\n=== VERIFICATION ===")
checks = {
    'Row count = 3,093':          len(compounds) == 3093,
    'is_covalent True = 107':     compounds['is_covalent'].sum() == 107,
    'fragment_flag True = 5':     compounds['fragment_flag'].sum() == 5,
    'mechanism_class unknown = 0':(compounds['mechanism_class']=='unknown').sum()==0,
    'mechanism_class covalent=21':(compounds['mechanism_class']=='covalent').sum()==21,
    'No null warhead for covalent':
        compounds[compounds['is_covalent']]['warhead_class'].notna().all(),
    'Columns added = 4':
        all(c in compounds.columns for c in
            ['is_covalent','warhead_class','mechanism_class','fragment_flag']),
}

all_pass = True
for desc, result in checks.items():
    status = "✅ PASS" if result else "❌ FAIL"
    if not result:
        all_pass = False
    print(f"  {status}  {desc}")

if not all_pass:
    print("\nFAILURES DETECTED — NOT writing parquet. Fix issues first.")
    raise SystemExit(1)

# ── Write ─────────────────────────────────────────────────────────────────────
compounds.to_parquet(src, index=False)
print(f"\n✅ Written: {src}")
print(f"   Rows: {len(compounds):,}")
print(f"   Columns: {len(compounds.columns)} (was 21, now 25)")

# ── Print final column list ───────────────────────────────────────────────────
print(f"\nFinal columns:")
for i, col in enumerate(compounds.columns):
    tag = " ← NEW" if col in (
        'is_covalent','warhead_class','mechanism_class','fragment_flag') else ""
    print(f"  {i+1:>2}. {col}{tag}")

# ── Distribution summary ──────────────────────────────────────────────────────
print(f"\nwarhead_class distribution:")
print(compounds['warhead_class'].value_counts(dropna=False).to_string())
print(f"\nmechanism_class distribution:")
print(compounds['mechanism_class'].value_counts().to_string())
print(f"\nfragment_flag: {compounds['fragment_flag'].sum()} True, "
      f"{(~compounds['fragment_flag']).sum()} False")