#!/usr/bin/env python3
"""
PAD4-DB v2 — Step 03c: SMILES integrity check
scripts/03_aggregate/03c_smiles_integrity.py

Validates three invariants on potency_space.parquet:
  1. Many-InChIKey-to-one-SMILES:  smiles_std → n_distinct_inchi_key ≤ 1
  2. One-InChIKey-to-many-SMILES:  inchi_key → n_distinct_smiles_std = 1 (hard requirement)
  3. Null SMILES count = 0

Output: outputs/tables/03c_smiles_integrity.json
"""

import json
import sys
from pathlib import Path

import pandas as pd

POT_PATH = Path("data/interim/normalized/potency_space.parquet")
OUT_PATH = Path("outputs/tables/03c_smiles_integrity.json")

if not POT_PATH.exists():
    sys.exit(f"ERROR: required input not found: {POT_PATH}")

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

pot = pd.read_parquet(POT_PATH)

# ── Check 3: null SMILES ───────────────────────────────────────────────────
n_null = int(pot["smiles_std"].isna().sum())
print(f"Check 3 — Null smiles_std: {n_null}")
if n_null > 0:
    print("  WARNING: null SMILES present:")
    print(pot[pot["smiles_std"].isna()][["inchi_key", "source", "aid"]].head(10).to_string())

# ── Check 1: SMILES → multiple InChIKeys ──────────────────────────────────
# (Indicates canonicalization failure: two distinct structures share a SMILES)
smiles_to_ik = (
    pot[pot["smiles_std"].notna()]
    .groupby("smiles_std")["inchi_key"]
    .nunique()
)
bad_smiles = smiles_to_ik[smiles_to_ik > 1].sort_values(ascending=False)
n_smiles_conflict = int(len(bad_smiles))

print(f"\nCheck 1 — SMILES mapping to >1 InChIKey: {n_smiles_conflict}")
if n_smiles_conflict > 0:
    print("  WARNING: canonicalization errors detected — SMILES maps to multiple InChIKeys:")
    top10 = bad_smiles.head(10).reset_index()
    top10.columns = ["smiles_std", "n_inchi_keys"]
    # Show which InChIKeys collide for each conflicted SMILES
    for _, row in top10.iterrows():
        iks = pot[pot["smiles_std"] == row["smiles_std"]]["inchi_key"].unique()
        print(f"    SMILES: {row['smiles_std'][:60]}...")
        print(f"    InChIKeys ({row['n_inchi_keys']}): {list(iks[:5])}")
else:
    print("  PASS — every SMILES maps to exactly one InChIKey")

top10_smiles_conflicts = (
    bad_smiles.head(10).reset_index()
    .rename(columns={"inchi_key": "n_inchi_keys"})
    .to_dict(orient="records")
) if n_smiles_conflict > 0 else []

# ── Check 2: InChIKey → multiple SMILES ───────────────────────────────────
# (Must be 0 — would indicate a bug in Step 01 standardization)
ik_to_smiles = (
    pot[pot["smiles_std"].notna()]
    .groupby("inchi_key")["smiles_std"]
    .nunique()
)
bad_ik = ik_to_smiles[ik_to_smiles > 1].sort_values(ascending=False)
n_ik_conflict = int(len(bad_ik))

print(f"\nCheck 2 — InChIKey mapping to >1 SMILES: {n_ik_conflict}")
if n_ik_conflict > 0:
    print("  ERROR: Step 01 standardization bug — same InChIKey has multiple SMILES:")
    top10_ik = bad_ik.head(10).reset_index()
    top10_ik.columns = ["inchi_key", "n_smiles"]
    for _, row in top10_ik.iterrows():
        smiles_list = pot[pot["inchi_key"] == row["inchi_key"]]["smiles_std"].unique()
        print(f"    InChIKey: {row['inchi_key']}")
        print(f"    SMILES ({row['n_smiles']}): {list(smiles_list[:3])}")
else:
    print("  PASS — every InChIKey maps to exactly one SMILES")

top10_ik_conflicts = (
    bad_ik.head(10).reset_index()
    .rename(columns={"smiles_std": "n_smiles"})
    .to_dict(orient="records")
) if n_ik_conflict > 0 else []

# ── Overall verdict ────────────────────────────────────────────────────────
print("\n=== SMILES INTEGRITY VERDICT ===")
if n_null == 0 and n_ik_conflict == 0:
    if n_smiles_conflict == 0:
        print("ALL CHECKS PASS")
    else:
        print(f"WARN: {n_smiles_conflict} SMILES map to >1 InChIKey (canonicalization note)")
        print("      These are likely salt/tautomer variants with identical canonical SMILES.")
        print("      Review top10 above before proceeding.")
else:
    if n_null > 0:
        print(f"FAIL: {n_null} null SMILES in potency_space")
    if n_ik_conflict > 0:
        print(f"FAIL: {n_ik_conflict} InChIKeys map to >1 SMILES — Step 01 bug")

# ── Write JSON ─────────────────────────────────────────────────────────────
result = {
    "smiles_null_count": n_null,
    "smiles_to_multiple_inchikeys": n_smiles_conflict,
    "inchikey_to_multiple_smiles": n_ik_conflict,
    "top10_smiles_conflicts": top10_smiles_conflicts,
    "top10_inchikey_conflicts": top10_ik_conflicts,
}
OUT_PATH.write_text(json.dumps(result, indent=2))
print(f"\nWritten → {OUT_PATH}")
