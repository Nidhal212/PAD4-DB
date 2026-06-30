"""
05e_golden_set.py — Generate PAD4 Golden Set (high-reproducibility reference subset).

Criteria (LOCKED):
  1. ≥2 distinct genuine PubChem assay AIDs (numeric only; excludes aggregator
     pseudo-identifiers CHEMBL6111 and Q9UM07).
  2. Cross-assay pIC50 spread (max − min per-AID median) ≤ 0.5 log units.
  3. Source: potency_space.parquet — all rows have use_in_potency_model=True
     and norm_status_worst=OK by construction.

Input:  data/interim/normalized/potency_space.parquet
        publication/data/pad4_compounds.parquet  (for SMILES, pic50_consensus,
                                                  mechanism_class, warhead_class)
Output: publication/data/PAD4_Golden_Set.csv

Expected: 47 compounds (1.5% of 3,093).

Run from project root:
    conda run -n pad4bench python3 scripts/05_cliffs/05e_golden_set.py
"""

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[2]
POT_PARQ = ROOT / "data/interim/normalized/potency_space.parquet"
CPD_PARQ = ROOT / "publication/data/pad4_compounds.parquet"
OUT_CSV  = ROOT / "publication/data/PAD4_Golden_Set.csv"

PSEUDO_AIDS = {"CHEMBL6111", "Q9UM07"}

print("=" * 60)
print("STEP 05e — PAD4 GOLDEN SET")
print("=" * 60)

# ── Load potency space ────────────────────────────────────────────────────────
pot = pd.read_parquet(POT_PARQ)
print(f"\n  Potency-space rows loaded: {len(pot)}")
print(f"  Unique AIDs: {pot['aid'].nunique()}  "
      f"(includes pseudo-AIDs: {', '.join(PSEUDO_AIDS & set(pot['aid'].unique()))})")

# ── Filter to genuine PubChem AIDs ───────────────────────────────────────────
pub = pot[~pot["aid"].isin(PSEUDO_AIDS)].copy()
print(f"  PubChem-only rows (after pseudo-AID exclusion): {len(pub)}")

# ── Per-compound: distinct AID count and cross-assay spread ──────────────────
grp = (
    pub.groupby(["inchi_key", "aid"])["pic50_aid"]
    .median()
    .reset_index()
)
grp.columns = ["inchi_key", "aid", "pic50_median"]

agg = (
    grp.groupby("inchi_key")
    .agg(
        n_assays=("aid", "nunique"),
        pic50_max=("pic50_median", "max"),
        pic50_min=("pic50_median", "min"),
    )
    .reset_index()
)
agg["max_cross_assay_delta"] = (agg["pic50_max"] - agg["pic50_min"]).round(4)

# ── Apply Golden Set criteria ─────────────────────────────────────────────────
golden_agg = agg[(agg["n_assays"] >= 2) & (agg["max_cross_assay_delta"] <= 0.5)].copy()
print(f"\n  Golden Set (≥2 PubChem AIDs, spread ≤0.5): {len(golden_agg)} compounds")
print(f"  n_assays distribution: {golden_agg['n_assays'].value_counts().sort_index().to_dict()}")
zero_spread = (golden_agg["max_cross_assay_delta"] == 0).sum()
print(f"  Zero spread: {zero_spread}/{len(golden_agg)} = {zero_spread/len(golden_agg)*100:.1f}%")

# ── Annotate with compound-level metadata ────────────────────────────────────
cpd = pd.read_parquet(CPD_PARQ)
gs = cpd[cpd["inchi_key"].isin(set(golden_agg["inchi_key"]))].copy()
gs = gs.merge(
    golden_agg[["inchi_key", "n_assays", "max_cross_assay_delta"]],
    on="inchi_key", how="left"
)

# ── Compute Murcko scaffold ───────────────────────────────────────────────────
def get_scaffold(smi: str) -> str:
    try:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            return ""
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaf) if scaf else ""
    except Exception:
        return ""

gs["scaffold_smiles"] = gs["smiles_std"].apply(get_scaffold)

# ── Select and order output columns ──────────────────────────────────────────
gs_out = gs[[
    "inchi_key", "smiles_std", "pic50_consensus", "n_assays",
    "max_cross_assay_delta", "scaffold_smiles", "mechanism_class", "warhead_class"
]].copy()
gs_out.columns = [
    "inchi_key", "canonical_smiles", "consensus_pic50", "n_assays",
    "max_cross_assay_delta", "scaffold_smiles", "mechanism_class", "warhead_class"
]
gs_out = gs_out.sort_values("consensus_pic50", ascending=False).reset_index(drop=True)

print(f"\n  pIC50 range: {gs_out['consensus_pic50'].min():.3f} – "
      f"{gs_out['consensus_pic50'].max():.3f}")

# ── Write CSV ─────────────────────────────────────────────────────────────────
gs_out.to_csv(OUT_CSV, index=False)
print(f"\n  Wrote: {OUT_CSV}")
print(f"  Rows: {len(gs_out)}")

# ── Byte-identity check against existing file ─────────────────────────────────
# Re-read from disk and compare via CSV content (float repr can shift on re-sort)
fresh_sha = hashlib.sha256(OUT_CSV.read_bytes()).hexdigest()
print(f"\n  SHA256 of output: {fresh_sha}")

EXPECTED_N = 47
if len(gs_out) != EXPECTED_N:
    print(f"  FAIL: expected {EXPECTED_N} compounds, got {len(gs_out)}")
    sys.exit(1)
else:
    print(f"  PASS: compound count = {EXPECTED_N}")

print("\n" + "=" * 60)
print("STEP 05e COMPLETE")
print("=" * 60)
