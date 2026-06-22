#!/usr/bin/env python3
"""
PAD4-DB v2 — Step 04b: Add source independence scores to pad4_compounds
scripts/04_dedup/04b_add_independence_scores.py

Adds two columns to data/processed/pad4_compounds.parquet:
  source_independence_score  float  (0.3 / 0.5 / 0.6 / 0.7 / 1.0)
  is_true_multi_source       bool   (score >= 0.7)

Overwrites data/processed/pad4_compounds.parquet in-place.
"""

import sys
from pathlib import Path

import pandas as pd

COMP_PATH = Path("data/processed/pad4_compounds.parquet")

if not COMP_PATH.exists():
    sys.exit(f"ERROR: required input not found: {COMP_PATH}")

df = pd.read_parquet(COMP_PATH)
print(f"Loaded {COMP_PATH}: {len(df):,} rows, {df.shape[1]} cols")

# ── Source independence score mapping ───────────────────────────────────────
# source_list is "|".join(sorted(unique sources)); use frozenset for lookup.

_SCORE_MAP: dict[frozenset, float] = {
    frozenset({"bindingdb", "chembl", "pubchem_confirmatory"}): 0.3,
    frozenset({"bindingdb", "pubchem_confirmatory"}):           0.5,
    frozenset({"bindingdb", "chembl"}):                        0.6,
}


def _independence_score(source_list: str) -> float:
    sources = frozenset(source_list.split("|"))
    if len(sources) == 1:
        return 1.0
    return _SCORE_MAP.get(sources, 0.7)


df["source_independence_score"] = df["source_list"].map(_independence_score)
df["is_true_multi_source"]      = df["source_independence_score"] >= 0.7

# ── Overwrite ────────────────────────────────────────────────────────────────
df.to_parquet(COMP_PATH, index=False)
print(f"Overwritten → {COMP_PATH}  ({len(df):,} rows, {df.shape[1]} cols)")

# ── Report ───────────────────────────────────────────────────────────────────
print("\nsource_independence_score distribution:")
score_dist = df["source_independence_score"].value_counts().sort_index()
for score, count in score_dist.items():
    print(f"  {score:.1f}: {count:,} compounds")

n_true  = df["is_true_multi_source"].sum()
n_false = (~df["is_true_multi_source"]).sum()
print(f"\nis_true_multi_source: {n_true:,} True / {n_false:,} False")

# ── Source list breakdown for context ───────────────────────────────────────
print("\nsource_list breakdown (top 10):")
for sl, cnt in df["source_list"].value_counts().head(10).items():
    sc = _independence_score(sl)
    print(f"  {sl!r:55s}  score={sc:.1f}  n={cnt:,}")
