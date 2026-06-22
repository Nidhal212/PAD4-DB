#!/usr/bin/env python3
"""
PAD4-DB v2 — Step 03a: Hard split into potency_space vs hts_space
scripts/03_aggregate/03a_split_spaces.py

Splits replicate_aggregated.parquet by use_in_potency_model flag.
  potency_space = use_in_potency_model == True
  hts_space     = use_in_potency_model == False

Output:
  data/interim/normalized/potency_space.parquet
  data/interim/normalized/hts_space.parquet
"""

import sys
from pathlib import Path

import pandas as pd

IN_PATH      = Path("data/interim/normalized/replicate_aggregated.parquet")
POT_PATH     = Path("data/interim/normalized/potency_space.parquet")
HTS_PATH     = Path("data/interim/normalized/hts_space.parquet")
TOTAL_EXPECTED = 339_687

for p in [IN_PATH]:
    if not p.exists():
        sys.exit(f"ERROR: required input not found: {p}")

df = pd.read_parquet(IN_PATH)

potency = df[df["use_in_potency_model"] == True].copy()
hts     = df[df["use_in_potency_model"] == False].copy()

total = len(potency) + len(hts)
if total != TOTAL_EXPECTED:
    sys.exit(
        f"ASSERTION FAILED: potency({len(potency)}) + hts({len(hts)}) = {total} "
        f"!= expected {TOTAL_EXPECTED}"
    )

potency.to_parquet(POT_PATH, index=False)
hts.to_parquet(HTS_PATH, index=False)

print(f"potency_space rows:              {len(potency):,}")
print(f"hts_space rows:                  {len(hts):,}")
print(f"sum:                             {total:,}  ✓ matches {TOTAL_EXPECTED:,}")
print(f"unique InChIKeys in potency_space: {potency['inchi_key'].nunique():,}")
print(f"unique InChIKeys in hts_space:     {hts['inchi_key'].nunique():,}")
print(f"Written → {POT_PATH}")
print(f"Written → {HTS_PATH}")
