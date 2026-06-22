#!/usr/bin/env python3
"""
PAD4-DB v2 — Step 03b: Log-space aggregation validation (IC50 rows)
scripts/03_aggregate/03b_logspace_qc.py

Independently recomputes mean(log10(value_nM)) per group from raw
normalized_activities.parquet and compares to potency_space.parquet.
Also reports pIC50 distribution and outlier fractions.

Output: outputs/tables/03b_logspace_qc.json
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

POT_PATH  = Path("data/interim/normalized/potency_space.parquet")
NORM_PATH = Path("data/interim/normalized/normalized_activities.parquet")
OUT_PATH  = Path("outputs/tables/03b_logspace_qc.json")
GROUP_COLS = ["inchi_key", "source", "aid", "endpoint_type"]

for p in [POT_PATH, NORM_PATH]:
    if not p.exists():
        sys.exit(f"ERROR: required input not found: {p}")

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

# ── Load IC50 rows from potency_space ──────────────────────────────────────
pot = pd.read_parquet(POT_PATH)
ic50_pot = pot[pot["endpoint_type"] == "IC50"].copy()
n_ic50_groups = len(ic50_pot)
print(f"IC50 groups in potency_space: {n_ic50_groups:,}")

# ── Independent recompute from normalized_activities ───────────────────────
norm = pd.read_parquet(NORM_PATH)
ic50_norm = norm[
    (norm["endpoint_type"] == "IC50") &
    norm["value_nM"].notna() &
    (norm["value_nM"] > 0)
].copy()
ic50_norm["log_nM"] = np.log10(ic50_norm["value_nM"])

recomputed = (
    ic50_norm
    .groupby(GROUP_COLS, sort=False)["log_nM"]
    .mean()
    .reset_index()
    .rename(columns={"log_nM": "log_value_mean_check"})
)
print(f"IC50 groups recomputed from normalized: {len(recomputed):,}")

# ── Join and compare ───────────────────────────────────────────────────────
merged = ic50_pot[GROUP_COLS + ["log_value_mean", "pic50_aid"]].merge(
    recomputed, on=GROUP_COLS, how="left"
)

# Rows in potency_space with a log_value_mean but no recomputed value
# (can happen if all replicates had value_nM <= 0 or null — rare)
has_both = merged["log_value_mean"].notna() & merged["log_value_mean_check"].notna()
diff = (merged.loc[has_both, "log_value_mean"] - merged.loc[has_both, "log_value_mean_check"]).abs()

max_diff  = float(diff.max()) if len(diff) > 0 else 0.0
mean_diff = float(diff.mean()) if len(diff) > 0 else 0.0

print(f"\nLog-space validation:")
print(f"  Compared groups:           {has_both.sum():,}")
print(f"  Max absolute difference:   {max_diff:.6f}")
print(f"  Mean absolute difference:  {mean_diff:.6f}")

if max_diff > 0.001:
    bad = merged[has_both & (diff > 0.001)]
    print(f"\n  ERROR: {len(bad)} groups exceed tolerance of 0.001:")
    print(bad[GROUP_COLS + ["log_value_mean", "log_value_mean_check"]].head(10).to_string())
else:
    print("  PASS — all diffs within tolerance")

# ── pIC50 distribution stats ───────────────────────────────────────────────
pic50_valid = ic50_pot["pic50_aid"].dropna()
percentiles = np.percentile(pic50_valid, [5, 25, 50, 75, 95])
dist = {
    "mean": float(pic50_valid.mean()),
    "std":  float(pic50_valid.std()),
    "min":  float(pic50_valid.min()),
    "max":  float(pic50_valid.max()),
    "p5":   float(percentiles[0]),
    "p25":  float(percentiles[1]),
    "p50":  float(percentiles[2]),
    "p75":  float(percentiles[3]),
    "p95":  float(percentiles[4]),
}
print(f"\npIC50 distribution (IC50 groups, n={len(pic50_valid):,}):")
for k, v in dist.items():
    print(f"  {k}: {v:.3f}")

# ── Outlier fractions ──────────────────────────────────────────────────────
n_below_3  = int((pic50_valid < 3).sum())
n_above_12 = int((pic50_valid > 12).sum())
frac_below = n_below_3  / len(pic50_valid) if len(pic50_valid) > 0 else 0.0
frac_above = n_above_12 / len(pic50_valid) if len(pic50_valid) > 0 else 0.0

print(f"\nOutliers:")
print(f"  pIC50 < 3  (>1 mM):  {n_below_3} rows  ({frac_below:.4%})")
print(f"  pIC50 > 12 (<1 pM):  {n_above_12} rows  ({frac_above:.4%})")

if n_below_3 > 0:
    print("  Rows below pIC50=3:")
    print(ic50_pot[ic50_pot["pic50_aid"] < 3][
        GROUP_COLS + ["pic50_aid", "value_nM_min", "value_nM_max", "norm_status_worst"]
    ].to_string())

if n_above_12 > 0:
    print("  Rows above pIC50=12:")
    print(ic50_pot[ic50_pot["pic50_aid"] > 12][
        GROUP_COLS + ["pic50_aid", "value_nM_min", "value_nM_max", "norm_status_worst"]
    ].to_string())

# ── Write JSON ─────────────────────────────────────────────────────────────
result = {
    "n_ic50_groups": n_ic50_groups,
    "n_compared": int(has_both.sum()),
    "max_abs_diff_log_mean": max_diff,
    "mean_abs_diff_log_mean": mean_diff,
    "pic50_distribution": dist,
    "n_outlier_below_3": n_below_3,
    "n_outlier_above_12": n_above_12,
    "outlier_fraction_below_3": frac_below,
    "outlier_fraction_above_12": frac_above,
}
OUT_PATH.write_text(json.dumps(result, indent=2))
print(f"\nWritten → {OUT_PATH}")
