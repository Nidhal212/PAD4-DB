"""
00_validate_canonical_numbers.py — Pre-flight guard for all figure/table scripts.

Run this FIRST before regenerating any figure or table.
If any assertion fails, STOP and investigate the discrepancy.
"""
import sys
from pathlib import Path
import pandas as pd

ROOT = Path('/home/nidhal/PAD4-db_V2')

CANON = {
    'n_compounds':    3093,
    'n_severe_cliffs': 94,
    'n_hubs':          4,
    'n_scaffolds':    1244,
}

HUB_IKS = {
    'A1': 'SMADULGDNOCLOP-GISFHXKWSA-N',
    'A2': 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
    'B1': 'UDCDEKJNAMHBFH-HSZRJFAPSA-N',
    'B2': 'DVCKJOQIVOGXEI-XMMPIXPASA-N',
}

errors = []

def check(condition, msg):
    if not condition:
        errors.append(f"  FAIL: {msg}")
        print(f"  [FAIL] {msg}")
    else:
        print(f"  [OK]   {msg}")

print("=" * 60)
print("PAD4-DB v2 — CANONICAL NUMBER VALIDATION")
print("=" * 60)

# ── Compound dataset ──────────────────────────────────────────────────────────
print("\n[1] pad4_compounds.parquet")
try:
    pad = pd.read_parquet(ROOT / 'publication/data/pad4_compounds.parquet')
    check(len(pad) == CANON['n_compounds'],
          f"n_compounds == {CANON['n_compounds']}  (got {len(pad)})")

    # Hub compounds present
    for lbl, ik in HUB_IKS.items():
        check(ik in pad['inchi_key'].values,
              f"Hub {lbl} ({ik[:14]}...) present in pad4_compounds")

    # hub_class column present and 4 hubs labelled
    if 'hub_class' in pad.columns:
        n_hub_a = (pad['hub_class'] == 'A').sum()
        n_hub_b = (pad['hub_class'] == 'B').sum()
        check(n_hub_a == 2, f"hub_class A == 2  (got {n_hub_a})")
        check(n_hub_b == 2, f"hub_class B == 2  (got {n_hub_b})")
        check(n_hub_a + n_hub_b == CANON['n_hubs'],
              f"n_hubs == {CANON['n_hubs']}  (got {n_hub_a + n_hub_b})")
    else:
        errors.append("  FAIL: hub_class column missing from pad4_compounds.parquet")
        print("  [FAIL] hub_class column missing from pad4_compounds.parquet")

    # Scaffold count via shared_assets
    print("\n[2] shared_assets.parquet (scaffold count)")
    assets = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
    check(len(assets) == CANON['n_compounds'],
          f"shared_assets n_compounds == {CANON['n_compounds']}  (got {len(assets)})")

except FileNotFoundError as e:
    errors.append(f"  FAIL: {e}")
    print(f"  [FAIL] {e}")

# ── Scaffold summary ──────────────────────────────────────────────────────────
print("\n[3] 05_scaffold_summary.csv")
try:
    scaff = pd.read_csv(ROOT / 'outputs/tables/05_scaffold_summary.csv')
    check(len(scaff) == CANON['n_scaffolds'],
          f"n_scaffolds == {CANON['n_scaffolds']}  (got {len(scaff)})")
    top1 = scaff.sort_values('n_compounds', ascending=False).iloc[0]['n_compounds']
    check(top1 == 174, f"Largest scaffold series == 174  (got {top1})")
except FileNotFoundError as e:
    errors.append(f"  FAIL: {e}")
    print(f"  [FAIL] {e}")

# ── Activity cliffs ───────────────────────────────────────────────────────────
print("\n[4] activity_cliffs.parquet")
try:
    ac = pd.read_parquet(ROOT / 'publication/data/activity_cliffs.parquet')
    sev = ac[ac['cliff_tier'] == 'severe']
    check(len(sev) == CANON['n_severe_cliffs'],
          f"n_severe_cliffs == {CANON['n_severe_cliffs']}  (got {len(sev)})")

    sev_iks = set(sev['inchi_key_a'].tolist() + sev['inchi_key_b'].tolist())
    check(len(sev_iks) == 99,
          f"Severe cliff compounds == 99  (got {len(sev_iks)})")

    # ecfp4_only_cliff column
    if 'ecfp4_only_cliff' in ac.columns:
        n_ecfp4_only = sev['ecfp4_only_cliff'].sum()
        check(n_ecfp4_only == 13,
              f"ecfp4_only_cliff == 13  (got {n_ecfp4_only})")
    else:
        errors.append("  FAIL: ecfp4_only_cliff column missing")
        print("  [FAIL] ecfp4_only_cliff column missing")

    # Max delta
    max_delta = sev['delta_pic50'].abs().max()
    check(abs(max_delta - 3.045) < 0.01,
          f"Max |ΔpIC50| ≈ 3.045  (got {max_delta:.3f})")

except FileNotFoundError as e:
    errors.append(f"  FAIL: {e}")
    print(f"  [FAIL] {e}")

# ── MMP pairs ────────────────────────────────────────────────────────────────
print("\n[5] mmp_pairs_cliff99.csv")
try:
    mmp = pd.read_csv(ROOT / 'outputs/mmp/mmp_pairs_cliff99.csv')
    check(len(mmp) == 707, f"MMP pairs == 707  (got {len(mmp)})")
    n_cores = mmp['shared_core'].nunique()
    check(n_cores == 24, f"Unique shared cores == 24  (got {n_cores})")
except FileNotFoundError as e:
    errors.append(f"  FAIL: {e}")
    print(f"  [FAIL] {e}")

# ── Final verdict ─────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
if errors:
    print(f"VALIDATION FAILED — {len(errors)} check(s) failed:")
    for e in errors:
        print(e)
    print("\nDo NOT regenerate figures until all checks pass.")
    sys.exit(1)
else:
    print("ALL CHECKS PASSED — canonical numbers confirmed.")
    print(f"  Compounds:        {CANON['n_compounds']}")
    print(f"  Severe cliffs:    {CANON['n_severe_cliffs']}")
    print(f"  Hub compounds:    {CANON['n_hubs']}")
    print(f"  Murcko scaffolds: {CANON['n_scaffolds']}")
    print("\nProceed with figure generation.")
