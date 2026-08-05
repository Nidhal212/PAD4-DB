"""
check_fig5_data.py — Data Sanity Check for Figure 5
Run this before fig05_cliff_network.py to ensure data integrity.
"""

import pandas as pd
from pathlib import Path

ROOT = Path('/home/nidhal/PAD4-db_V2')
print("=" * 60)
print("FIGURE 5 — DATA SANITY CHECK")
print("=" * 60)

# Replicating the exact Hub InChI keys from your figure script
HUB_IKS = {
    'A1': 'SMADULGDNOCLOP-GISFHXKWSA-N',
    'A2': 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
    'B1': 'UDCDEKJNAMHBFH-HSZRJFAPSA-N',
    'B2': 'DVCKJOQIVOGXEI-XMMPIXPASA-N',
}

# --- Helper: Auto-detect columns (matches the logic in your figure code) ---
def get_col(df, candidates):
    for col in candidates:
        if col in df.columns:
            return col
    return None

# ========================================
# CHECK 1: shared_assets.parquet (Nodes)
# ========================================
print("\n[1] Checking shared_assets.parquet (Node data)...")
path_shared = ROOT / 'data/interim/shared_assets.parquet'
if path_shared.exists():
    df_shared = pd.read_parquet(path_shared)
    print(f"  ✅ File found: {len(df_shared)} rows.")
    
    # Check columns
    mech_col = get_col(df_shared, ['mechanism_class', 'assay_mechanism_classes'])
    if mech_col:
        print(f"  ✅ Mechanism column found: '{mech_col}'")
        print(f"     - Values present: {df_shared[mech_col].nunique()} unique classes.")
    else:
        print(f"  ❌ FAIL: Could not find 'mechanism_class' or 'assay_mechanism_classes'.")
    
    # Check pIC50
    if 'pIC50' in df_shared.columns:
        print(f"  ✅ pIC50 column found. Range: {df_shared['pIC50'].min():.2f} - {df_shared['pIC50'].max():.2f}")
    else:
        print(f"  ❌ FAIL: 'pIC50' column missing from shared_assets.")
        
    # Check InChI Keys
    if 'inchi_key' in df_shared.columns:
        print(f"  ✅ 'inchi_key' column found.")
    else:
        print(f"  ❌ FAIL: 'inchi_key' column missing from shared_assets.")
else:
    print(f"  ❌ FAIL: File not found at {path_shared}")

# ========================================
# CHECK 2: activity_cliffs.parquet (Edges)
# ========================================
print("\n[2] Checking activity_cliffs.parquet (Cliff edges)...")
path_cliffs = ROOT / 'data/processed/activity_cliffs.parquet'
if path_cliffs.exists():
    df_cliffs = pd.read_parquet(path_cliffs)
    print(f"  ✅ File found: {len(df_cliffs)} total rows.")
    
    # Check tier column
    tier_col = get_col(df_cliffs, ['cliff_tier', 'Cliff_Tier'])
    if tier_col:
        print(f"  ✅ Tier column found: '{tier_col}'")
        counts = df_cliffs[tier_col].value_counts()
        print(f"     - Severe count: {counts.get('severe', 0)} (Expected canonical: 94)")
        if counts.get('severe', 0) != 94:
            print(f"     ⚠️  WARNING: Severe count does not match canonical 94!")
    else:
        print(f"  ❌ FAIL: Could not find 'cliff_tier' or 'Cliff_Tier'.")
        
    # Check edge columns
    req_cols = ['inchi_key_a', 'inchi_key_b', 'delta_pic50']
    missing = [c for c in req_cols if c not in df_cliffs.columns]
    if not missing:
        print(f"  ✅ Required edge columns found: {req_cols}")
        print(f"     - Delta pIC50 range: {df_cliffs['delta_pic50'].abs().min():.2f} - {df_cliffs['delta_pic50'].abs().max():.2f}")
    else:
        print(f"  ❌ FAIL: Missing columns: {missing}")
else:
    print(f"  ❌ FAIL: File not found at {path_cliffs}")

# ========================================
# CHECK 3: mmp_discontinuity_scores.csv (Degrees)
# ========================================
print("\n[3] Checking mmp_discontinuity_scores.csv (Degree data)...")
path_disc = ROOT / 'outputs/mmp/mmp_discontinuity_scores.csv'
if path_disc.exists():
    df_disc = pd.read_csv(path_disc)
    print(f"  ✅ File found: {len(df_disc)} rows.")
    
    if 'inchi_key' in df_disc.columns and 'severe_cliff_degree' in df_disc.columns:
        print(f"  ✅ Required columns found: 'inchi_key', 'severe_cliff_degree'")
    else:
        missing = [c for c in ['inchi_key', 'severe_cliff_degree'] if c not in df_disc.columns]
        print(f"  ❌ FAIL: Missing columns: {missing}")
else:
    print(f"  ❌ FAIL: File not found at {path_disc}")

# ========================================
# CHECK 4: Critical Hub Validation (A1, A2, B1, B2)
# ========================================
print("\n[4] Checking HUB Compounds (Critical check)...")
print("   These 4 InChI keys must exist across all data sources for the figure to be valid.")

# Check in shared_assets
if 'df_shared' in locals() and 'inchi_key' in df_shared.columns:
    print(f"  > In SHARED_ASSETS:")
    for label, ik in HUB_IKS.items():
        if ik in df_shared['inchi_key'].values:
            row = df_shared[df_shared['inchi_key'] == ik].iloc[0]
            print(f"      ✅ {label} found. pIC50: {row['pIC50']:.2f}, Mech: {row[mech_col] if mech_col else 'N/A'}")
        else:
            print(f"      ❌ {label} MISSING from shared_assets!")

# Check in activity_cliffs (do they participate in severe cliffs?)
if 'df_cliffs' in locals() and 'inchi_key_a' in df_cliffs.columns:
    print(f"\n  > In ACTIVITY_CLIFFS (Severe edges):")
    # Filter severe first if possible
    if tier_col:
        severe_df = df_cliffs[df_cliffs[tier_col].str.lower() == 'severe']
        all_severe_nodes = set(severe_df['inchi_key_a'].tolist() + severe_df['inchi_key_b'].tolist())
    else:
        all_severe_nodes = set(df_cliffs['inchi_key_a'].tolist() + df_cliffs['inchi_key_b'].tolist())
        
    for label, ik in HUB_IKS.items():
        if ik in all_severe_nodes:
            print(f"      ✅ {label} found in severe cliff network.")
        else:
            print(f"      ❌ {label} MISSING from severe cliff network!")

# Check in discontinuity scores (get their degree)
if 'df_disc' in locals() and 'inchi_key' in df_disc.columns:
    print(f"\n  > In DISCONTINUITY SCORES (Degrees):")
    degree_map = dict(zip(df_disc['inchi_key'], df_disc['severe_cliff_degree']))
    for label, ik in HUB_IKS.items():
        deg = degree_map.get(ik)
        if deg is not None:
            print(f"      ✅ {label} found. Degree: {deg}")
        else:
            print(f"      ❌ {label} MISSING from discontinuity scores!")

print("\n" + "=" * 60)
print("DATA CHECK COMPLETE.")
print("If you see all '✅' marks, your data is perfectly ready to plot!")
print("If you see '❌' or '⚠️', check your data pipeline before running the plot script.")
print("=" * 60)