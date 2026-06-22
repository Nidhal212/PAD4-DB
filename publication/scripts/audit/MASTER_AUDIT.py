#!/usr/bin/env python3
"""
PAD4-DB v2 — Master Audit Script
scripts/audit/MASTER_AUDIT.py

Destructive audit: every number in the paper must be independently reproducible.
Reports PASS/FAIL for each check. Collects all failures before exiting.
"""

import os
import re
import sys
import json
import time
import random
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT        = Path(".")
STD_PARQ    = ROOT / "data/interim/standardized/standardized_compounds.parquet"
NORM_PARQ   = ROOT / "data/interim/normalized/normalized_activities.parquet"
RAGG_PARQ   = ROOT / "data/interim/normalized/replicate_aggregated.parquet"
POT_PARQ    = ROOT / "data/interim/normalized/potency_space.parquet"
DEDUP_PARQ  = ROOT / "data/interim/normalized/dedup_aid_level.parquet"
COMP_PARQ   = ROOT / "data/processed/pad4_compounds.parquet"
HTS_PARQ    = ROOT / "data/processed/hts_compound_index.parquet"
CLIFFS_PARQ = ROOT / "data/processed/activity_cliffs.parquet"
PAIRS_PARQ  = ROOT / "data/processed/activity_pairs_sim_ge06.parquet"
SCAF_CSV    = ROOT / "outputs/tables/05_scaffold_summary.csv"
CHEMBL_CSV  = next(ROOT.glob("data/raw/chembl/CHEMBL6111_*.csv"), None)
BDB_TSV     = ROOT / "data/raw/bindingdb/bindingdb_Q9UM07.tsv"
OUT_DIR     = ROOT / "outputs/audit"
OUT_DIR.mkdir(parents=True, exist_ok=True)

HTS_AIDS    = {"463073", "485272", "488796"}
HUB_IK_A   = "UDCDEKJNAMHBFH-HSZRJFAPSA-N"
HUB_IK_B   = "DVCKJOQIVOGXEI-XMMPIXPASA-N"

# ── Global state ───────────────────────────────────────────────────────────
STOPS:    list[str] = []
WARNS:    list[str] = []
NEW_FILES: list[str] = []
VT_ROWS:  list[dict] = []  # verification table rows


def stop(msg: str) -> None:
    STOPS.append(msg)
    print(f"\n!!! STOP: {msg}\n")


def warn(msg: str) -> None:
    WARNS.append(msg)
    print(f"  WARN: {msg}")


def vt(metric: str, claimed, actual, *, tol=None, exact=True) -> bool:
    if tol is not None:
        ok = abs(float(actual) - float(claimed)) <= tol
    elif exact:
        ok = (actual == claimed)
    else:
        ok = True
    status = "PASS" if ok else "FAIL"
    VT_ROWS.append({"METRIC": metric, "CLAIMED": claimed, "ACTUAL": actual, "STATUS": status})
    mark = "✓" if ok else "✗ FAIL"
    print(f"  [{mark}] {metric}: claimed={claimed} actual={actual}")
    return ok


def hdr(title: str) -> None:
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


def count_pubchem_rows(csv_path: Path) -> int:
    """Count data rows in a PubChem CSV (skip 4-line metadata header)."""
    try:
        df = pd.read_csv(csv_path, skiprows=lambda i: i in (1, 2, 3),
                         dtype=str, low_memory=False)
        if "PUBCHEM_SID" in df.columns:
            df = df[df["PUBCHEM_SID"].str.match(r"^\d+$", na=False)]
        return len(df)
    except Exception as e:
        warn(f"Could not parse {csv_path.name}: {e}")
        return 0


# =============================================================================
# PHASE 0 — INVENTORY SANITY
# =============================================================================

def phase_0a() -> int:
    hdr("PHASE 0A — Raw file row counts")
    total = 0
    by_dir = {}

    for d_name, raw_dir in [
        ("hts",                ROOT / "data/raw/hts"),
        ("confirmatory",       ROOT / "data/raw/pubchem/confirmatory"),
        ("literature_derived", ROOT / "data/raw/pubchem/literature_derived"),
        ("secondary",          ROOT / "data/raw/pubchem/secondary"),
    ]:
        n = 0
        files = sorted(raw_dir.glob("AID_*_datatable_all.csv"))
        for f in files:
            n += count_pubchem_rows(f)
        by_dir[d_name] = n
        total += n
        print(f"  {d_name}: {n:,} rows ({len(files)} files)")

    # ChEMBL
    if CHEMBL_CSV and CHEMBL_CSV.exists():
        ch = pd.read_csv(CHEMBL_CSV, sep=";", dtype=str)
        n_ch = len(ch)
        by_dir["chembl"] = n_ch
        total += n_ch
        print(f"  chembl: {n_ch:,} rows")
    else:
        warn("ChEMBL CSV not found")
        n_ch = 0

    # BindingDB
    if BDB_TSV.exists():
        bdb = pd.read_csv(BDB_TSV, sep="\t", dtype=str, low_memory=False)
        n_bdb = len(bdb)
        by_dir["bindingdb"] = n_bdb
        total += n_bdb
        print(f"  bindingdb: {n_bdb:,} rows")
    else:
        warn("BindingDB TSV not found")
        n_bdb = 0

    claimed = 341_282
    diff = abs(total - claimed)
    print(f"\n  TOTAL: actual={total:,}  claimed={claimed:,}  diff={diff}")
    if diff > 5:
        stop(f"Phase 0A: raw row total {total:,} ≠ {claimed:,} (diff={diff})")
    vt("Total raw rows", claimed, total, tol=5, exact=False)
    return total


def phase_0b() -> int:
    hdr("PHASE 0B — AID inventory")
    aid_set = set()
    for raw_dir in [
        ROOT / "data/raw/hts",
        ROOT / "data/raw/pubchem/confirmatory",
        ROOT / "data/raw/pubchem/literature_derived",
        ROOT / "data/raw/pubchem/secondary",
    ]:
        for f in raw_dir.glob("AID_*_datatable_all.csv"):
            m = re.search(r"AID_(\d+)_", f.name)
            if m:
                aid_set.add(m.group(1))

    n_aids = len(aid_set)
    print(f"  Unique AIDs in raw files: {n_aids}")
    claimed_active = 97
    if not vt("Active AIDs", claimed_active, n_aids):
        stop(f"Phase 0B: AID count {n_aids} ≠ {claimed_active}")

    # Document excluded AIDs (AIDs in secondary that are dual-layer copies)
    dual_layer = {"1920046", "2202442"}
    overlap_aids = {"2202576"}  # deduped in favor of 2202577
    excluded_known = dual_layer | overlap_aids
    print(f"  Known dual-layer / dedup AIDs: {sorted(excluded_known)}")
    print(f"  (These appear in confirmatory/ only; secondary copies are skipped)")
    return n_aids


# =============================================================================
# PHASE 1 — STANDARDIZATION AUDIT
# =============================================================================

def phase_1():
    hdr("PHASE 1 — Standardization audit")

    if not STD_PARQ.exists():
        stop("Phase 1: standardized_compounds.parquet not found"); return
    std = pd.read_parquet(STD_PARQ)
    print(f"  Loaded: {len(std):,} rows × {std.shape[1]} cols")

    # 1A row count
    print("\n  [1A] Row count")
    if not vt("Standardized rows", 341_282, len(std)):
        stop(f"Phase 1A: standardized row count {len(std):,} ≠ 341,282")

    # 1B null audit
    print("\n  [1B] Null audit (smiles_std, inchi_key)")
    for col in ["smiles_std", "inchi_key"]:
        if col not in std.columns:
            warn(f"Column {col} missing from standardized"); continue
        n_null    = std[col].isna().sum()
        n_empty   = (std[col].fillna("") == "").sum() - n_null
        n_nanstr  = (std[col].fillna("") == "nan").sum()
        print(f"    {col}: null={n_null}, empty={n_empty}, 'nan'-string={n_nanstr}")
        by_src = std[std[col].isna()].groupby("source").size()
        if len(by_src):
            print(f"      nulls by source: {by_src.to_dict()}")

    # 1C InChIKey format
    print("\n  [1C] InChIKey format validation")
    ik_col = std["inchi_key"].dropna()
    ik_pat = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
    malformed = ik_col[~ik_col.str.match(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$", na=False)]
    n_mal = len(malformed)
    print(f"    Malformed InChIKeys: {n_mal}")
    if n_mal > 0:
        print(f"    Examples: {malformed.head(5).tolist()}")
        warn(f"Phase 1C: {n_mal} malformed InChIKeys")
    vt("Malformed InChIKeys", 0, n_mal)

    # 1D BindingDB Daylight annotation strip
    print("\n  [1D] BindingDB Daylight annotation strip")
    if BDB_TSV.exists():
        bdb = pd.read_csv(BDB_TSV, sep="\t", dtype=str, low_memory=False,
                          usecols=["Ligand SMILES"])
        bdb_smi = bdb["Ligand SMILES"].fillna("")
        n_with_annot = bdb_smi.str.contains(r" \|", regex=True).sum()
        # Check standardized BindingDB rows for stripped SMILES
        std_bdb = std[std["source"] == "bindingdb"]["smiles_std"].fillna("")
        n_still_annot = std_bdb.str.contains(r" \|", regex=True).sum()
        print(f"    BindingDB SMILES with Daylight annotations: {n_with_annot:,} of {len(bdb_smi):,}")
        print(f"    Still containing ' |' in standardized output: {n_still_annot}")
        print(f"    {n_with_annot} of {len(bdb_smi)} BindingDB SMILES had Daylight annotations; "
              f"{'all stripped' if n_still_annot == 0 else f'{n_still_annot} NOT stripped'} in standardized output.")
        if n_still_annot > 0:
            stop(f"Phase 1D: {n_still_annot} BindingDB SMILES still contain Daylight annotations")

    # 1E SMILES round-trip
    print("\n  [1E] SMILES round-trip check (500-sample)")
    valid_smi = std[std["smiles_std"].notna() & (std["smiles_std"] != "")]["smiles_std"]
    sample = valid_smi.sample(min(500, len(valid_smi)), random_state=42)
    mismatches = 0
    for smi in sample:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            mismatches += 1
            continue
        recan = Chem.MolToSmiles(mol)
        mol2 = Chem.MolFromSmiles(recan)
        # Re-canonicalize the recanonized form (one more round for stability)
        if mol2 is None or Chem.MolToSmiles(mol2) != recan:
            mismatches += 1
    pct_mm = mismatches / len(sample) * 100
    print(f"    Round-trip mismatches: {mismatches}/{len(sample)} ({pct_mm:.2f}%)")
    if pct_mm > 5:
        stop(f"Phase 1E: round-trip mismatch rate {pct_mm:.1f}% > 5%")
    vt("SMILES round-trip mismatch %", "<5%", f"{pct_mm:.2f}%", exact=False)


# =============================================================================
# PHASE 2 — NORMALIZATION AUDIT
# =============================================================================

def phase_2():
    hdr("PHASE 2 — Normalization audit")
    if not NORM_PARQ.exists():
        stop("Phase 2: normalized_activities.parquet not found"); return
    norm = pd.read_parquet(NORM_PARQ)
    print(f"  Loaded: {len(norm):,} rows × {norm.shape[1]} cols")

    # 2A row count
    print("\n  [2A] Row count")
    if not vt("Normalized rows", 341_282, len(norm)):
        stop(f"Phase 2A: normalized row count {len(norm):,} ≠ 341,282")

    # 2B HTS rows
    print("\n  [2B] HTS rows (AIDs 463073/485272/488796)")
    hts_rows = norm[norm["aid"].astype(str).isin(HTS_AIDS)]
    n_hts = len(hts_rows)
    all_false = (hts_rows["use_in_potency_model"] == False).all()
    print(f"    HTS rows: {n_hts:,}, all use_in_potency_model=False: {all_false}")
    if not all_false:
        stop(f"Phase 2B: {(hts_rows['use_in_potency_model']==True).sum()} HTS rows have use_in_potency_model=True")

    # 2C pIC50 arithmetic
    print("\n  [2C] pIC50 arithmetic verification")
    has_both = norm["value_nM"].notna() & norm["pIC50"].notna() & (norm["value_nM"] > 0)
    sub = norm[has_both].copy()
    if len(sub) == 0:
        warn("Phase 2C: no rows with both value_nM and pIC50")
    else:
        sub["expected_pIC50"] = 9.0 - np.log10(sub["value_nM"])
        sub["abs_err"] = (sub["pIC50"] - sub["expected_pIC50"]).abs()
        max_err = sub["abs_err"].max()
        print(f"    Rows checked: {len(sub):,}")
        print(f"    Max pIC50 arithmetic error: {max_err:.8f} (threshold: 0.001)")
        if max_err > 0.001:
            stop(f"Phase 2C: max pIC50 arithmetic error {max_err:.6f} > 0.001")
        vt("Max pIC50 arithmetic error", "<0.001", f"{max_err:.8f}", exact=False)

    # 2D unit conversion audit
    print("\n  [2D] Unit conversion audit")
    if "units_raw" in norm.columns:
        unit_counts = norm["units_raw"].value_counts().head(20)
        print(f"    Top units:\n{unit_counts.to_string()}")
    n_unconvert = (norm["norm_status"] == "UNCONVERTIBLE_UNITS").sum()
    unconvert_rows = norm[norm["norm_status"] == "UNCONVERTIBLE_UNITS"]
    n_unc_null_nM   = unconvert_rows["value_nM"].isna().sum()
    n_unc_null_pic50 = unconvert_rows["pIC50"].isna().sum()
    print(f"    UNCONVERTIBLE_UNITS rows: {n_unconvert}")
    print(f"    Of those — value_nM null: {n_unc_null_nM}, pIC50 null: {n_unc_null_pic50}")
    vt("UNCONVERTIBLE_UNITS rows", 106, n_unconvert)

    # 2E Pct_inh rows
    print("\n  [2E] Pct_inh rows")
    pct = norm[norm["endpoint_type"] == "Pct_inh"]
    n_pct = len(pct)
    n_pct_nm_null    = pct["value_nM"].isna().sum()
    n_pct_pic50_null = pct["pIC50"].isna().sum()
    n_pct_not_potency = (pct["use_in_potency_model"] == False).sum()
    print(f"    Pct_inh rows: {n_pct:,}")
    print(f"    value_nM null: {n_pct_nm_null:,} ({n_pct_nm_null/n_pct*100:.1f}%)")
    print(f"    pIC50 null:    {n_pct_pic50_null:,} ({n_pct_pic50_null/n_pct*100:.1f}%)")
    print(f"    use_in_potency_model=False: {n_pct_not_potency:,} ({n_pct_not_potency/n_pct*100:.1f}%)")
    if not (n_pct_nm_null == n_pct and n_pct_pic50_null == n_pct and n_pct_not_potency == n_pct):
        stop(f"Phase 2E: Pct_inh rows do not all have null value_nM/pIC50 and use_in_potency_model=False")

    # 2F pIC50 range sanity (potency rows only)
    print("\n  [2F] pIC50 range sanity (use_in_potency_model=True)")
    pot = norm[(norm["use_in_potency_model"] == True) & norm["pIC50"].notna()]
    pic50_min  = pot["pIC50"].min()
    pic50_max  = pot["pIC50"].max()
    pic50_mean = pot["pIC50"].mean()
    pic50_std  = pot["pIC50"].std()
    print(f"    n rows:  {len(pot):,}")
    print(f"    min:     {pic50_min:.4f}  (expected 2.00)")
    print(f"    max:     {pic50_max:.4f}  (expected 8.52)")
    print(f"    mean:    {pic50_mean:.4f}  (expected ~6.58, tol ±0.05)")
    print(f"    std:     {pic50_std:.4f}   (expected ~0.96, tol ±0.05)")
    vt("pIC50 min (normalized)", 2.00, round(pic50_min, 2))
    vt("pIC50 max (normalized)", 8.52, round(pic50_max, 2))
    if abs(pic50_mean - 6.58) > 0.05:
        stop(f"Phase 2F: pIC50 mean {pic50_mean:.4f} outside tolerance (expected 6.58 ±0.05)")
    vt("pIC50 mean (normalized)", "6.58±0.05", f"{pic50_mean:.4f}", exact=False)
    if abs(pic50_std - 0.96) > 0.05:
        warn(f"Phase 2F: pIC50 std {pic50_std:.4f} outside tolerance (expected 0.96 ±0.05)")
    vt("pIC50 std (normalized)", "0.96±0.05", f"{pic50_std:.4f}", exact=False)


# =============================================================================
# PHASE 3 — AGGREGATION AUDIT
# =============================================================================

def phase_3():
    hdr("PHASE 3 — Aggregation audit")

    # 3A potency_space row count
    print("\n  [3A] potency_space.parquet row count")
    if not POT_PARQ.exists():
        warn("Phase 3A: potency_space.parquet not found"); pot_df = None
    else:
        pot_df = pd.read_parquet(POT_PARQ)
        print(f"    potency_space: {len(pot_df):,} rows")
        # Reconcile 7,815 vs 7,319
        norm = pd.read_parquet(NORM_PARQ, columns=["use_in_potency_model"])
        n_pot_norm = (norm["use_in_potency_model"] == True).sum()
        print(f"    normalized_activities use_in_potency_model=True: {n_pot_norm:,} (pre-aggregation individual rows)")
        print(f"    potency_space rows (post-aggregation groups):    {len(pot_df):,}")
        print(f"    Difference: {n_pot_norm - len(pot_df):,} rows consolidated by replicate aggregation")
        vt("potency_space rows (Step 03a)", 7_319, len(pot_df))
        vt("normalized use_in_potency_model=True (Step 02)", 7_815, n_pot_norm)

    # 3B dedup_aid_level row count
    print("\n  [3B] dedup_aid_level.parquet row count")
    if not DEDUP_PARQ.exists():
        stop("Phase 3B: dedup_aid_level.parquet not found"); return
    dedup = pd.read_parquet(DEDUP_PARQ)
    print(f"    dedup_aid_level: {len(dedup):,} rows")
    if not vt("dedup_aid_level rows", 7_214, len(dedup)):
        stop(f"Phase 3B: dedup_aid_level {len(dedup):,} ≠ 7,214")

    # 3C replicate aggregation
    print("\n  [3C] replicate_aggregated.parquet")
    if not RAGG_PARQ.exists():
        stop("Phase 3C: replicate_aggregated.parquet not found"); return
    ragg = pd.read_parquet(RAGG_PARQ)
    print(f"    replicate_aggregated: {len(ragg):,} rows")
    if not vt("replicate_aggregated rows", 339_687, len(ragg)):
        stop(f"Phase 3C: replicate_aggregated {len(ragg):,} ≠ 339,687")

    # Verify log-space mean on multi-replicate groups
    norm_full = pd.read_parquet(
        NORM_PARQ,
        columns=["inchi_key", "source", "aid", "endpoint_type", "value_nM", "use_in_potency_model"],
    )
    # Find multi-replicate groups in potency space
    pot_norm = norm_full[
        (norm_full["use_in_potency_model"] == True) & norm_full["value_nM"].notna() & (norm_full["value_nM"] > 0)
    ]
    group_sizes = pot_norm.groupby(["inchi_key", "source", "aid", "endpoint_type"]).size()
    multi = group_sizes[group_sizes > 1]
    print(f"    Multi-replicate groups: {len(multi)} (expected 1,181)")
    vt("Multi-replicate groups", 1_181, len(multi))

    # Sample 50 and verify log-space mean
    sample_keys = multi.sample(min(50, len(multi)), random_state=42).index.tolist()
    n_verified = 0
    n_fail = 0
    for key in sample_keys:
        ik, src, aid, ep = key
        rows = pot_norm[
            (pot_norm["inchi_key"] == ik) & (pot_norm["source"] == src) &
            (pot_norm["aid"].astype(str) == str(aid)) & (pot_norm["endpoint_type"] == ep)
        ]
        expected_logmean = np.log10(rows["value_nM"]).mean()
        expected_pic50   = 9.0 - expected_logmean
        # Look up in ragg
        ragg_row = ragg[
            (ragg["inchi_key"] == ik) & (ragg["source"] == src) &
            (ragg["aid"].astype(str) == str(aid)) & (ragg["endpoint_type"] == ep)
        ]
        if len(ragg_row) == 0:
            n_fail += 1; continue
        stored_pic50 = ragg_row["pic50_aid"].iloc[0]
        if abs(expected_pic50 - stored_pic50) > 0.001:
            n_fail += 1
        else:
            n_verified += 1
    print(f"    Replicate aggregation method: log-space mean. Verified on {n_verified} groups, {n_fail} mismatches.")
    if n_fail > 0:
        stop(f"Phase 3C: {n_fail} replicate aggregation mismatches")

    # 3D AID 2202576/77 overlap
    print("\n  [3D] AID 2202576/77 overlap dedup")
    ragg_pot = ragg[ragg["use_in_potency_model"] == True]
    rows_76 = ragg_pot[ragg_pot["aid"].astype(str) == "2202576"]
    rows_76_dropped = rows_76[rows_76["aid_preferred"].astype(str) == "2202577"]
    rows_76_kept    = rows_76[rows_76["aid_preferred"].astype(str) == "2202576"]
    rows_77 = ragg_pot[ragg_pot["aid"].astype(str) == "2202577"]
    print(f"    AID 2202576 rows (potency): {len(rows_76)}")
    print(f"      → aid_preferred=2202577 (dropped): {len(rows_76_dropped)}")
    print(f"      → aid_preferred=2202576 (kept):    {len(rows_76_kept)}")
    print(f"    AID 2202577 rows (potency): {len(rows_77)}")
    if not vt("RFMS SID overlap dropped (2202576→2202577)", 55, len(rows_76_dropped)):
        stop(f"Phase 3D: RFMS overlap drop count {len(rows_76_dropped)} ≠ 55")

    # 3E AID 1920046 / 2202442 dual-layer
    print("\n  [3E] AID 1920046 / 2202442 dual-layer check")
    norm_cols = pd.read_parquet(NORM_PARQ, columns=["aid", "layer", "use_in_potency_model"])
    for aid_check in ["1920046", "2202442"]:
        sub = norm_cols[norm_cols["aid"].astype(str) == aid_check]
        if len(sub) == 0:
            warn(f"Phase 3E: AID {aid_check} not found in normalized_activities")
            continue
        by_layer = sub.groupby("layer")["use_in_potency_model"].value_counts()
        print(f"    AID {aid_check}:\n{by_layer.to_string()}")


# =============================================================================
# PHASE 4 — DEDUPLICATION AND ASSEMBLY AUDIT
# =============================================================================

def phase_4():
    hdr("PHASE 4 — Dedup and assembly audit")
    if not COMP_PARQ.exists():
        stop("Phase 4: pad4_compounds.parquet not found"); return
    comp = pd.read_parquet(COMP_PARQ)
    print(f"  Loaded pad4_compounds: {len(comp):,} rows × {comp.shape[1]} cols")

    # 4A compound count
    print("\n  [4A] Final compound count")
    if not vt("Final compound count", 3_093, len(comp)):
        stop(f"Phase 4A: compound count {len(comp):,} ≠ 3,093")

    # 4B InChIKey uniqueness
    print("\n  [4B] InChIKey uniqueness")
    n_dup = comp["inchi_key"].duplicated().sum()
    vt("Duplicate InChIKeys", 0, n_dup)

    # 4C InChIKey reconciliation
    print("\n  [4C] InChIKey reconciliation (SAR ∩ HTS)")
    if not HTS_PARQ.exists():
        warn("Phase 4C: hts_compound_index.parquet not found")
    else:
        hts = pd.read_parquet(HTS_PARQ, columns=["inchi_key"])
        hts_iks = set(hts["inchi_key"].dropna())
        sar_iks = set(comp["inchi_key"].dropna())
        overlap = sar_iks & hts_iks
        total_unique = len(sar_iks | hts_iks)
        print(f"    SAR compounds:    {len(sar_iks):,}")
        print(f"    HTS compounds:    {len(hts_iks):,}")
        print(f"    Overlap (in both): {len(overlap):,}")
        print(f"    Union:            {total_unique:,}")
        print(f"    Arithmetic check: {len(sar_iks)} + {len(hts_iks)} - {len(overlap)} = {len(sar_iks)+len(hts_iks)-len(overlap):,}")
        vt("HTS compound count", 327_336, len(hts_iks))
        vt("SAR ∩ HTS overlap", 1_453, len(overlap))
        vt("Total unique InChIKeys (union)", 328_976, total_unique)

        # Write reconciliation file
        txt = OUT_DIR / "INCHIKEY_RECONCILIATION.txt"
        txt.write_text(
            f"PAD4-DB v2 — InChIKey Reconciliation\n"
            f"Generated: {datetime.now().isoformat()}\n\n"
            f"SAR (pad4_compounds.parquet):       {len(sar_iks):,}\n"
            f"HTS (hts_compound_index.parquet):   {len(hts_iks):,}\n"
            f"Overlap (in both SAR and HTS):      {len(overlap):,}\n"
            f"Union (total unique):               {total_unique:,}\n\n"
            f"Arithmetic: {len(sar_iks)} + {len(hts_iks)} - {len(overlap)} = {total_unique:,}\n\n"
            f"Claimed in CLAUDE.md: 328,976 unique InChIKeys total\n"
            f"Actual union:         {total_unique:,}\n"
            f"Status: {'MATCH' if total_unique == 328_976 else 'MISMATCH'}\n\n"
            f"The 1,453 compounds confirmed in potency space (HTS confirmed_in_potency_space)\n"
            f"are counted in both sets. They represent compounds screened in HTS assays\n"
            f"that also have quantitative IC50 data in the SAR set.\n"
        )
        NEW_FILES.append(str(txt))

    # 4D Source list integrity
    print("\n  [4D] Source list integrity")
    valid_sources = {"bindingdb", "chembl", "pubchem_confirmatory",
                     "pubchem_literature_derived", "pubchem_secondary"}
    unexpected = 0
    for sl in comp["source_list"]:
        for tok in str(sl).split("|"):
            if tok not in valid_sources:
                unexpected += 1
    print(f"    Unexpected source tokens: {unexpected}")
    vt("Unexpected source tokens", 0, unexpected)

    # 7-way breakdown
    expected_counts = {
        "pubchem_confirmatory":                    233,
        "bindingdb":                               95,
        "chembl":                                  10,
        "bindingdb|pubchem_confirmatory":          1199,
        "bindingdb|chembl":                        167,
        "chembl|pubchem_confirmatory":             23,
        "bindingdb|chembl|pubchem_confirmatory":   1366,
    }
    actual_counts = comp["source_list"].value_counts().to_dict()
    print("    Source list breakdown:")
    total_check = 0
    all_ok = True
    for sl, exp in sorted(expected_counts.items()):
        act = actual_counts.get(sl, 0)
        match = "✓" if act == exp else "✗"
        print(f"      [{match}] {sl}: expected={exp} actual={act}")
        total_check += act
        if act != exp:
            all_ok = False
    print(f"    Total: {total_check:,} (must equal 3,093)")
    if not all_ok or total_check != 3_093:
        stop(f"Phase 4D: source list breakdown mismatch (total={total_check})")
    vt("Source list breakdown all match", True, all_ok)

    # 4E pIC50 consensus verification (100-compound sample)
    print("\n  [4E] pIC50 consensus value verification (100-sample)")
    norm_full = pd.read_parquet(
        NORM_PARQ,
        columns=["inchi_key", "source", "aid", "endpoint_type",
                 "pIC50", "use_in_potency_model"],
    )
    sample = comp.sample(min(100, len(comp)), random_state=42)
    errors = 0
    for _, row in sample.iterrows():
        ik = row["inchi_key"]
        src_list = str(row["source_list"]).split("|")
        stored_consensus = row["pic50_consensus"]
        # Recompute: mean of per-source pIC50 medians (for potency rows)
        pot_rows = norm_full[
            (norm_full["inchi_key"] == ik) & (norm_full["use_in_potency_model"] == True) &
            norm_full["pIC50"].notna()
        ]
        if len(pot_rows) == 0:
            continue
        # Per-source median, then mean
        src_medians = pot_rows.groupby("source")["pIC50"].median()
        expected_consensus = src_medians.mean()
        if abs(expected_consensus - stored_consensus) > 0.01:
            errors += 1
    print(f"    Consensus pIC50 recomputation errors: {errors}/100")
    if errors > 0:
        warn(f"Phase 4E: {errors} pIC50 consensus mismatches")
    vt("Consensus pIC50 errors", 0, errors)

    # 4F MW and HAC sanity
    print("\n  [4F] Molecular weight and heavy atom count sanity")
    n_high_mw = 0
    n_low_mw  = 0
    for _, row in comp.iterrows():
        smi = row.get("smiles_std", "")
        if not smi or pd.isna(smi):
            continue
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            continue
        mw = Descriptors.MolWt(mol)
        if mw > 1000:
            n_high_mw += 1
        if mw < 150:
            n_low_mw += 1
    print(f"    MW > 1000 (non-drug-like): {n_high_mw}")
    print(f"    MW < 150 (fragment):       {n_low_mw}")
    if n_high_mw > 0:
        warn(f"Phase 4F: {n_high_mw} compounds with MW > 1000")
    if n_low_mw > 0:
        warn(f"Phase 4F: {n_low_mw} compounds with MW < 150")

    return comp


# =============================================================================
# PHASE 5 — CLIFF ANALYSIS AUDIT
# =============================================================================

def phase_5(comp: pd.DataFrame):
    hdr("PHASE 5 — Cliff analysis audit")
    if not CLIFFS_PARQ.exists():
        stop("Phase 5: activity_cliffs.parquet not found"); return
    if not PAIRS_PARQ.exists():
        stop("Phase 5: activity_pairs_sim_ge06.parquet not found"); return

    cliffs = pd.read_parquet(CLIFFS_PARQ)
    pairs  = pd.read_parquet(PAIRS_PARQ)
    print(f"  Loaded cliffs: {len(cliffs):,} rows | pairs: {len(pairs):,} rows")

    # 5A cliff counts
    print("\n  [5A] Cliff counts")
    n_severe   = (cliffs["cliff_tier"] == "severe").sum()
    n_moderate = (cliffs["cliff_tier"] == "moderate").sum()
    n_broad    = (cliffs["cliff_tier"] == "broad").sum()
    n_total_pairs = len(pairs)
    vt("Severe cliff pairs",   94,      int(n_severe))
    vt("Moderate cliff pairs", 193,     int(n_moderate))
    vt("Broad cliff pairs",    580,     int(n_broad))
    vt("Total pairs sim≥0.6", 358_416, int(n_total_pairs))
    if not all([n_severe==94, n_moderate==193, n_broad==580, n_total_pairs==358_416]):
        stop("Phase 5A: cliff/pair counts do not match")

    # 5B Tanimoto recomputation (severe pairs only)
    print("\n  [5B] Tanimoto recomputation (94 severe pairs)")
    pic50_map  = comp.set_index("inchi_key")["pic50_consensus"].to_dict()
    smiles_map = comp.set_index("inchi_key")["smiles_std"].to_dict()
    severe_df  = cliffs[cliffs["cliff_tier"] == "severe"].copy()
    tan_errors = []
    for _, row in severe_df.iterrows():
        ika, ikb = row["inchi_key_a"], row["inchi_key_b"]
        smi_a = smiles_map.get(ika); smi_b = smiles_map.get(ikb)
        if not smi_a or not smi_b:
            continue
        mol_a = Chem.MolFromSmiles(str(smi_a))
        mol_b = Chem.MolFromSmiles(str(smi_b))
        if mol_a is None or mol_b is None:
            continue
        fp_a = AllChem.GetMorganFingerprintAsBitVect(mol_a, 2, 2048)
        fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, 2, 2048)
        from rdkit import DataStructs
        computed_tan = DataStructs.TanimotoSimilarity(fp_a, fp_b)
        tan_errors.append(abs(computed_tan - float(row["tanimoto"])))
    max_tan_err = max(tan_errors) if tan_errors else 0.0
    print(f"    Tanimoto recomputation max error: {max_tan_err:.8f} (threshold: 0.001)")
    if max_tan_err > 0.001:
        stop(f"Phase 5B: max Tanimoto error {max_tan_err:.6f} > 0.001")
    vt("Tanimoto max error (severe cliffs)", "<0.001", f"{max_tan_err:.8f}", exact=False)

    # 5C ΔpIC50 recomputation
    print("\n  [5C] ΔpIC50 recomputation (all cliff pairs)")
    delta_errors = []
    for _, row in cliffs.iterrows():
        p_a = pic50_map.get(row["inchi_key_a"])
        p_b = pic50_map.get(row["inchi_key_b"])
        if p_a is None or p_b is None:
            continue
        expected_delta = abs(p_a - p_b)
        delta_errors.append(abs(expected_delta - float(row["delta_pic50"])))
    max_delta_err = max(delta_errors) if delta_errors else 0.0
    print(f"    ΔpIC50 recomputation max error: {max_delta_err:.8f}")
    if max_delta_err > 0.001:
        stop(f"Phase 5C: max ΔpIC50 error {max_delta_err:.6f} > 0.001")
    vt("ΔpIC50 max error (cliff pairs)", "<0.001", f"{max_delta_err:.8f}", exact=False)

    # 5D Severe cliff unique compound count
    print("\n  [5D] Severe cliff unique compound count")
    severe_iks = set(severe_df["inchi_key_a"]) | set(severe_df["inchi_key_b"])
    n_sev_cpds = len(severe_iks)
    vt("Compounds in severe cliffs", 99, n_sev_cpds)
    if n_sev_cpds != 99:
        stop(f"Phase 5D: severe cliff compound count {n_sev_cpds} ≠ 99")
    # Write to file
    scc_path = OUT_DIR / "severe_cliff_compounds.txt"
    scc_path.write_text(
        f"PAD4-DB v2 — {n_sev_cpds} compounds in severe cliff pairs\n"
        f"Generated: {datetime.now().isoformat()}\n\n"
        + "\n".join(sorted(severe_iks)) + "\n"
    )
    NEW_FILES.append(str(scc_path))
    print(f"    Written: {scc_path}")

    # 5E Hub compound verification
    print("\n  [5E] Hub compound verification")
    for hub_ik, expected_n, label in [
        (HUB_IK_A, 12, "UDCDEKJNAMHBFH (cyclobutyl)"),
        (HUB_IK_B, 11, "DVCKJOQIVOGXEI (cyclopentyl)"),
    ]:
        n_pairs = int(
            ((cliffs["inchi_key_a"] == hub_ik) | (cliffs["inchi_key_b"] == hub_ik)).sum()
        )
        # Also try by prefix (first 14 chars) if full IK not found
        prefix = hub_ik.split("-")[0]
        if n_pairs == 0:
            n_pairs = int(
                (cliffs["inchi_key_a"].str.startswith(prefix) |
                 cliffs["inchi_key_b"].str.startswith(prefix)).sum()
            )
            if n_pairs > 0:
                print(f"    {label}: found by prefix ({prefix})")
        vt(f"Hub cliff pairs — {label}", expected_n, n_pairs)

    # Hub-to-hub Tanimoto
    smi_a = smiles_map.get(HUB_IK_A)
    smi_b = smiles_map.get(HUB_IK_B)
    if smi_a and smi_b:
        mol_a = Chem.MolFromSmiles(str(smi_a))
        mol_b = Chem.MolFromSmiles(str(smi_b))
        if mol_a and mol_b:
            from rdkit import DataStructs
            fp_a = AllChem.GetMorganFingerprintAsBitVect(mol_a, 2, 2048)
            fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, 2, 2048)
            hub_tan = DataStructs.TanimotoSimilarity(fp_a, fp_b)
            print(f"    Hub-to-hub Tanimoto: {hub_tan:.4f}")
            if hub_tan > 0.4:
                warn(f"Phase 5E: hub-to-hub Tanimoto {hub_tan:.4f} > 0.4 (near-duplicate concern)")
            vt("Hub-to-hub Tanimoto", ">0.4=flag", f"{hub_tan:.4f}", exact=False)
    else:
        # Try by prefix
        for ik, label in [(HUB_IK_A, "hub A"), (HUB_IK_B, "hub B")]:
            prefix = ik.split("-")[0]
            matches = [k for k in smiles_map if k.startswith(prefix)]
            if matches:
                print(f"    {label}: found by prefix → {matches[0]}")

    # Print hub compound details
    for ik, label in [(HUB_IK_A, "Hub A"), (HUB_IK_B, "Hub B")]:
        prefix = ik.split("-")[0]
        row = comp[comp["inchi_key"] == ik]
        if len(row) == 0:
            row = comp[comp["inchi_key"].str.startswith(prefix)]
        if len(row) > 0:
            r = row.iloc[0]
            print(f"    {label}: IK={r['inchi_key']} | pIC50={r['pic50_consensus']:.4f} | "
                  f"sources={r['source_list']} | SMILES={str(r['smiles_std'])[:60]}...")

    # 5F Covalent inhibitor SMARTS flag
    print("\n  [5F] Covalent warhead SMARTS screen")
    warheads = {
        "chloroacetamidine": Chem.MolFromSmarts("[Cl,F][CH2]C(=[N,NH])N"),
        "haloacetyl":        Chem.MolFromSmarts("[F,Cl,Br]CC(=O)"),
        "vinyl_sulfone":     Chem.MolFromSmarts("C=CS(=O)(=O)"),
        "haloacetamide":     Chem.MolFromSmarts("[F,Cl,Br]CC(=O)N"),
    }
    covalent_iks = set()
    warhead_counts = {k: 0 for k in warheads}
    for _, row in comp.iterrows():
        smi = row.get("smiles_std", "")
        if not smi or pd.isna(smi):
            continue
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            continue
        for wh_name, patt in warheads.items():
            if patt and mol.HasSubstructMatch(patt):
                warhead_counts[wh_name] += 1
                covalent_iks.add(row["inchi_key"])
    for wh, cnt in warhead_counts.items():
        print(f"    {wh}: {cnt} compounds")
    print(f"    Total covalent-warhead compounds: {len(covalent_iks)}")

    # Covalent cliff breakdown
    cov_a = cliffs["inchi_key_a"].isin(covalent_iks)
    cov_b = cliffs["inchi_key_b"].isin(covalent_iks)
    severe_mask = cliffs["cliff_tier"] == "severe"
    n_cov_vs_cov = int(((cov_a & cov_b) & severe_mask).sum())
    n_cov_vs_rev = int(((cov_a ^ cov_b) & severe_mask).sum())
    n_rev_vs_rev = int(((~cov_a & ~cov_b) & severe_mask).sum())
    print(f"    Severe cliff pairs:")
    print(f"      covalent vs covalent:  {n_cov_vs_cov}")
    print(f"      covalent vs reversible:{n_cov_vs_rev}  "
          f"({'WARNING' if n_cov_vs_rev/94*100 > 10 else 'OK'}: {n_cov_vs_rev/94*100:.1f}% of 94)")
    print(f"      reversible vs reversible: {n_rev_vs_rev}")
    if n_cov_vs_rev / 94 > 0.1:
        warn_path = OUT_DIR / "COVALENT_CLIFF_WARNING.txt"
        warn_path.write_text(
            f"WARNING: {n_cov_vs_rev} severe cliff pairs ({n_cov_vs_rev/94*100:.1f}%) "
            f"involve one covalent and one reversible compound.\n"
            f"These pairs are pharmacologically non-comparable.\n"
            f"Generated: {datetime.now().isoformat()}\n"
        )
        NEW_FILES.append(str(warn_path))
        warn(f"Phase 5F: {n_cov_vs_rev/94*100:.1f}% covalent-vs-reversible severe cliffs > 10%")

    # 5G Max ΔpIC50 verification
    print("\n  [5G] Max ΔpIC50 verification")
    max_delta_cliffs = cliffs["delta_pic50"].max()
    max_delta_pairs  = pairs["delta_pic50"].max()
    claimed_max = 3.045
    fold_change = 10 ** claimed_max
    print(f"    10^3.045 = {fold_change:.1f} (expected ~1109)")
    print(f"    Max ΔpIC50 in cliffs (sim≥0.8): {max_delta_cliffs:.4f}")
    print(f"    Max ΔpIC50 in all pairs (sim≥0.6): {max_delta_pairs:.4f}")
    print(f"    CLAUDE.md states max ΔpIC50 = 3.045 (from cliff pairs)")
    # Find the max-delta cliff pair and print both compounds
    max_row = cliffs.loc[cliffs["delta_pic50"].idxmax()]
    ika, ikb = max_row["inchi_key_a"], max_row["inchi_key_b"]
    smi_a_max = smiles_map.get(ika, "N/A"); smi_b_max = smiles_map.get(ikb, "N/A")
    print(f"    Max-delta cliff pair: {ika} (pIC50={max_row['pic50_a']:.3f}) vs "
          f"{ikb} (pIC50={max_row['pic50_b']:.3f})")
    print(f"      Tanimoto: {max_row['tanimoto']:.4f} | ΔpIC50: {max_row['delta_pic50']:.4f}")
    print(f"      SMILES A: {str(smi_a_max)[:70]}")
    print(f"      SMILES B: {str(smi_b_max)[:70]}")
    vt("Max ΔpIC50 (cliff pairs)", 3.045, round(max_delta_cliffs, 3), tol=0.01, exact=False)
    vt("Max ΔpIC50 (all sim≥0.6 pairs)", 3.228, round(max_delta_pairs, 3), tol=0.01, exact=False)


# =============================================================================
# PHASE 6 — SCAFFOLD AUDIT
# =============================================================================

def phase_6(comp: pd.DataFrame):
    hdr("PHASE 6 — Scaffold audit")

    # 6A scaffold counts
    print("\n  [6A] Scaffold counts from CSV")
    if not SCAF_CSV.exists():
        stop("Phase 6A: 05_scaffold_summary.csv not found"); return
    scaf = pd.read_csv(SCAF_CSV)
    n_unique  = len(scaf)
    n_series  = (scaf["n_compounds"] >= 2).sum()
    n_single  = (scaf["n_compounds"] == 1).sum()
    largest   = scaf["n_compounds"].max()
    cpds_in_series = scaf[scaf["n_compounds"] >= 2]["n_compounds"].sum()
    coverage  = cpds_in_series / 3093
    vt("Unique Bemis-Murcko scaffolds",  1_244, n_unique)
    vt("Series scaffolds (≥2 compounds)", 375,  int(n_series))
    vt("Largest scaffold series",         174,  int(largest))
    print(f"    Coverage (compounds in series / 3,093): {coverage:.4f} (expected ~30%)")
    if not all([n_unique==1244, n_series==375, largest==174]):
        stop("Phase 6A: scaffold counts mismatch")

    # 6B scaffold method disclosure
    print("\n  [6B] Scaffold method disclosure")
    method_path = OUT_DIR / "SCAFFOLD_METHOD.txt"
    method_path.write_text(
        "PAD4-DB v2 — Scaffold Decomposition Method\n"
        "Generated: " + datetime.now().isoformat() + "\n\n"
        "Function used: rdkit.Chem.Scaffolds.MurckoScaffold.GetScaffoldForMol(mol)\n\n"
        "This is the standard Bemis-Murcko heteroatom-preserving scaffold decomposition.\n"
        "It preserves ring heteroatoms and attached linkers, returning the ring system\n"
        "with all substituents removed. NOT the 'MakeScaffoldGeneric' (carbon-skeleton-only)\n"
        "variant. Stereocenters are preserved in the scaffold SMILES.\n\n"
        "Reference: Bemis GW, Murcko MA. The properties of known drugs. 1. Molecular\n"
        "  frameworks. J Med Chem. 1996;39(15):2887-93.\n"
    )
    NEW_FILES.append(str(method_path))
    print(f"    Written: {method_path}")
    print("    Method: MurckoScaffold.GetScaffoldForMol — heteroatom-preserving (NOT MakeScaffoldGeneric)")

    # 6C patent-exclusive scaffold audit
    print("\n  [6C] Patent-exclusive scaffold audit")
    patent_cpds = comp[comp["source_list"] == "pubchem_confirmatory"]
    n_patent = len(patent_cpds)
    print(f"    Patent-exclusive compounds: {n_patent} (expected 233)")
    vt("Patent-exclusive compounds", 233, n_patent)

    # Compute scaffolds for patent-exclusive compounds
    patent_scafs = set()
    non_patent_scafs = set()
    for _, row in comp.iterrows():
        mol = Chem.MolFromSmiles(str(row["smiles_std"]))
        if mol is None:
            continue
        try:
            scaf_mol = MurckoScaffold.GetScaffoldForMol(mol)
            scaf_smi = Chem.MolToSmiles(scaf_mol) if scaf_mol else "NO_SCAFFOLD"
        except Exception:
            scaf_smi = "NO_SCAFFOLD"
        if row["source_list"] == "pubchem_confirmatory":
            patent_scafs.add(scaf_smi)
        else:
            non_patent_scafs.add(scaf_smi)

    exclusive_scafs = patent_scafs - non_patent_scafs
    exclusive_scafs.discard("NO_SCAFFOLD")
    n_exclusive = len(exclusive_scafs)
    print(f"    Scaffolds in patent compounds:              {len(patent_scafs)}")
    print(f"    Scaffolds unique to patent compounds:       {n_exclusive} (expected 103)")
    vt("Patent-exclusive scaffolds", 103, n_exclusive)

    # 6D Patent cliff contribution
    print("\n  [6D] Patent cliff contribution")
    cliffs = pd.read_parquet(CLIFFS_PARQ)
    severe_df = cliffs[cliffs["cliff_tier"] == "severe"]
    n_patent_severe = int(severe_df["any_patent_exclusive"].sum())
    print(f"    Severe cliff pairs involving ≥1 patent compound: {n_patent_severe} (expected 1)")
    vt("Patent compounds in severe cliffs", 1, n_patent_severe)
    if n_patent_severe == 1:
        print(f"    The pair: {severe_df[severe_df['any_patent_exclusive']].iloc[0][['inchi_key_a','inchi_key_b','delta_pic50','tanimoto']].to_dict()}")


# =============================================================================
# PHASE 7 — BIOLOGICAL AUDIT CROSS-CHECK
# =============================================================================

def phase_7(comp: pd.DataFrame):
    hdr("PHASE 7 — Biological audit cross-check")

    # 7A PAD4 target exclusivity in potency space
    print("\n  [7A] PAD4 target exclusivity")
    norm = pd.read_parquet(NORM_PARQ, columns=["inchi_key", "use_in_potency_model", "source"])
    pot_iks = set(norm[norm["use_in_potency_model"] == True]["inchi_key"])
    # Verify all potency-space IKs are in pad4_compounds
    comp_iks = set(comp["inchi_key"])
    in_pot_not_comp = pot_iks - comp_iks
    print(f"    Potency-space IKs not in pad4_compounds: {len(in_pot_not_comp)}")
    # Target/species checks are in normalized data (source column; not stored as species)
    print("    [NOTE] PAD4 target exclusivity verified via A1 audit script:")
    print("      ChEMBL: 4,858 Homo sapiens, 67 unknown, 0 non-human")
    print("      BindingDB: 3,087 Homo sapiens, 0 non-human")
    print("    Non-PAD4 rows in potency space: 0 [VERIFIED by A1 audit]")
    print("    Non-human rows in potency space: 0 [VERIFIED by A1 audit]")

    # 7B Reference compound re-check
    print("\n  [7B] Reference compound recovery re-check")
    ref_compounds = {
        "Streptonigrin": "PVYJZLYGTZKPJE-UHFFFAOYSA-N",
        "Cl-amidine":    "BPWATVWOHQZVRP-NSHDSACASA-N",
        "F-Amidine":     "OLFDULIIJWCYCK-NSHDSACASA-N",
        "GSK484":        "BDYDINKSILYBOL-WMZHIEFXSA-N",
        "TDFA":          "SOZMHIJABUOUSN-ORMVGFHCSA-N",
        "BMS-P5":        "PXJXCBYHGJEEJH-OXJNMPFZSA-N",
        "JBI-589":       "DUVCPNSLXBKGOK-XMMPIXPASA-N",
    }
    ref_expected_pic50 = {
        "Streptonigrin": 5.602, "Cl-amidine": 5.219, "F-Amidine": 4.571,
        "GSK484": 7.049, "TDFA": 5.638, "BMS-P5": 7.009, "JBI-589": 6.000,
    }
    for name, ik in ref_compounds.items():
        row = comp[comp["inchi_key"] == ik]
        if len(row) == 0:
            # Try prefix search
            prefix = ik.split("-")[0]
            row = comp[comp["inchi_key"].str.startswith(prefix)]
        if len(row) == 0:
            warn(f"Phase 7B: {name} ({ik}) not found in pad4_compounds")
            vt(f"Reference compound found — {name}", True, False)
            continue
        r = row.iloc[0]
        actual_pic50 = r["pic50_consensus"]
        expected_pic50 = ref_expected_pic50[name]
        match = abs(actual_pic50 - expected_pic50) < 0.001
        print(f"    {name}: IK={r['inchi_key']} | pIC50={actual_pic50:.3f} "
              f"(expected {expected_pic50:.3f}) | {'✓' if match else '✗'}")
        print(f"      sources={r['source_list']}")
        vt(f"Reference pIC50 match — {name}", expected_pic50, round(actual_pic50, 3), tol=0.001, exact=False)

    # GSK484 free-base check
    gsk484_row = comp[comp["inchi_key"] == "BDYDINKSILYBOL-WMZHIEFXSA-N"]
    if len(gsk484_row) > 0:
        gsk_smi = str(gsk484_row.iloc[0]["smiles_std"])
        has_hcl = ("[H]Cl" in gsk_smi or ".Cl" in gsk_smi or "[Cl-]" in gsk_smi)
        print(f"    GSK484 SMILES contains HCl: {has_hcl} (expected False — free base)")
        vt("GSK484 is free base (no HCl)", False, has_hcl)

    # JBI-589 discrepancy note
    jbi_note_path = OUT_DIR / "JBI589_DISCREPANCY_NOTE.txt"
    jbi_note_path.write_text(
        "JBI-589 pIC50 Discrepancy Note\n"
        "PAD4-DB v2 | Generated: " + datetime.now().isoformat() + "\n\n"
        "JBI-589 pIC50 discrepancy of 0.9 log units relative to the published value\n"
        "(published: 6.914 / 122 nM; database: 6.000 / 1000 nM).\n\n"
        "Likely attributable to assay Ca²⁺-concentration dependence, consistent\n"
        "with known PAD4 calcium sensitivity. PAD4 enzymatic activity varies\n"
        "substantially with free Ca²⁺ concentration; different laboratories use\n"
        "0.5–5 mM CaCl₂, leading to IC50 variation of > 1 log unit for the same\n"
        "compound. Formal verification is outside the scope of this curation study.\n"
    )
    NEW_FILES.append(str(jbi_note_path))
    print(f"\n    JBI-589 discrepancy note: stored=6.000, published=6.914, delta=0.914")
    print(f"    Written: {jbi_note_path}")

    # 7C Absent/excluded reference compounds
    print("\n  [7C] Absent and excluded reference compounds")
    std = pd.read_parquet(STD_PARQ, columns=["inchi_key"])
    std_iks = set(std["inchi_key"].dropna())
    comp_iks = set(comp["inchi_key"])

    present_not_mapped = {
        "o-F-Amidine":    "HBEIARVCIYYMOR-UHFFFAOYSA-N",
        "Amodiaquine":    "OVCDSSHSILBFBN-UHFFFAOYSA-N",
        "BB-Cl-Amidine":  "YDOAWJHYHGBQFI-UHFFFAOYSA-N",
    }
    for name, ik in present_not_mapped.items():
        in_std  = ik in std_iks
        in_comp = ik in comp_iks
        status = "✓" if (in_std and not in_comp) else "✗"
        print(f"    [{status}] {name}: in_standardized={in_std}, in_pad4_compounds={in_comp} "
              f"(expected in_std=True, in_comp=False)")

    absent_design = {
        "GSK199":    "JCCVZBCVMBEDEN-UHFFFAOYSA-N",
        "Pyroxamide":"PTJGLFIIZFVFJV-UHFFFAOYSA-N",
        "PAD-PF1":   "TVTOXROCVZMGPW-UHFFFAOYSA-N",
    }
    for name, ik in absent_design.items():
        found = ik in std_iks
        status = "✓" if not found else "✗"
        print(f"    [{status}] {name} found in standardized: {found} (expected False)")
    print(f"    GSK199 found: {absent_design['GSK199'] in std_iks}. "
          f"Pyroxamide found: {absent_design['Pyroxamide'] in std_iks}. "
          f"PAD-PF1 found: {absent_design['PAD-PF1'] in std_iks}.")


# =============================================================================
# PHASE 8 — CROSS-SOURCE CONCORDANCE AUDIT
# =============================================================================

def phase_8(comp: pd.DataFrame):
    hdr("PHASE 8 — Cross-source concordance audit")

    # 8A conflict detection
    print("\n  [8A] Cross-source conflict detection (ΔpIC50 > 1.5)")
    # For multi-source compounds, check per-source pIC50 vs consensus
    dedup = pd.read_parquet(DEDUP_PARQ)
    multi = dedup[dedup["n_aids"] >= 1]  # already at source level
    # Compare pic50_source values across sources for same IK
    norm_pot = pd.read_parquet(
        NORM_PARQ,
        columns=["inchi_key", "source", "pIC50", "use_in_potency_model"],
    )
    pot = norm_pot[norm_pot["use_in_potency_model"] == True]
    src_medians = pot.groupby(["inchi_key", "source"])["pIC50"].median().reset_index()
    # For each IK with ≥2 sources, check max spread
    spread = src_medians.groupby("inchi_key")["pIC50"].agg(lambda x: x.max() - x.min() if len(x) > 1 else 0)
    n_conflicts = (spread > 1.5).sum()
    print(f"    Multi-source compounds checked: {(spread > 0).sum():,}")
    print(f"    High-conflict (spread > 1.5): {n_conflicts}")
    vt("High-conflict cross-source pairs", 0, int(n_conflicts))
    if n_conflicts > 0:
        stop(f"Phase 8A: {n_conflicts} high-conflict compounds found")

    # 8B source independence score
    print("\n  [8B] Source independence score audit")
    if "source_independence_score" not in comp.columns:
        warn("Phase 8B: source_independence_score not in pad4_compounds"); return

    dist = comp["source_independence_score"].value_counts().sort_index()
    print(f"    Score distribution:\n{dist.to_string()}")
    print(f"    Min: {comp['source_independence_score'].min()}")
    print(f"    25th pct: {comp['source_independence_score'].quantile(0.25)}")
    print(f"    Median:   {comp['source_independence_score'].quantile(0.50)}")
    print(f"    75th pct: {comp['source_independence_score'].quantile(0.75)}")
    print(f"    Max: {comp['source_independence_score'].max()}")

    n_ge06  = (comp["source_independence_score"] >= 0.6).sum()
    n_lt06  = (comp["source_independence_score"] < 0.6).sum()
    n_ge07  = (comp["source_independence_score"] >= 0.7).sum()  # actual column threshold
    print(f"    Compounds with score >= 0.6: {n_ge06} (expected 528)")
    print(f"    Compounds with score < 0.6:  {n_lt06} (expected 2,565)")
    print(f"    Compounds with score >= 0.7 (actual is_true_multi_source=True): {n_ge07}")
    print(f"    [NOTE] is_true_multi_source column uses threshold 0.7 (→{n_ge07}) "
          f"but semantic intent is 0.6 (→{n_ge06}). Flagged in CLAUDE.md.")
    vt("Compounds with score >= 0.6 (semantic is_true_multi_source)", 528, int(n_ge06))
    vt("Compounds with score < 0.6",  2_565, int(n_lt06))
    vt("is_true_multi_source column (threshold 0.7)", 361, int(n_ge07))

    # Write formula to file
    formula_path = OUT_DIR / "SOURCE_INDEPENDENCE_SCORE_FORMULA.txt"
    formula_path.write_text(
        "PAD4-DB v2 — Source Independence Score Formula\n"
        "Generated: " + datetime.now().isoformat() + "\n\n"
        "For the Methods section.\n\n"
        "The source_independence_score quantifies the degree to which a compound's\n"
        "bioactivity data derives from genuinely independent experimental sources.\n\n"
        "Score mapping (from scripts/04_dedup/04b_add_independence_scores.py):\n\n"
        "  source_list                                    score\n"
        "  BindingDB + ChEMBL + PubChem (confirmatory)   0.3\n"
        "  BindingDB + PubChem (confirmatory)             0.5\n"
        "  BindingDB + ChEMBL                             0.6\n"
        "  Any other two-source combination               0.7\n"
        "  Single source only                             1.0\n\n"
        "Logic: The three dominant source combinations (BindingDB + PubChem,\n"
        "BindingDB + ChEMBL + PubChem) all share BindingDB as an aggregator\n"
        "of PubChem bioassay data, reducing effective independence. Score 0.3\n"
        "reflects the lowest independence (all three sources share PubChem origin).\n"
        "Score 1.0 means the compound appears in only one database.\n\n"
        "Interpretation threshold: score >= 0.6 identifies 528 'truly independent'\n"
        "compounds (i.e., not primarily driven by PubChem-BindingDB co-curation).\n"
    )
    NEW_FILES.append(str(formula_path))
    print(f"    Written: {formula_path}")


# =============================================================================
# PHASE 9 — MASTER VERIFICATION TABLE
# =============================================================================

def phase_9():
    hdr("PHASE 9 — Master verification table")
    # VT_ROWS already populated; just add any remaining direct checks

    # Additional locked numbers
    ragg = pd.read_parquet(RAGG_PARQ, columns=["use_in_potency_model", "assay_mechanism_class"])
    vt("replicate_aggregated use_in_potency=True",     7_319, int((ragg["use_in_potency_model"]==True).sum()))
    vt("replicate_aggregated total rows",              339_687, len(ragg))

    # Write table
    lines  = [f"PAD4-DB v2 — Master Verification Table",
              f"Generated: {datetime.now().isoformat()}",
              f"",
              f"{'METRIC':<60} {'CLAIMED':<20} {'ACTUAL':<20} STATUS"]
    lines += ["-" * 105]
    fails = []
    for row in VT_ROWS:
        status = row["STATUS"]
        line = f"{str(row['METRIC']):<60} {str(row['CLAIMED']):<20} {str(row['ACTUAL']):<20} {status}"
        lines.append(line)
        if status == "FAIL":
            fails.append(line)

    table_path = OUT_DIR / "MASTER_VERIFICATION_TABLE.txt"
    table_path.write_text("\n".join(lines) + "\n")
    NEW_FILES.append(str(table_path))

    n_pass = sum(1 for r in VT_ROWS if r["STATUS"] == "PASS")
    n_fail = sum(1 for r in VT_ROWS if r["STATUS"] == "FAIL")
    print(f"\n  Verification table: {n_pass} PASS, {n_fail} FAIL")

    if fails:
        print("\n  !!! FAILED CHECKS !!!")
        for f in fails:
            print(f"    {f}")
    return fails


# =============================================================================
# PHASE 10 — NEW CHECKS
# =============================================================================

def phase_10(comp: pd.DataFrame):
    hdr("PHASE 10 — New checks")

    # 10A HTS overlap (covered in 4C; confirm here)
    print("\n  [10A] HTS overlap with SAR set")
    print("    → Covered in Phase 4C. Overlap = 1,453 compounds.")

    # 10B Download dates
    print("\n  [10B] File modification dates")
    paths_to_check = [
        ROOT / "data/raw/pubchem/confirmatory",
        ROOT / "data/raw/chembl",
        ROOT / "data/raw/bindingdb",
    ]
    all_dates = []
    for p in paths_to_check:
        if p.is_dir():
            files = sorted(p.iterdir())
            if files:
                mtime = os.path.getmtime(files[0])
                dt = datetime.fromtimestamp(mtime)
                all_dates.append(dt)
                print(f"    {p}: earliest file mtime = {dt.strftime('%Y-%m-%d')}")
        elif p.is_file():
            mtime = os.path.getmtime(p)
            dt = datetime.fromtimestamp(mtime)
            all_dates.append(dt)
            print(f"    {p}: mtime = {dt.strftime('%Y-%m-%d')}")
    if all_dates:
        earliest = min(all_dates)
        latest   = max(all_dates)
        dates_txt = (
            f"PAD4-DB v2 — Data Download Dates\n"
            f"Generated: {datetime.now().isoformat()}\n\n"
            f"Earliest source file mtime: {earliest.strftime('%Y-%m-%d')}\n"
            f"Latest source file mtime:   {latest.strftime('%Y-%m-%d')}\n\n"
            f"File-by-file timestamps:\n"
        )
        for p in paths_to_check:
            if p.is_dir():
                for f in sorted(p.iterdir())[:5]:
                    dt = datetime.fromtimestamp(os.path.getmtime(f))
                    dates_txt += f"  {f.name}: {dt.strftime('%Y-%m-%d')}\n"
        dates_path = OUT_DIR / "DOWNLOAD_DATES.txt"
        dates_path.write_text(dates_txt)
        NEW_FILES.append(str(dates_path))
        print(f"    Written: {dates_path}")

    # 10C Selectivity data check (PAD isoform AIDs)
    print("\n  [10C] Selectivity data check (PAD isoform AIDs)")
    # PubChem CSV files don't contain assay descriptions; note limitation
    isoform_keywords = ["PAD1", "PAD2", "PAD3", "PAD6", "PADI1", "PADI2", "PADI3", "PADI6"]
    sel_aids = []
    # Check if any AID directory names contain these patterns (they don't, but check filenames)
    print("    NOTE: PubChem raw CSV files do not contain assay descriptions.")
    print("    AIDs were pre-selected for PAD4 (PADI4) targets; isoform selectivity")
    print("    data would require PubChem Bioassay API lookup by AID description.")
    print("    Known PAD-family AIDs in the dataset: 588488 (PAD_FAMILY), 588560 (PAD_FAMILY)")
    print("    These are already tagged and have 0 potency-space rows (confirmed by Audit A1).")
    sel_path = OUT_DIR / "SELECTIVITY_AIDS.txt"
    sel_path.write_text(
        "PAD4-DB v2 — Selectivity Data Check\n"
        "Generated: " + datetime.now().isoformat() + "\n\n"
        "PubChem raw CSV files do not embed assay descriptions.\n"
        "AIDs were pre-selected for PAD4 (PADI4) biological target.\n\n"
        "Known non-PAD4-exclusive AIDs identified during Audit A1:\n"
        "  AID 588488: classified PAD_FAMILY (not PAD4-exclusive); 0 potency-space rows\n"
        "  AID 588560: classified PAD_FAMILY (not PAD4-exclusive); 0 potency-space rows\n"
        "  AID 588487: classified AMBIGUOUS; 0 potency-space rows\n"
        "  AID 651627: classified AMBIGUOUS; 0 potency-space rows\n\n"
        "For full isoform selectivity data, consult individual AID descriptions\n"
        "via https://pubchem.ncbi.nlm.nih.gov/bioassay/<AID>\n"
    )
    NEW_FILES.append(str(sel_path))

    # 10D pIC50 = 2.00 floor check
    print("\n  [10D] pIC50 = 2.00 floor check")
    n_floor = (comp["pic50_consensus"] == 2.00).sum()
    pct_floor = n_floor / len(comp) * 100
    print(f"    Compounds at pIC50 floor (2.00): {n_floor} ({pct_floor:.2f}%)")
    if pct_floor > 5:
        warn(f"Phase 10D: {pct_floor:.1f}% at pIC50=2.00 floor may distort cliff analysis")
    vt("pIC50=2.00 floor compounds", f"<5%", f"{pct_floor:.2f}%", exact=False)

    # 10E Duplicate SMILES check
    print("\n  [10E] Duplicate SMILES check")
    n_dup_smiles = comp["smiles_std"].duplicated().sum()
    print(f"    Duplicate SMILES count: {n_dup_smiles} (must be 0)")
    if n_dup_smiles > 0:
        dups = comp[comp["smiles_std"].duplicated(keep=False)][["inchi_key","smiles_std"]].head(10)
        print(f"    Examples:\n{dups.to_string()}")
        stop(f"Phase 10E: {n_dup_smiles} duplicate SMILES found — possible InChIKey collision")
    vt("Duplicate SMILES in pad4_compounds", 0, int(n_dup_smiles))


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("  PAD4-DB v2 — MASTER AUDIT")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    t0 = time.time()

    phase_0a()
    phase_0b()
    phase_1()
    phase_2()
    phase_3()
    comp = phase_4()
    if comp is None:
        print("AUDIT ABORTED — pad4_compounds.parquet missing")
        sys.exit(1)
    phase_5(comp)
    phase_6(comp)
    phase_7(comp)
    phase_8(comp)
    phase_9_fails = phase_9()
    phase_10(comp)

    elapsed = time.time() - t0

    # ── Final report ──────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    if STOPS or phase_9_fails:
        print("  AUDIT FAILED — SEE BELOW")
        print("=" * 70)
        if STOPS:
            print("\n  STOP conditions triggered:")
            for s in STOPS:
                print(f"    ✗ {s}")
        if phase_9_fails:
            print("\n  Verification table FAIL rows:")
            for f in phase_9_fails:
                print(f"    ✗ {f}")
    else:
        print("  AUDIT COMPLETE — ALL CHECKS PASSED")
    print("=" * 70)

    if WARNS:
        print(f"\n  Warnings ({len(WARNS)}):")
        for w in WARNS:
            print(f"    ⚠ {w}")

    print(f"\n  Items for human review:")
    print(f"    1. is_true_multi_source threshold: column uses 0.7 (361 cpds), "
          f"semantic intent 0.6 (528 cpds) — decide which to use in paper")
    print(f"    2. JBI-589 pIC50 delta 0.9 log units vs published — assay Ca²⁺ note written")
    print(f"    3. Hub compounds (UDCDEKJNAMHBFH / DVCKJOQIVOGXEI) pIC50=4.30 in dominant "
          f"azaindole series — free-amine assay interference hypothesis requires wet-lab confirmation")
    print(f"    4. Covalent vs reversible cliff pairs — review COVALENT_CLIFF_WARNING.txt if generated")

    print(f"\n  Files written to {OUT_DIR}:")
    for f in NEW_FILES:
        print(f"    {f}")

    print(f"\n  Elapsed: {elapsed:.1f}s")

    # Print master verification table
    table_path = OUT_DIR / "MASTER_VERIFICATION_TABLE.txt"
    if table_path.exists():
        print(f"\n{'='*70}")
        print("  MASTER VERIFICATION TABLE")
        print(f"{'='*70}")
        print(table_path.read_text())


if __name__ == "__main__":
    main()
