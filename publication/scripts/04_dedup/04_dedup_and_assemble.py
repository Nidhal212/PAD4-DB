#!/usr/bin/env python3
"""
PAD4-DB v2 — Step 04: Cross-AID Dedup + Compound-Level Assembly
scripts/04_dedup/04_dedup_and_assemble.py

Inputs:
  data/interim/normalized/potency_space.parquet      (7,319 rows)
  data/interim/normalized/hts_space.parquet          (332,368 rows)
  data/interim/standardized/standardized_compounds.parquet

Outputs:
  data/interim/normalized/dedup_aid_level.parquet    (InChIKey × source × endpoint_type)
  data/processed/pad4_compounds.parquet              (InChIKey — one row per compound)
  data/processed/hts_compound_index.parquet          (InChIKey — all unique HTS compounds)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────
POT_PATH   = Path("data/interim/normalized/potency_space.parquet")
HTS_PATH   = Path("data/interim/normalized/hts_space.parquet")
STD_PATH   = Path("data/interim/standardized/standardized_compounds.parquet")
NORM_PATH  = Path("data/interim/normalized/normalized_activities.parquet")

OUT_AID    = Path("data/interim/normalized/dedup_aid_level.parquet")
OUT_COMP   = Path("data/processed/pad4_compounds.parquet")
OUT_HTS    = Path("data/processed/hts_compound_index.parquet")

SRC_GROUP  = ["inchi_key", "source", "endpoint_type"]


# =============================================================================
# parse_args
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="PAD4-DB Step 04: Cross-AID Dedup + Assembly")
    p.add_argument("--dry-run", action="store_true",
                   help="Use first 500 potency rows + 1000 HTS rows")
    p.add_argument("--debug", action="store_true",
                   help="Print group sizes at each aggregation step")
    return p.parse_args()


# =============================================================================
# load_inputs
# =============================================================================

def load_inputs(dry_run: bool, debug: bool):
    for path in [POT_PATH, HTS_PATH, STD_PATH, NORM_PATH]:
        if not path.exists():
            sys.exit(f"ERROR: required input not found: {path}")

    pot = pd.read_parquet(POT_PATH)
    hts = pd.read_parquet(HTS_PATH)
    std = pd.read_parquet(STD_PATH, columns=["inchi_key", "smiles_std", "n_heavy_atoms", "mol_weight"])

    if dry_run:
        pot = pot.head(500).copy()
        hts = hts.head(1000).copy()
        print(f"DRY-RUN: potency={len(pot)} rows, hts={len(hts)} rows")

    if debug:
        print(f"Loaded potency_space: {len(pot):,} rows")
        print(f"Loaded hts_space:     {len(hts):,} rows")

    return pot, hts, std


# =============================================================================
# apply_aid_preferred_filter
# =============================================================================

def apply_aid_preferred_filter(df: pd.DataFrame, debug: bool, dry_run: bool = False):
    """
    Drop rows where aid != aid_preferred (AID 2202576/77 overlap dedup).
    Returns (filtered_df, aid_preferred_iks) where aid_preferred_iks is the
    set of InChIKeys that had a preference applied in the original df.
    """
    # Track BEFORE filtering — these InChIKeys had the preference applied
    aid_preferred_iks = set(df[df["aid"] != df["aid_preferred"]]["inchi_key"])

    filtered = df[df["aid"] == df["aid_preferred"]].copy()

    if not dry_run:
        assert 7_200 <= len(filtered) <= 7_319, (
            f"ASSERTION FAILED: post-filter row count {len(filtered)} outside [7200, 7319]"
        )

    dropped = len(df) - len(filtered)
    if debug or dropped > 0:
        print(f"apply_aid_preferred_filter: dropped {dropped} rows "
              f"({len(aid_preferred_iks)} InChIKeys affected) → {len(filtered):,} remain")

    return filtered, aid_preferred_iks


# =============================================================================
# aggregate_to_source_level  →  Output 1
# =============================================================================

def aggregate_to_source_level(df: pd.DataFrame, aid_preferred_iks: set,
                               debug: bool) -> pd.DataFrame:
    """
    Group by InChIKey × source × endpoint_type.
    All aggregation done in LOG SPACE (pic50_aid is already pIC50 = -log10(IC50_M)).
    """
    # ── Numeric aggregations ────────────────────────────────────────────────
    base = df.groupby(SRC_GROUP, sort=False).agg(
        n_aids=("aid", "nunique"),
        n_total_measurements=("n_replicates", "sum"),
        _pic50_mean=("pic50_aid", "mean"),
        _pic50_std=("pic50_aid", "std"),   # ddof=1; NaN for n=1 groups
        _pic50_max=("pic50_aid", "max"),
        _pic50_min=("pic50_aid", "min"),
        value_nM_min=("value_nM_min", "min"),
        value_nM_max=("value_nM_max", "max"),
        use_in_potency_model=("use_in_potency_model", "any"),
    ).reset_index()

    base["log_value_mean"]    = base["_pic50_mean"]
    base["log_value_std"]     = base["_pic50_std"]
    base["pic50_source"]      = base["log_value_mean"]
    base["cross_aid_spread"]  = base["_pic50_max"] - base["_pic50_min"]
    base["high_variance_aid"] = base["cross_aid_spread"] > 1.0
    base = base.drop(columns=["_pic50_mean", "_pic50_std", "_pic50_max", "_pic50_min"])

    # ── String / categorical aggregations ───────────────────────────────────
    aid_agg = (
        df.groupby(SRC_GROUP, sort=False)["aid"]
        .agg(lambda x: "|".join(sorted(x.unique())))
        .reset_index()
        .rename(columns={"aid": "aid_list"})
    )

    mech_agg = (
        df.groupby(SRC_GROUP, sort=False)["assay_mechanism_class"]
        .agg(lambda x: x.mode().iloc[0])
        .reset_index()
    )

    out1 = base.merge(aid_agg, on=SRC_GROUP).merge(mech_agg, on=SRC_GROUP)

    # ── aid_preferred_used ──────────────────────────────────────────────────
    out1["aid_preferred_used"] = out1["inchi_key"].isin(aid_preferred_iks)

    # ── Assertions ──────────────────────────────────────────────────────────
    assert len(out1) <= len(df), (
        f"ASSERTION FAILED: source-level aggregation expanded rows ({len(out1)} > {len(df)})"
    )

    if debug:
        print(f"aggregate_to_source_level: {len(out1):,} groups")
        print(f"  n_aids distribution: {out1.n_aids.value_counts().to_dict()}")
        print(f"  high_variance_aid: {out1.high_variance_aid.sum()}")

    COLS = [
        "inchi_key", "source", "endpoint_type", "assay_mechanism_class",
        "n_aids", "n_total_measurements", "pic50_source", "log_value_mean",
        "log_value_std", "value_nM_min", "value_nM_max", "cross_aid_spread",
        "high_variance_aid", "aid_list", "aid_preferred_used", "use_in_potency_model",
    ]
    return out1[COLS]


# =============================================================================
# join_mol_properties
# =============================================================================

def join_mol_properties(df: pd.DataFrame, std_df: pd.DataFrame) -> pd.DataFrame:
    """Join smiles_std, n_heavy_atoms, mol_weight from standardized_compounds."""
    lut = (
        std_df[std_df["smiles_std"].notna()]
        .drop_duplicates("inchi_key")
        .set_index("inchi_key")[["smiles_std", "n_heavy_atoms", "mol_weight"]]
    )
    df["smiles_std"]    = df["inchi_key"].map(lut["smiles_std"])
    df["n_heavy_atoms"] = df["inchi_key"].map(lut["n_heavy_atoms"])
    df["mol_weight"]    = df["inchi_key"].map(lut["mol_weight"])
    return df


# =============================================================================
# aggregate_to_compound_level  →  Output 2
# =============================================================================

def aggregate_to_compound_level(out1: pd.DataFrame, std_df: pd.DataFrame,
                                 debug: bool, dry_run: bool = False) -> pd.DataFrame:
    """
    Group by InChIKey. Input = Output 1 filtered to use_in_potency_model == True.
    All aggregation in LOG SPACE (pic50_source is pIC50 units).
    """
    src_in = out1[out1["use_in_potency_model"] == True].copy()

    if debug:
        print(f"aggregate_to_compound_level: input {len(src_in):,} source-level rows")

    base = src_in.groupby("inchi_key", sort=False).agg(
        n_sources=("source", "nunique"),
        n_total_measurements=("n_total_measurements", "sum"),
        _pic50_mean=("pic50_source", "mean"),
        _pic50_std=("pic50_source", "std"),   # ddof=1; NaN for n=1 source
        _pic50_max=("pic50_source", "max"),
        _pic50_min=("pic50_source", "min"),
        use_in_potency_model=("use_in_potency_model", "any"),
    ).reset_index()

    base["log_value_mean_global"] = base["_pic50_mean"]
    base["log_value_std_global"]  = base["_pic50_std"]
    base["pic50_consensus"]       = base["log_value_mean_global"]
    base["source_spread"]         = base["_pic50_max"] - base["_pic50_min"]
    base["pic50_min"]             = base["_pic50_min"]
    base["pic50_max"]             = base["_pic50_max"]
    base["multi_source"]          = base["n_sources"] > 1
    base["concordant"]            = base["source_spread"] < 0.5
    base["discordant"]            = base["source_spread"] >= 1.0
    base["high_confidence"]       = base["multi_source"] & base["concordant"]
    base = base.drop(columns=["_pic50_mean", "_pic50_std", "_pic50_max", "_pic50_min"])

    # source_list and assay_mechanism_classes
    src_list = (
        src_in.groupby("inchi_key", sort=False)["source"]
        .agg(lambda x: "|".join(sorted(x.unique())))
        .reset_index()
        .rename(columns={"source": "source_list"})
    )
    mech_list = (
        src_in.groupby("inchi_key", sort=False)["assay_mechanism_class"]
        .agg(lambda x: "|".join(sorted(x.unique())))
        .reset_index()
        .rename(columns={"assay_mechanism_class": "assay_mechanism_classes"})
    )

    out2 = base.merge(src_list, on="inchi_key").merge(mech_list, on="inchi_key")

    # Join mol properties
    out2 = join_mol_properties(out2, std_df)

    # ── Assertions ──────────────────────────────────────────────────────────
    unique_compounds = len(out2)
    if not dry_run:
        assert 2_800 <= unique_compounds <= 3_093, (
            f"ASSERTION FAILED: compound count {unique_compounds} outside [2800, 3093]"
        )
        pic50_range_ok = out2["pic50_consensus"].between(2, 12).all()
        assert pic50_range_ok, (
            f"ASSERTION FAILED: pic50_consensus out of [2, 12] range\n"
            f"{out2[~out2['pic50_consensus'].between(2, 12)][['inchi_key','pic50_consensus']]}"
        )

    print(f"FINAL COMPOUND COUNT: {unique_compounds}")

    if debug:
        print(f"  multi_source: {out2.multi_source.sum()}")
        print(f"  high_confidence: {out2.high_confidence.sum()}")
        print(f"  concordant: {out2.concordant.sum()}")
        print(f"  discordant: {out2.discordant.sum()}")

    COLS = [
        "inchi_key", "smiles_std", "n_sources", "source_list", "n_total_measurements",
        "pic50_consensus", "log_value_mean_global", "log_value_std_global",
        "source_spread", "multi_source", "concordant", "discordant", "high_confidence",
        "assay_mechanism_classes", "use_in_potency_model",
        "n_heavy_atoms", "mol_weight", "pic50_min", "pic50_max",
    ]
    return out2[COLS]


# =============================================================================
# build_hts_index  →  Output 3
# =============================================================================

def build_hts_index(hts_df: pd.DataFrame, compound_inchikeys: set,
                    std_df: pd.DataFrame, debug: bool) -> pd.DataFrame:
    """
    Build HTS compound index with activity scoring.

    value_nM_min is null for all Pct_inh rows (% inhibition, not a concentration).
    Load normalized_activities.parquet to recover raw % inhibition values from value_raw.
    Activity threshold: mean % inhibition > 50% across replicates within an AID.
    """
    # ── Base: n_hts_assays per inchi_key ────────────────────────────────────
    hts_base = (
        hts_df.groupby("inchi_key", sort=False)["aid"]
        .nunique()
        .reset_index()
        .rename(columns={"aid": "n_hts_assays"})
    )

    # ── Load raw % inhibition from normalized_activities ────────────────────
    norm = pd.read_parquet(
        NORM_PATH,
        columns=["inchi_key", "aid", "layer", "endpoint_type", "value_raw", "qualifier", "norm_status"]
    )
    pct = norm[
        (norm["endpoint_type"] == "Pct_inh") &
        (norm["layer"] == "HTS") &
        (norm["norm_status"] == "OK")
    ].copy()

    pct["pct_inh"] = pd.to_numeric(pct["value_raw"], errors="coerce")
    pct = pct[pct["pct_inh"].notna()]

    if debug:
        print(f"build_hts_index: {len(pct):,} Pct_inh rows with parseable values")

    # ── max_pct_inh per compound ────────────────────────────────────────────
    max_pct = (
        pct.groupby("inchi_key", sort=False)["pct_inh"]
        .max()
        .reset_index()
        .rename(columns={"pct_inh": "max_pct_inh"})
    )

    # ── hts_activity_score: n AIDs with mean pct_inh > 50 (no qualifier) ───
    pct_unqualified = pct[pct["qualifier"].isin(["", "="])].copy()
    aid_mean = (
        pct_unqualified.groupby(["inchi_key", "aid"], sort=False)["pct_inh"]
        .mean()
        .reset_index()
    )
    active_per_ik = (
        aid_mean[aid_mean["pct_inh"] > 50]
        .groupby("inchi_key", sort=False)
        .size()
        .reset_index()
        .rename(columns={0: "hts_activity_score"})
    )

    # ── Merge all ────────────────────────────────────────────────────────────
    result = hts_base.merge(max_pct, on="inchi_key", how="left")
    result = result.merge(active_per_ik, on="inchi_key", how="left")
    result["hts_activity_score"] = result["hts_activity_score"].fillna(0).astype(int)
    result["any_active"]              = result["hts_activity_score"] > 0
    result["hts_consensus_confidence"] = (
        result["hts_activity_score"] / result["n_hts_assays"]
    )
    result["hts_outcome"] = result["any_active"].map({True: "Active", False: "Inactive"})
    result["confirmed_in_potency_space"] = result["inchi_key"].isin(compound_inchikeys)

    # ── Join smiles_std ───────────────────────────────────────────────────────
    smiles_lut = (
        std_df[std_df["smiles_std"].notna()]
        .drop_duplicates("inchi_key")
        .set_index("inchi_key")["smiles_std"]
    )
    result["smiles_std"] = result["inchi_key"].map(smiles_lut)

    # ── Assert ────────────────────────────────────────────────────────────────
    assert len(result) <= 327_336, (
        f"ASSERTION FAILED: HTS compound count {len(result)} > 327,336"
    )
    print(f"HTS UNIQUE COMPOUNDS: {len(result):,}")
    print(f"HTS confirmed in potency space: {result.confirmed_in_potency_space.sum():,}")

    COLS = [
        "inchi_key", "smiles_std", "n_hts_assays", "max_pct_inh",
        "any_active", "hts_activity_score", "hts_consensus_confidence",
        "hts_outcome", "confirmed_in_potency_space",
    ]
    return result[COLS]


# =============================================================================
# write_outputs
# =============================================================================

def write_outputs(out1: pd.DataFrame, out2: pd.DataFrame, out3: pd.DataFrame,
                  dry_run: bool):
    suffix = "_dryrun" if dry_run else ""

    def _path(p: Path) -> Path:
        return p.parent / (p.stem + suffix + p.suffix)

    for df, path in [(out1, _path(OUT_AID)), (out2, _path(OUT_COMP)), (out3, _path(OUT_HTS))]:
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(path, index=False)
        print(f"Written → {path}  ({len(df):,} rows, {df.shape[1]} cols)")


# =============================================================================
# print_summary
# =============================================================================

def print_summary(out1: pd.DataFrame, out2: pd.DataFrame, out3: pd.DataFrame):
    print("\n=== STEP 04 SUMMARY ===")

    print("\n--- dedup_aid_level (Output 1) ---")
    print(f"  Rows:             {len(out1):,}")
    print(f"  high_variance_aid:{out1.high_variance_aid.sum():,}")
    print(f"  aid_preferred_used:{out1.aid_preferred_used.sum():,}")
    print(f"  Sources:          {out1.source.value_counts().to_dict()}")
    print(f"  n_aids dist:      {out1.n_aids.value_counts().sort_index().to_dict()}")

    print("\n--- pad4_compounds (Output 2) ---")
    print(f"  Rows (compounds): {len(out2):,}")
    print(f"  pic50_consensus:  {out2.pic50_consensus.min():.3f} – {out2.pic50_consensus.max():.3f}")
    print(f"  multi_source:     {out2.multi_source.sum():,} / {len(out2):,}")
    print(f"  high_confidence:  {out2.high_confidence.sum():,} / {len(out2):,}")
    print(f"  concordant:       {out2.concordant.sum():,} / {len(out2):,}")
    print(f"  discordant:       {out2.discordant.sum():,} / {len(out2):,}")
    print(f"  smiles_std null:  {out2.smiles_std.isna().sum():,}")
    print(f"  assay_mech_classes: {out2.assay_mechanism_classes.value_counts().head(6).to_dict()}")

    print("\n--- hts_compound_index (Output 3) ---")
    print(f"  Rows:                       {len(out3):,}")
    print(f"  any_active:                 {out3.any_active.sum():,}")
    print(f"  confirmed_in_potency_space: {out3.confirmed_in_potency_space.sum():,}")
    print(f"  hts_outcome dist:           {out3.hts_outcome.value_counts().to_dict()}")
    print(f"  max_pct_inh null:           {out3.max_pct_inh.isna().sum():,}")


# =============================================================================
# main
# =============================================================================

def main():
    args = parse_args()

    pot, hts, std = load_inputs(args.dry_run, args.debug)

    # Step 1: AID 2202576/77 dedup
    pot_filtered, aid_preferred_iks = apply_aid_preferred_filter(pot, args.debug, args.dry_run)

    # Step 2: Source-level aggregation → Output 1
    print("\nAggregating to source level ...", flush=True)
    out1 = aggregate_to_source_level(pot_filtered, aid_preferred_iks, args.debug)

    # Step 3: Compound-level aggregation → Output 2
    print("Aggregating to compound level ...", flush=True)
    out2 = aggregate_to_compound_level(out1, std, args.debug, args.dry_run)

    # Step 4: HTS index → Output 3
    print("Building HTS compound index ...", flush=True)
    compound_inchikeys = set(out2["inchi_key"])
    out3 = build_hts_index(hts, compound_inchikeys, std, args.debug)

    # Write
    write_outputs(out1, out2, out3, args.dry_run)
    print_summary(out1, out2, out3)


if __name__ == "__main__":
    main()
