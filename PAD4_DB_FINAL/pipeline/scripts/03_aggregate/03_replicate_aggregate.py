#!/usr/bin/env python3
"""
PAD4-DB v2 — Step 03: Replicate Aggregation
scripts/03_aggregate/03_replicate_aggregate.py

Aggregates normalized_activities.parquet to grain:
  InChIKey × source × aid × endpoint_type

Joins smiles_std from standardized_compounds.parquet.
Output: data/interim/normalized/replicate_aggregated.parquet
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

NORM_PARQUET = Path("data/interim/normalized/normalized_activities.parquet")
STD_PARQUET  = Path("data/interim/standardized/standardized_compounds.parquet")
OUT_PATH     = Path("data/interim/normalized/replicate_aggregated.parquet")

# ── Snakemake integration ─────────────────────────────────────────────────
if "snakemake" in dir():
    _sm = snakemake  # noqa
    NORM_PARQUET = Path(_sm.input.normalized)
    STD_PARQUET  = Path(_sm.input.standardized)
    OUT_PATH     = Path(_sm.output.parquet)

GROUP_COLS = ["inchi_key", "source", "aid", "endpoint_type"]

# =============================================================================
# Assay mechanism classification
# =============================================================================

_COVALENT_IRREVERSIBLE_AIDS = {
    627371, 627428, 627432, 626724, 626728,
    1069613, 1069614, 1069618, 1069619, 1069623,
    626735, 626738, 712876, 725671, 1069608,
    1196521, 1364668, 1422898, 1422904, 2076402, 651867,
}

_FP_BINDING_AIDS = {
    1069597, 1069598, 1069599, 1069600, 1069601, 1069604, 1069605, 1069606,
    588487, 651627,
}

AID_ASSAY_CLASS: dict[int, str] = {
    # HTS screens
    463073: "HTS", 485272: "HTS", 488796: "HTS",
    # RFMS (Reactant-Free Mass Spectrometry)
    2202576: "RFMS", 2202577: "RFMS",
    2202596: "RFMS", 2202597: "RFMS", 2202717: "RFMS",
    # FP (Fluorescence Polarization) binding / confirmatory
    1919095: "FP", 1920200: "FP", 1920046: "FP", 2202442: "FP",
    2041348: "FP", 2041349: "FP", 2053867: "FP", 2053917: "FP",
    1806765: "FP", 1625405: "FP",
    # ABPP (Activity-Based Protein Profiling)
    1069605: "ABPP", 1069606: "ABPP", 588487: "ABPP",
    # MALDI mass shift
    1069597: "MALDI", 1069598: "MALDI", 1069599: "MALDI",
    1069600: "MALDI", 1069601: "MALDI", 1069604: "MALDI",
    # BAEE / Colorimetric enzymatic
    492970: "BAEE_Colorimetric", 588559: "BAEE_Colorimetric", 651627: "BAEE_Colorimetric",
    # Kinetics
    627432: "Kinetics",
    # Cellular
    2039667: "Cellular", 2053916: "Cellular", 2053915: "Cellular", 2193457: "Cellular",
}


def assign_assay_mechanism_class(aid_str: str, endpoint_type: str,
                                 layer: str, assay_class: str) -> str:
    """5-step priority: AID override → endpoint_type → layer → assay_class keywords → default."""
    try:
        aid = int(aid_str)
    except (ValueError, TypeError):
        aid = None

    # Step 1 — AID-specific overrides
    if aid is not None:
        if aid in _COVALENT_IRREVERSIBLE_AIDS:
            return "covalent_irreversible"
        if aid in _FP_BINDING_AIDS:
            return "fp_binding"

    # Step 2 — endpoint_type
    ep = (endpoint_type or "").strip()
    if ep == "Kinact_Ki":
        return "covalent_irreversible"
    if ep in ("Ki", "Kd"):
        return "fp_binding"
    if ep == "EC50":
        return "cellular"
    if ep == "Pct_inh":
        return "screening_single_conc"

    # Step 3 — layer
    if layer == "HTS":
        return "hts_screen"
    if layer == "E":
        return "cellular"
    if layer == "D":
        return "fp_binding"

    # Step 4 — assay_class keyword matching (case-insensitive)
    ac = (assay_class or "").upper()
    if "RFMS" in ac:
        return "rfms_enzymatic"
    if "FP" in ac or "FLUORESCENCE" in ac:
        return "fp_binding"
    if "BAEE" in ac or "COLORIMETRIC" in ac:
        return "baee_colorimetric"
    if "ABPP" in ac:
        return "fp_binding"
    if "MALDI" in ac or "MASS" in ac:
        return "fp_binding"
    if "CELLULAR" in ac or "WESTERN" in ac:
        return "cellular"
    if "KINETICS" in ac:
        return "covalent_irreversible"

    # Step 5 — default
    return "baee_colorimetric"


# =============================================================================
# Aggregation helpers
# =============================================================================

_STATUS_RANK = {
    "OK": 0,
    "NO_VALUE": 1,
    "QUALIFIER_ONLY": 2,
    "UNCONVERTIBLE_UNITS": 3,
    "PARSE_ERROR": 4,
}


def _worst_status(statuses: pd.Series) -> str:
    return max(statuses, key=lambda s: _STATUS_RANK.get(s, -1))


def _agg_qualifier(qualifiers: pd.Series) -> str:
    """Map '<=' → '<'; keep '<', '>', '~'; return most restrictive present."""
    mapped = set()
    for q in qualifiers:
        if q in ("<", "<="):
            mapped.add("<")
        elif q == ">":
            mapped.add(">")
        elif q == "~":
            mapped.add("~")
    return "<" if "<" in mapped else ">" if ">" in mapped else "~" if "~" in mapped else ""


# =============================================================================
# Main aggregation
# =============================================================================

def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    # ── Base aggregations (fast built-ins) ─────────────────────────────────
    base = df.groupby(GROUP_COLS, sort=False).agg(
        n_replicates=("value_nM", "size"),
        layer=("layer", "first"),
        aid_preferred=("aid_preferred", "first"),
        use_in_potency_model=("use_in_potency_model", "any"),
        value_nM_min=("value_nM", "min"),
        value_nM_max=("value_nM", "max"),
    ).reset_index()

    # ── Custom: worst norm_status ───────────────────────────────────────────
    status_agg = (
        df.groupby(GROUP_COLS, sort=False)["norm_status"]
        .agg(_worst_status)
        .reset_index()
        .rename(columns={"norm_status": "norm_status_worst"})
    )

    # ── Custom: any_qualifier ───────────────────────────────────────────────
    qual_agg = (
        df.groupby(GROUP_COLS, sort=False)["qualifier"]
        .agg(_agg_qualifier)
        .reset_index()
        .rename(columns={"qualifier": "any_qualifier"})
    )

    # ── Log-space stats (valid rows only) ──────────────────────────────────
    valid = df[df["value_nM"].notna() & (df["value_nM"] > 0)].copy()
    valid["log_nM"] = np.log10(valid["value_nM"])
    log_agg = (
        valid.groupby(GROUP_COLS, sort=False)["log_nM"]
        .agg(log_value_mean="mean", log_value_std="std")
        .reset_index()
    )

    # ── Merge ───────────────────────────────────────────────────────────────
    result = base.merge(status_agg, on=GROUP_COLS)
    result = result.merge(qual_agg, on=GROUP_COLS)
    result = result.merge(log_agg, on=GROUP_COLS, how="left")

    # ── Derived columns ─────────────────────────────────────────────────────
    result["pic50_aid"] = np.where(
        result["log_value_mean"].notna(),
        9.0 - result["log_value_mean"],
        np.nan,
    )
    result["hts_flag"] = result["layer"] == "HTS"

    # ── assay_mechanism_class ───────────────────────────────────────────────
    assay_cls_map = result["aid"].map(
        lambda a: AID_ASSAY_CLASS.get(int(a), "") if str(a).lstrip("-").isdigit() else ""
    )
    result["assay_mechanism_class"] = [
        assign_assay_mechanism_class(aid, ep, layer, ac)
        for aid, ep, layer, ac in zip(
            result["aid"], result["endpoint_type"], result["layer"], assay_cls_map
        )
    ]

    return result


def load_smiles_lookup() -> pd.Series:
    std = pd.read_parquet(STD_PARQUET, columns=["inchi_key", "smiles_std"])
    # One smiles_std per inchi_key — take first non-null
    valid = std[std["smiles_std"].notna()]
    return valid.drop_duplicates("inchi_key").set_index("inchi_key")["smiles_std"]


def main():
    parser = argparse.ArgumentParser(description="PAD4-DB Step 03: Replicate Aggregation")
    parser.add_argument("--dry-run", action="store_true",
                        help="Sample 5000 rows and write to _dryrun.parquet")
    args = parser.parse_args()

    for p in [NORM_PARQUET, STD_PARQUET]:
        if not p.exists():
            sys.exit(f"ERROR: required input not found: {p}")

    print(f"Loading {NORM_PARQUET} ...", flush=True)
    df = pd.read_parquet(NORM_PARQUET)
    print(f"  {len(df):,} rows loaded", flush=True)

    if args.dry_run:
        df = df.sample(n=min(5000, len(df)), random_state=42)
        print(f"  DRY-RUN: sampled {len(df):,} rows", flush=True)

    print("Aggregating ...", flush=True)
    result = aggregate(df)
    print(f"  {len(result):,} groups after aggregation", flush=True)

    print(f"Joining smiles_std from {STD_PARQUET} ...", flush=True)
    smiles_lut = load_smiles_lookup()
    result["smiles_std"] = result["inchi_key"].map(smiles_lut)
    print(f"  smiles_std filled: {result['smiles_std'].notna().sum():,} / {len(result):,}", flush=True)

    # ── Final column order ──────────────────────────────────────────────────
    FINAL_COLS = [
        "inchi_key", "smiles_std",
        "source", "aid", "aid_preferred", "layer",
        "endpoint_type", "assay_mechanism_class",
        "n_replicates",
        "pic50_aid", "log_value_mean", "log_value_std",
        "value_nM_min", "value_nM_max",
        "any_qualifier",
        "use_in_potency_model", "hts_flag", "norm_status_worst",
    ]
    result = result[FINAL_COLS]

    # ── Write ───────────────────────────────────────────────────────────────
    out_path = OUT_PATH.parent / (OUT_PATH.stem + "_dryrun.parquet") if args.dry_run else OUT_PATH
    out_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(out_path, index=False)
    print(f"Written → {out_path}", flush=True)

    # ── Summary ─────────────────────────────────────────────────────────────
    print("\n=== SUMMARY ===")
    print(f"Input rows:          {len(df):,}")
    print(f"Output groups:       {len(result):,}")
    print(f"Multi-replicate:     {(result.n_replicates > 1).sum():,}")
    print(f"use_in_potency=True: {result.use_in_potency_model.sum():,}")
    print(f"hts_flag=True:       {result.hts_flag.sum():,}")
    print(f"smiles_std missing:  {result.smiles_std.isna().sum():,}")
    print("\nAssay mechanism class breakdown:")
    print(result.assay_mechanism_class.value_counts().to_string())
    print("\nnorm_status_worst breakdown:")
    print(result.norm_status_worst.value_counts().to_string())
    print("\nendpoint_type breakdown:")
    print(result.endpoint_type.value_counts().to_string())


if __name__ == "__main__":
    main()
