#!/usr/bin/env python3
"""
PAD4-DB v2 — Step 02: Activity Normalization
scripts/02_normalize/02_normalize_activities.py

Run from project root:
    conda activate pad4bench
    python scripts/02_normalize/02_normalize_activities.py [--dry-run] [--source NAME] [--aid INT]

Outputs:
    data/interim/normalized/normalized_activities.parquet
    outputs/tables/02_normalization_report.csv
"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path

import pandas as pd

# =============================================================================
# PATHS
# =============================================================================

BASE_DIR    = Path("data/raw")
INTERIM_DIR = Path("data/interim/normalized")
OUT_DIR     = Path("outputs/tables")
STD_PARQUET = Path("data/interim/standardized/standardized_compounds.parquet")

# ── Snakemake integration ─────────────────────────────────────────────────
if "snakemake" in dir():
    _sm = snakemake  # noqa
    STD_PARQUET = Path(_sm.input.standardized)
    INTERIM_DIR = Path(_sm.output.parquet).parent
    OUT_DIR     = Path(_sm.output.report).parent

BINDINGDB_ENDPOINT_COLS = ["IC50 (nM)", "Ki (nM)", "Kd (nM)", "EC50 (nM)"]

# =============================================================================
# REGISTRY  — identical to 01_standardize_smiles.py
# Dual-layer AIDs 1920046, 2202442: canonical = confirmatory/; skipped in secondary
# =============================================================================

REGISTRY = {
    "HTS": {
        463073: {"subdir": "hts", "layer": "HTS"},
        485272: {"subdir": "hts", "layer": "HTS"},
        488796: {"subdir": "hts", "layer": "HTS"},
    },
    "A": {
        492970:  {"subdir": "pubchem/confirmatory", "layer": "A"},
        320707:  {"subdir": "pubchem/confirmatory", "layer": "A"},
        1804546: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1804627: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1805620: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1806182: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1806183: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1806764: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1919095: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1920200: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1963715: {"subdir": "pubchem/confirmatory", "layer": "A"},
        2202576: {"subdir": "pubchem/confirmatory", "layer": "A"},
        2202577: {"subdir": "pubchem/confirmatory", "layer": "A"},
        2202596: {"subdir": "pubchem/confirmatory", "layer": "A"},
        2202597: {"subdir": "pubchem/confirmatory", "layer": "A"},
        2202717: {"subdir": "pubchem/confirmatory", "layer": "A"},
        2202442: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1330527: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1471656: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1474465: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1474486: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1511938: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1632998: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1651501: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1651502: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1813082: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1813806: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1875531: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1920046: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1973686: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1993481: {"subdir": "pubchem/confirmatory", "layer": "A"},
        2006936: {"subdir": "pubchem/confirmatory", "layer": "A"},
        2034405: {"subdir": "pubchem/confirmatory", "layer": "A"},
        2034406: {"subdir": "pubchem/confirmatory", "layer": "A"},
        2034407: {"subdir": "pubchem/confirmatory", "layer": "A"},
        2071731: {"subdir": "pubchem/confirmatory", "layer": "A"},
        2134413: {"subdir": "pubchem/confirmatory", "layer": "A"},
        2200614: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1422897: {"subdir": "pubchem/confirmatory", "layer": "A"},
        725673:  {"subdir": "pubchem/confirmatory", "layer": "A"},
        725596:  {"subdir": "pubchem/confirmatory", "layer": "A"},
        725597:  {"subdir": "pubchem/confirmatory", "layer": "A"},
        725598:  {"subdir": "pubchem/confirmatory", "layer": "A"},
        1196525: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1069613: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1069619: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1069623: {"subdir": "pubchem/confirmatory", "layer": "A"},
        627371:  {"subdir": "pubchem/confirmatory", "layer": "A"},
        627428:  {"subdir": "pubchem/confirmatory", "layer": "A"},
        627432:  {"subdir": "pubchem/confirmatory", "layer": "A"},
        626724:  {"subdir": "pubchem/confirmatory", "layer": "A"},
        626728:  {"subdir": "pubchem/confirmatory", "layer": "A"},
        1069614: {"subdir": "pubchem/confirmatory", "layer": "A"},
        1069618: {"subdir": "pubchem/confirmatory", "layer": "A"},
        588488:  {"subdir": "pubchem/confirmatory", "layer": "A"},
        588559:  {"subdir": "pubchem/confirmatory", "layer": "A"},
        588560:  {"subdir": "pubchem/confirmatory", "layer": "A"},
    },
    "C": {
        626735:  {"subdir": "pubchem/literature_derived", "layer": "C"},
        626738:  {"subdir": "pubchem/literature_derived", "layer": "C"},
        712876:  {"subdir": "pubchem/literature_derived", "layer": "C"},
        725671:  {"subdir": "pubchem/literature_derived", "layer": "C"},
        1069608: {"subdir": "pubchem/literature_derived", "layer": "C"},
        1196521: {"subdir": "pubchem/literature_derived", "layer": "C"},
        1364668: {"subdir": "pubchem/literature_derived", "layer": "C"},
        1422898: {"subdir": "pubchem/literature_derived", "layer": "C"},
        1422904: {"subdir": "pubchem/literature_derived", "layer": "C"},
        2076402: {"subdir": "pubchem/literature_derived", "layer": "C"},
        651867:  {"subdir": "pubchem/literature_derived", "layer": "C"},
    },
    "D": {
        2041348: {"subdir": "pubchem/secondary", "layer": "D"},
        2041349: {"subdir": "pubchem/secondary", "layer": "D"},
        2053867: {"subdir": "pubchem/secondary", "layer": "D"},
        2053915: {"subdir": "pubchem/secondary", "layer": "D"},
        2053917: {"subdir": "pubchem/secondary", "layer": "D"},
        2193457: {"subdir": "pubchem/secondary", "layer": "D"},
        1625405: {"subdir": "pubchem/secondary", "layer": "D"},
        1806765: {"subdir": "pubchem/secondary", "layer": "D"},
        1069597: {"subdir": "pubchem/secondary", "layer": "D"},
        1069598: {"subdir": "pubchem/secondary", "layer": "D"},
        1069599: {"subdir": "pubchem/secondary", "layer": "D"},
        1069600: {"subdir": "pubchem/secondary", "layer": "D"},
        1069601: {"subdir": "pubchem/secondary", "layer": "D"},
        1069604: {"subdir": "pubchem/secondary", "layer": "D"},
        1069605: {"subdir": "pubchem/secondary", "layer": "D"},
        1069606: {"subdir": "pubchem/secondary", "layer": "D"},
        588487:  {"subdir": "pubchem/secondary", "layer": "D"},
        651627:  {"subdir": "pubchem/secondary", "layer": "D"},
        1346144: {"subdir": "pubchem/secondary", "layer": "D"},
    },
    "E": {
        2039667: {"subdir": "pubchem/secondary", "layer": "E"},
        2053916: {"subdir": "pubchem/secondary", "layer": "E"},
        627442:  {"subdir": "pubchem/secondary", "layer": "E"},
        627443:  {"subdir": "pubchem/secondary", "layer": "E"},
        1422917: {"subdir": "pubchem/secondary", "layer": "E"},
    },
}

# =============================================================================
# UNIT CONVERSION
# =============================================================================

_UNIT_FACTORS = {
    "nm": 1.0, "nanomolar": 1.0, "nmol/l": 1.0,
    "um": 1e3, "µm": 1e3, "micromolar": 1e3, "umol/l": 1e3,
    "mm": 1e6, "millimolar": 1e6, "mmol/l": 1e6,
    "pm": 1e-3, "picomolar": 1e-3, "pmol/l": 1e-3,
    "m": 1e9, "mol/l": 1e9,
}


def convert_to_nM(value, units_raw: str):
    """
    Returns (value_nM float or None, norm_status str).
    """
    if value is None:
        return None, "NO_VALUE"
    units_key = units_raw.strip().lower().replace(" ", "")
    factor = _UNIT_FACTORS.get(units_key)
    if factor is None:
        return None, "UNCONVERTIBLE_UNITS"
    try:
        return float(value) * factor, "OK"
    except (ValueError, TypeError):
        return None, "PARSE_ERROR"


# =============================================================================
# QUALIFIER PARSING
# =============================================================================

_QUALIFIER_RE = re.compile(r"^([<>~=]{1,2})\s*")


def parse_value_qualifier(s: str):
    """
    Parse a raw value string that may have a leading qualifier.
    Returns (qualifier str, value float or None, norm_status str).
    """
    s = s.strip()
    if not s:
        return "", None, "NO_VALUE"

    m = _QUALIFIER_RE.match(s)
    qualifier = m.group(1) if m else ""
    # Normalise ">=" → ">"  "<=" → "<"
    qualifier = qualifier.replace("=", "").replace(">=", ">").replace("<=", "<")
    # Re-clean after replace
    qualifier = qualifier.strip()

    num_str = s[m.end():].strip() if m else s
    if not num_str:
        return qualifier, None, "QUALIFIER_ONLY"

    try:
        return qualifier, float(num_str), "OK"
    except ValueError:
        return qualifier, None, "PARSE_ERROR"


# =============================================================================
# ENDPOINT TYPE NORMALISATION
# =============================================================================

_ENDPOINT_MAP = {
    "ic50": "IC50", "ic 50": "IC50", "pic50": "IC50",
    "ki": "Ki", "kd": "Kd",
    "ec50": "EC50",
    "ratio": "Kinact_Ki", "kinact/ki": "Kinact_Ki", "kinact_ki": "Kinact_Ki",
    "inhibition": "Pct_inh", "% inhibition": "Pct_inh",
    "pct_inh": "Pct_inh", "activity": "Pct_inh",
    "potency": "IC50",
}


def normalize_endpoint_type(raw: str) -> str:
    if not raw:
        return ""
    return _ENDPOINT_MAP.get(raw.strip().lower(), raw.strip())


# =============================================================================
# pIC50
# =============================================================================

def compute_pic50(value_nM):
    if value_nM is None or value_nM <= 0:
        return None
    return -math.log10(value_nM * 1e-9)


# =============================================================================
# PUBCHEM CSV PARSING  (4-line metadata header)
# =============================================================================

def parse_pubchem_csv(filepath: Path, max_rows=None):
    """Returns (columns list, data_rows list of dicts)."""
    if not filepath.exists():
        print(f"  ERROR: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)

    with open(filepath, newline="", encoding="utf-8", errors="replace") as f:
        raw = f.read()

    lines = raw.splitlines()
    header_idx = next(
        (i for i, l in enumerate(lines[:10]) if "PUBCHEM_SID" in l), 0
    )
    reader = csv.DictReader(lines[header_idx:])
    cols = list(reader.fieldnames or [])
    rows = []
    for row in reader:
        sid = row.get("PUBCHEM_SID", "").strip()
        if not sid.lstrip("-").isdigit():
            continue
        rows.append(row)
        if max_rows and len(rows) >= max_rows:
            break
    return cols, rows


# =============================================================================
# ACTIVITY EXTRACTION — PubChem
# =============================================================================

def _extract_activity_from_row(row: dict, col_set: set):
    """
    Determine activity value from a PubChem row using column priority.

    Returns (endpoint_type, endpoint_source_col, value_raw, qualifier,
             value_nM, units_raw, norm_status, pic50_direct).
    pic50_direct is only set for pIC50_min/max rows.
    """
    # Strategy 1: Standard Value + Standard Units + Standard Type
    if "Standard Value" in col_set and "Standard Units" in col_set:
        ep = normalize_endpoint_type(row.get("Standard Type", ""))
        qual_col = (row.get("Standard Relation", "") or "").strip().strip("'")
        val_raw  = (row.get("Standard Value", "") or "").strip()
        units    = (row.get("Standard Units", "") or "").strip()

        if not val_raw:
            return ep, "Standard Value", "", qual_col or "", None, units, "NO_VALUE", None

        q_parsed, val_f, st = parse_value_qualifier(val_raw)
        qualifier = qual_col if qual_col and qual_col not in ("=",) else q_parsed
        if st == "QUALIFIER_ONLY":
            return ep, "Standard Value", val_raw, qualifier, None, units, "QUALIFIER_ONLY", None
        if st == "PARSE_ERROR":
            return ep, "Standard Value", val_raw, qualifier, None, units, "PARSE_ERROR", None
        # Pct_inh (% units) → value has no nM equivalent; status is OK
        if ep == "Pct_inh" or units.strip() == "%":
            return ep, "Standard Value", val_raw, qualifier, None, units, "OK", None
        value_nM, norm_st = convert_to_nM(val_f, units)
        return ep, "Standard Value", val_raw, qualifier, value_nM, units, norm_st, None

    # Strategy 2: IC50 column (nM) — AID 1804546, 1919095, 1920200
    if "IC50" in col_set:
        ep = normalize_endpoint_type(row.get("Standard Type", "IC50") or "IC50")
        qual_col = (row.get("IC50 Qualifier", row.get("Standard Relation", "")) or "").strip().strip("'")
        val_raw  = (row.get("IC50", "") or "").strip()
        if not val_raw:
            return ep, "IC50", "", qual_col or "", None, "nM", "NO_VALUE", None
        q_parsed, val_f, st = parse_value_qualifier(val_raw)
        qualifier = qual_col or q_parsed
        if st not in ("OK", "QUALIFIER_ONLY"):
            return ep, "IC50", val_raw, qualifier, None, "nM", st, None
        value_nM, norm_st = convert_to_nM(val_f, "nM")
        return ep, "IC50", val_raw, qualifier, value_nM, "nM", norm_st, None

    # Strategy 3: PubChem Standard Value (µM) + Standard Type
    if "PubChem Standard Value" in col_set:
        ep = normalize_endpoint_type(row.get("Standard Type", "") or "")
        qual_col = (row.get("Standard Relation", "") or "").strip().strip("'")
        val_raw  = (row.get("PubChem Standard Value", "") or "").strip()
        if not val_raw:
            return ep, "PubChem Standard Value", "", qual_col or "", None, "uM", "NO_VALUE", None
        q_parsed, val_f, st = parse_value_qualifier(val_raw)
        qualifier = qual_col or q_parsed
        if st not in ("OK", "QUALIFIER_ONLY"):
            return ep, "PubChem Standard Value", val_raw, qualifier, None, "uM", st, None
        value_nM, norm_st = convert_to_nM(val_f, "uM")
        return ep, "PubChem Standard Value", val_raw, qualifier, value_nM, "uM", norm_st, None

    # Strategy 4: Average IC50 (µM) — AID 492970
    if "Average IC50" in col_set:
        qual_col = (row.get("Qualifier", "") or "").strip().strip("'")
        val_raw  = (row.get("Average IC50", "") or "").strip()
        if not val_raw:
            return "IC50", "Average IC50", "", qual_col or "", None, "uM", "NO_VALUE", None
        q_parsed, val_f, st = parse_value_qualifier(val_raw)
        qualifier = qual_col or q_parsed
        if st not in ("OK", "QUALIFIER_ONLY"):
            return "IC50", "Average IC50", val_raw, qualifier, None, "uM", st, None
        value_nM, norm_st = convert_to_nM(val_f, "uM")
        return "IC50", "Average IC50", val_raw, qualifier, value_nM, "uM", norm_st, None

    # Strategy 5: pIC50_min / pIC50_max — AID 1346144
    if "pIC50_min" in col_set or "pIC50_max" in col_set:
        val_raw = ((row.get("pIC50_min") or row.get("pIC50_max")) or "").strip()
        src_col = "pIC50_min" if "pIC50_min" in col_set else "pIC50_max"
        if not val_raw:
            return "IC50", src_col, "", "", None, "nM", "NO_VALUE", None
        try:
            pic50 = float(val_raw)
            value_nM = 10.0 ** (9.0 - pic50)
            return "IC50", src_col, val_raw, "=", value_nM, "nM", "OK", pic50
        except ValueError:
            return "IC50", src_col, val_raw, "", None, "nM", "PARSE_ERROR", None

    # Strategy 6: Inhibition columns (HTS, screening panels)
    for col in ("Inhibition", "Inhibition at 6.36 uM", "Inhibition at 10 uM",
                "Inhibition at 100 uM"):
        if col in col_set:
            val_raw = (row.get(col, "") or "").strip()
            if not val_raw:
                return "Pct_inh", col, "", "", None, "%", "NO_VALUE", None
            q, val_f, st = parse_value_qualifier(val_raw)
            return "Pct_inh", col, val_raw, q, None, "%", "OK", None

    # Strategy 7: Qualitative only — no numeric value (Activity Comment, mass shift)
    ep = normalize_endpoint_type(row.get("Standard Type", "") or "")
    outcome = (row.get("PUBCHEM_ACTIVITY_OUTCOME", "") or "").strip()
    return ep or "Unknown", "PUBCHEM_ACTIVITY_OUTCOME", outcome, "", None, "", "NO_VALUE", None


def _cid_from_row(row: dict) -> str:
    cid = (row.get("PUBCHEM_CID", "") or "").strip()
    if cid and cid not in ("0", ""):
        return cid
    return (row.get("PUBCHEM_SID", "") or "").strip()


def extract_pubchem_activities(aid: int, meta: dict, source: str,
                               std_lookup: dict, dry_run: bool) -> list:
    subdir = meta["subdir"]
    layer  = meta["layer"]
    fp = BASE_DIR / subdir / f"AID_{aid}_datatable_all.csv"

    cols, rows = parse_pubchem_csv(fp, max_rows=100 if dry_run else None)
    col_set = set(cols)
    out = []

    for row in rows:
        cid_or_id = _cid_from_row(row)
        lookup_key = (source, str(aid), cid_or_id)
        std_info = std_lookup.get(lookup_key, {})
        inchi_key   = std_info.get("inchi_key")
        aid_pref    = std_info.get("aid_preferred", str(aid))

        (ep_type, ep_src_col, val_raw, qualifier,
         value_nM, units_raw, norm_status, pic50_direct) = _extract_activity_from_row(row, col_set)

        if ep_type == "Pct_inh":
            pic50 = None
        elif pic50_direct is not None:
            pic50 = pic50_direct
        else:
            pic50 = compute_pic50(value_nM) if ep_type == "IC50" else None

        use_in_model = (
            ep_type == "IC50"
            and value_nM is not None
            and value_nM > 0
            and layer not in ("HTS", "D", "E")
        )

        out.append({
            "inchi_key":          inchi_key,
            "source":             source,
            "aid":                str(aid),
            "aid_preferred":      aid_pref,
            "layer":              layer,
            "endpoint_type":      ep_type,
            "endpoint_source_col": ep_src_col,
            "value_raw":          val_raw,
            "qualifier":          qualifier,
            "value_nM":           value_nM,
            "pIC50":              pic50,
            "units_raw":          units_raw,
            "norm_status":        norm_status,
            "use_in_potency_model": use_in_model,
            "pchembl_value":      None,
            "pchembl_mismatch":   False,
        })
    return out


# =============================================================================
# STANDARDIZED PARQUET LOOKUP
# =============================================================================

def load_standardized() -> dict:
    """
    Returns dict: (source, aid_str, cid_or_id_str) → {inchi_key, aid_preferred}.
    """
    if not STD_PARQUET.exists():
        print(f"ERROR: Standardized parquet not found: {STD_PARQUET}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_parquet(STD_PARQUET, columns=["source", "aid", "cid_or_id",
                                                "inchi_key", "aid_preferred"])
    lookup = {}
    for row in df.itertuples(index=False):
        key = (row.source, str(row.aid), str(row.cid_or_id))
        lookup[key] = {"inchi_key": row.inchi_key, "aid_preferred": row.aid_preferred}
    return lookup


# =============================================================================
# CHEMBL EXTRACTION
# =============================================================================

def extract_chembl_activities(std_lookup: dict, dry_run: bool) -> list:
    chembl_dir = BASE_DIR / "chembl"
    files = list(chembl_dir.glob("*.csv"))
    if not files:
        print("ERROR: No CSV in data/raw/chembl/", file=sys.stderr)
        sys.exit(1)

    fp = files[0]
    print(f"  [chembl] {fp.name} ...", end=" ", flush=True)
    out = []

    with open(fp, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter=";")
        for i, row in enumerate(reader):
            if dry_run and i >= 100:
                break

            cid_or_id = (row.get("Molecule ChEMBL ID", "") or "").strip()
            lookup_key = ("chembl", "CHEMBL6111", cid_or_id)
            std_info = std_lookup.get(lookup_key, {})
            inchi_key  = std_info.get("inchi_key")
            aid_pref   = std_info.get("aid_preferred", "CHEMBL6111")

            ep_raw   = (row.get("Standard Type", "") or "").strip()
            ep_type  = normalize_endpoint_type(ep_raw)
            qual_raw = (row.get("Standard Relation", "") or "").strip().strip("'")
            val_raw  = (row.get("Standard Value", "") or "").strip()
            units    = (row.get("Standard Units", "") or "").strip()

            pchembl_str = (row.get("pChEMBL Value", "") or "").strip()
            try:
                pchembl_value = float(pchembl_str) if pchembl_str else None
            except ValueError:
                pchembl_value = None

            if not val_raw:
                norm_status, value_nM, qualifier = "NO_VALUE", None, qual_raw
            else:
                q_parsed, val_f, st = parse_value_qualifier(val_raw)
                qualifier = qual_raw if qual_raw and qual_raw not in ("=",) else q_parsed
                if st == "QUALIFIER_ONLY":
                    norm_status, value_nM = "QUALIFIER_ONLY", None
                elif st == "PARSE_ERROR":
                    norm_status, value_nM = "PARSE_ERROR", None
                elif ep_type == "Pct_inh" or units.strip() == "%":
                    norm_status, value_nM = "OK", None
                else:
                    value_nM, norm_status = convert_to_nM(val_f, units)

            pic50 = compute_pic50(value_nM) if ep_type == "IC50" else None

            if pchembl_value is not None and pic50 is not None:
                pchembl_mismatch = abs(pic50 - pchembl_value) > 0.1
            else:
                pchembl_mismatch = False

            use_in_model = (
                ep_type == "IC50"
                and value_nM is not None
                and value_nM > 0
            )

            out.append({
                "inchi_key":           inchi_key,
                "source":              "chembl",
                "aid":                 "CHEMBL6111",
                "aid_preferred":       aid_pref,
                "layer":               "chembl",
                "endpoint_type":       ep_type,
                "endpoint_source_col": "Standard Value",
                "value_raw":           val_raw,
                "qualifier":           qualifier,
                "value_nM":            value_nM,
                "pIC50":               pic50,
                "units_raw":           units,
                "norm_status":         norm_status,
                "use_in_potency_model": use_in_model,
                "pchembl_value":       pchembl_value,
                "pchembl_mismatch":    pchembl_mismatch,
            })

    print(f"{len(out)} rows, OK")
    return out


# =============================================================================
# BINDINGDB EXTRACTION  (explode: up to 4 rows per input row)
# =============================================================================

_BDB_ENDPOINT_TYPE = {
    "IC50 (nM)": "IC50",
    "Ki (nM)":   "Ki",
    "Kd (nM)":   "Kd",
    "EC50 (nM)": "EC50",
}


def extract_bindingdb_activities(std_lookup: dict, dry_run: bool) -> list:
    bdb_dir = BASE_DIR / "bindingdb"
    files = list(bdb_dir.glob("*.tsv")) + list(bdb_dir.glob("*.csv"))
    if not files:
        print("ERROR: No TSV/CSV in data/raw/bindingdb/", file=sys.stderr)
        sys.exit(1)

    fp = files[0]
    sep = "\t" if fp.suffix == ".tsv" else ","
    print(f"  [bindingdb] {fp.name} ...", end=" ", flush=True)
    out = []

    with open(fp, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter=sep)
        for i, row in enumerate(reader):
            if dry_run and i >= 100:
                break

            cid_or_id = (row.get("BindingDB MonomerID", "") or "").strip()
            lookup_key = ("bindingdb", "Q9UM07", cid_or_id)
            std_info = std_lookup.get(lookup_key, {})
            inchi_key  = std_info.get("inchi_key")
            aid_pref   = std_info.get("aid_preferred", "Q9UM07")

            for ep_col in BINDINGDB_ENDPOINT_COLS:
                raw_cell = (row.get(ep_col, "") or "").strip()
                # BindingDB uses whitespace-only cells for missing values
                if not raw_cell or not raw_cell.strip():
                    continue

                ep_type  = _BDB_ENDPOINT_TYPE[ep_col]
                q_parsed, val_f, st = parse_value_qualifier(raw_cell)

                if st == "QUALIFIER_ONLY":
                    value_nM, norm_status = None, "QUALIFIER_ONLY"
                elif st == "PARSE_ERROR":
                    value_nM, norm_status = None, "PARSE_ERROR"
                else:
                    # BindingDB values already in nM
                    value_nM, norm_status = convert_to_nM(val_f, "nM")

                pic50 = compute_pic50(value_nM) if ep_type == "IC50" else None

                use_in_model = (
                    ep_type == "IC50"
                    and value_nM is not None
                    and value_nM > 0
                )

                out.append({
                    "inchi_key":           inchi_key,
                    "source":              "bindingdb",
                    "aid":                 "Q9UM07",
                    "aid_preferred":       aid_pref,
                    "layer":               "bindingdb",
                    "endpoint_type":       ep_type,
                    "endpoint_source_col": ep_col,
                    "value_raw":           raw_cell,
                    "qualifier":           q_parsed,
                    "value_nM":            value_nM,
                    "pIC50":               pic50,
                    "units_raw":           "nM",
                    "norm_status":         norm_status,
                    "use_in_potency_model": use_in_model,
                    "pchembl_value":       None,
                    "pchembl_mismatch":    False,
                })

    print(f"{len(out)} exploded rows from {i+1} input rows, OK")
    return out


# =============================================================================
# BUILD NORMALIZED DATAFRAME
# =============================================================================

_SCHEMA_COLS = [
    "inchi_key", "source", "aid", "aid_preferred", "layer",
    "endpoint_type", "endpoint_source_col", "value_raw", "qualifier",
    "value_nM", "pIC50", "units_raw", "norm_status", "use_in_potency_model",
    "pchembl_value", "pchembl_mismatch",
]


def build_normalized_df(rows: list) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=_SCHEMA_COLS)
    df["value_nM"]       = pd.to_numeric(df["value_nM"],       errors="coerce")
    df["pIC50"]          = pd.to_numeric(df["pIC50"],          errors="coerce")
    df["pchembl_value"]  = pd.to_numeric(df["pchembl_value"],  errors="coerce")
    df["qualifier"]      = df["qualifier"].fillna("")
    df["value_raw"]      = df["value_raw"].fillna("")
    df["units_raw"]      = df["units_raw"].fillna("")
    return df


# =============================================================================
# QC REPORT
# =============================================================================

def build_report(df: pd.DataFrame) -> pd.DataFrame:
    records = []
    for (source, aid, layer), grp in df.groupby(["source", "aid", "layer"], sort=False):
        # Count input rows (before BindingDB explosion they're all present)
        n_exploded = len(grp)
        n_ok            = (grp["norm_status"] == "OK").sum()
        n_no_value      = (grp["norm_status"] == "NO_VALUE").sum()
        n_unconvertible = (grp["norm_status"] == "UNCONVERTIBLE_UNITS").sum()
        n_parse_error   = (grp["norm_status"].isin(["PARSE_ERROR", "QUALIFIER_ONLY"])).sum()
        pct_ok          = round(100 * n_ok / n_exploded, 1) if n_exploded else 0.0
        n_use_model     = grp["use_in_potency_model"].sum()

        ok_pic50 = grp[grp["pIC50"].notna()]["pIC50"]
        records.append({
            "source":          source,
            "aid":             aid,
            "layer":           layer,
            "n_input_rows":    n_exploded,
            "n_exploded_rows": n_exploded,
            "n_ok":            n_ok,
            "n_no_value":      n_no_value,
            "n_unconvertible": n_unconvertible,
            "n_parse_error":   n_parse_error,
            "pct_ok":          pct_ok,
            "n_use_in_model":  n_use_model,
            "median_pic50":    round(ok_pic50.median(), 3) if len(ok_pic50) else None,
            "min_pic50":       round(ok_pic50.min(), 3)    if len(ok_pic50) else None,
            "max_pic50":       round(ok_pic50.max(), 3)    if len(ok_pic50) else None,
        })
    return pd.DataFrame(records)


# =============================================================================
# SUMMARY PRINT
# =============================================================================

def print_summary(df: pd.DataFrame, report_df: pd.DataFrame):
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)

    total           = len(df)
    n_ok            = (df["norm_status"] == "OK").sum()
    n_no_value      = (df["norm_status"] == "NO_VALUE").sum()
    n_unconvertible = (df["norm_status"] == "UNCONVERTIBLE_UNITS").sum()
    n_parse_error   = (df["norm_status"].isin(["PARSE_ERROR", "QUALIFIER_ONLY"])).sum()
    n_use_model     = df["use_in_potency_model"].sum()

    pic50_ok = df[df["pIC50"].notna()]["pIC50"]
    chembl_mismatch = df["pchembl_mismatch"].sum()

    print(f"  Total output rows      : {total:>8,}")
    print(f"  OK                     : {n_ok:>8,}  ({100*n_ok/total:.1f}%)")
    print(f"  NO_VALUE               : {n_no_value:>8,}  ({100*n_no_value/total:.1f}%)")
    print(f"  UNCONVERTIBLE_UNITS    : {n_unconvertible:>8,}  ({100*n_unconvertible/total:.1f}%)")
    print(f"  PARSE_ERROR/QUALIFIER  : {n_parse_error:>8,}  ({100*n_parse_error/total:.1f}%)")
    print(f"  use_in_potency_model   : {n_use_model:>8,}")
    if len(pic50_ok):
        print(f"  pIC50 range            : {pic50_ok.min():.2f} – {pic50_ok.max():.2f}")
    print(f"  ChEMBL pchembl_mismatch: {chembl_mismatch:>8,}")

    print(f"\n  Per-source pct_ok (non-HTS sources flagged if < 85%):")
    for _, r in report_df.iterrows():
        flag = "  ⚠ " if r["layer"] != "HTS" and r["pct_ok"] < 85 else "    "
        print(f"  {flag}{r['source']:30} AID {str(r['aid']):>10}  "
              f"{r['n_exploded_rows']:>7} rows  {r['pct_ok']:>5.1f}% OK  "
              f"model={r['n_use_in_model']}")

    warn = report_df[(report_df["pct_ok"] < 85) & (report_df["layer"] != "HTS")]
    if not warn.empty:
        print(f"\n  ⚠  {len(warn)} non-HTS source(s) below 85% OK:")
        for _, r in warn.iterrows():
            print(f"    {r['source']} AID {r['aid']}: {r['pct_ok']}%  "
                  f"(no_value={r['n_no_value']} unconvertible={r['n_unconvertible']} "
                  f"parse={r['n_parse_error']})")
    else:
        print(f"\n  ✅ All non-HTS sources ≥ 85% OK")


# =============================================================================
# ARGS
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="PAD4-DB v2 activity normalization")
    p.add_argument("--dry-run", action="store_true",
                   help="Process first 100 rows per source; suffix outputs with _dryrun")
    p.add_argument("--source", default=None,
                   help="Filter: pubchem, pubchem_hts, pubchem_confirmatory, "
                        "pubchem_literature_derived, pubchem_secondary, chembl, bindingdb")
    p.add_argument("--aid", type=int, default=None,
                   help="Filter to a single AID integer (for debugging)")
    return p.parse_args()


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()
    dry_run       = args.dry_run
    source_filter = args.source
    aid_filter    = args.aid

    suffix      = "_dryrun" if dry_run else ""
    parquet_out = INTERIM_DIR / f"normalized_activities{suffix}.parquet"
    report_out  = OUT_DIR    / f"02_normalization_report{suffix}.csv"

    INTERIM_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("PAD4-DB v2 — Step 02: Activity Normalization")
    if dry_run:
        print("  [DRY-RUN: first 100 rows per source]")
    if source_filter:
        print(f"  [SOURCE FILTER: {source_filter}]")
    if aid_filter:
        print(f"  [AID FILTER: {aid_filter}]")
    print("=" * 72)

    print("\nLoading standardized compounds lookup ...", end=" ", flush=True)
    std_lookup = load_standardized()
    print(f"{len(std_lookup):,} entries")

    all_rows = []

    print("\n── PubChem sources ───────────────────────────────────────────────────")
    for category, aids in REGISTRY.items():
        for aid, meta in aids.items():
            if aid_filter and aid != aid_filter:
                continue
            subdir = meta["subdir"]
            layer  = meta["layer"]
            if subdir == "hts":
                source = "pubchem_hts"
            else:
                part = subdir.replace("pubchem/", "").replace("/", "_")
                source = f"pubchem_{part}"

            if source_filter and source_filter not in (source, "pubchem"):
                continue

            print(f"  [{source}] AID {aid} ({layer}) ...", end=" ", flush=True)
            rows = extract_pubchem_activities(aid, meta, source, std_lookup, dry_run)
            all_rows.extend(rows)
            n_ok = sum(1 for r in rows if r["norm_status"] == "OK")
            print(f"{len(rows)} rows, {n_ok} OK ({100*n_ok/len(rows):.0f}%)" if rows else "0 rows WARNING")

    print("\n── ChEMBL ────────────────────────────────────────────────────────────")
    if not source_filter or source_filter == "chembl":
        all_rows.extend(extract_chembl_activities(std_lookup, dry_run))

    print("\n── BindingDB ─────────────────────────────────────────────────────────")
    if not source_filter or source_filter == "bindingdb":
        all_rows.extend(extract_bindingdb_activities(std_lookup, dry_run))

    if not all_rows:
        print("\nNo rows processed. Check --source / --aid filters.")
        sys.exit(1)

    print(f"\n── Building DataFrame ({len(all_rows):,} rows) ───────────────────────────────")
    df = build_normalized_df(all_rows)

    print(f"── Writing parquet: {parquet_out}")
    df.to_parquet(parquet_out, index=False)
    print(f"  {len(df):,} rows written")

    report_df = build_report(df)
    report_df.to_csv(report_out, index=False)
    print(f"  Report: {report_out}")

    print_summary(df, report_df)
    print(f"\n  Output : {parquet_out}")
    print(f"  Report : {report_out}")


if __name__ == "__main__":
    main()
