#!/usr/bin/env python3
"""
PAD4-DB v2 — Step 01: SMILES Standardization
scripts/01_standardize/01_standardize_smiles.py

Run from project root:
    conda activate pad4bench
    python scripts/01_standardize/01_standardize_smiles.py [--dry-run] [--source NAME]

Outputs:
    data/interim/standardized/standardized_compounds.parquet
    outputs/tables/01_standardization_report.csv
"""

import argparse
import csv
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

# =============================================================================
# PATHS
# =============================================================================

BASE_DIR    = Path("data/raw")
INTERIM_DIR = Path("data/interim/standardized")
OUT_DIR     = Path("outputs/tables")

# =============================================================================
# REGISTRY  — mirrors 00_check_raw_inventory.py, dual-layer AIDs load from
# confirmatory/ only (1920046, 2202442 excluded from secondary here)
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
        2202442: {"subdir": "pubchem/confirmatory", "layer": "A"},  # dual-layer: confirmatory is canonical
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
        1920046: {"subdir": "pubchem/confirmatory", "layer": "A"},  # dual-layer: confirmatory is canonical
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
        # 1920046 and 2202442 intentionally omitted here — loaded from confirmatory/ (Layer A)
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

SMILES_COLS = [
    "PUBCHEM_EXT_DATASOURCE_SMILES",
    "PUBCHEM_OPENEYE_ISO_SMILES",
    "PUBCHEM_CANONICAL_SMILES",
    "PUBCHEM_ISOMERIC_SMILES",
]

# AIDs where overlapping SIDs should prefer AID 2202577
RFMS_AID_PREFERRED = 2202577
RFMS_AIDS = {2202576, 2202577}


# =============================================================================
# CHEMISTRY
# =============================================================================

_fragment_chooser = rdMolStandardize.LargestFragmentChooser()
_uncharger = rdMolStandardize.Uncharger()


def _strip_extended_smiles(smiles: str) -> str:
    """Strip Daylight extended SMILES annotations (e.g. ' |r,THB:...|') before parsing."""
    idx = smiles.find(" |")
    if idx != -1:
        smiles = smiles[:idx]
    return smiles.strip()


def standardize_mol(smiles: str) -> tuple:
    """
    Returns (smiles_std, inchi_key, n_heavy_atoms, mol_weight, std_status).
    All values except std_status are None on failure.
    """
    if not smiles or not smiles.strip():
        return None, None, None, None, "NO_SMILES"

    mol = Chem.MolFromSmiles(_strip_extended_smiles(smiles))
    if mol is None:
        return None, None, None, None, "PARSE_FAIL"

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None, None, None, None, "SANITIZE_FAIL"

    try:
        mol = _fragment_chooser.choose(mol)
        mol = _uncharger.uncharge(mol)
        smiles_std = Chem.MolToSmiles(mol, isomericSmiles=True)
        inchi_key  = Chem.inchi.MolToInchiKey(mol)
        n_heavy    = mol.GetNumHeavyAtoms()
        mw         = round(Descriptors.MolWt(mol), 4)
    except Exception:
        return None, None, None, None, "SANITIZE_FAIL"

    return smiles_std, inchi_key, n_heavy, mw, "OK"


# =============================================================================
# PUBCHEM PARSING
# =============================================================================

def parse_pubchem_csv(filepath: Path, max_rows: int | None = None) -> list[dict]:
    """
    Parse a PubChem AID CSV with 4-line metadata header.
    Returns list of raw row dicts.
    """
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
    rows = []
    for row in reader:
        sid = row.get("PUBCHEM_SID", "").strip()
        if not sid.lstrip("-").isdigit():
            continue  # skip metadata description rows
        rows.append(row)
        if max_rows and len(rows) >= max_rows:
            break

    return rows


def extract_smiles_pubchem(row: dict) -> str:
    for col in SMILES_COLS:
        val = row.get(col, "").strip()
        if val:
            return val
    return ""


def extract_cid(row: dict) -> str:
    cid = row.get("PUBCHEM_CID", "").strip()
    if cid and cid != "0":
        return cid
    return row.get("PUBCHEM_SID", "").strip()


# =============================================================================
# SOURCE ITERATORS
# =============================================================================

def iter_pubchem_sources(dry_run: bool, source_filter: str | None) -> list[dict]:
    """
    Yields processed rows from all PubChem AIDs.
    Handles RFMS overlap tagging for AIDs 2202576/77.
    """
    max_rows = 100 if dry_run else None
    rows_out = []

    # Collect RFMS SIDs for both AIDs first (needed for overlap tagging)
    rfms_sids: dict[int, set] = {}
    if source_filter is None or source_filter in ("pubchem", "pubchem_confirmatory"):
        for aid in RFMS_AIDS:
            meta = REGISTRY["A"][aid]
            fp = BASE_DIR / meta["subdir"] / f"AID_{aid}_datatable_all.csv"
            if fp.exists():
                raw_rows = parse_pubchem_csv(fp)
                rfms_sids[aid] = {r.get("PUBCHEM_SID", "").strip() for r in raw_rows}

    overlap_sids: set = set()
    if len(rfms_sids) == 2:
        overlap_sids = rfms_sids[2202576] & rfms_sids[2202577]

    for category, aids in REGISTRY.items():
        for aid, meta in aids.items():
            subdir   = meta["subdir"]
            layer    = meta["layer"]
            source   = f"pubchem_{subdir.replace('pubchem/', '').replace('/', '_')}"
            if subdir == "hts":
                source = "pubchem_hts"

            if source_filter and source_filter not in (source, "pubchem"):
                continue

            fp = BASE_DIR / subdir / f"AID_{aid}_datatable_all.csv"
            print(f"  [{source}] AID {aid} ({layer}) ...", end=" ", flush=True)

            raw_rows = parse_pubchem_csv(fp, max_rows=max_rows)
            count = 0
            for row in raw_rows:
                sid        = row.get("PUBCHEM_SID", "").strip()
                smiles_raw = extract_smiles_pubchem(row)
                cid_or_id  = extract_cid(row)

                # RFMS overlap tagging
                if aid in RFMS_AIDS:
                    aid_preferred = str(RFMS_AID_PREFERRED) if sid in overlap_sids else str(aid)
                else:
                    aid_preferred = str(aid)

                smiles_std, inchi_key, n_heavy, mw, status = standardize_mol(smiles_raw)
                rows_out.append({
                    "source":       source,
                    "aid":          str(aid),
                    "layer":        layer,
                    "cid_or_id":    cid_or_id,
                    "smiles_raw":   smiles_raw,
                    "smiles_std":   smiles_std,
                    "inchi_key":    inchi_key,
                    "std_status":   status,
                    "n_heavy_atoms": n_heavy,
                    "mol_weight":   mw,
                    "aid_preferred": aid_preferred,
                })
                count += 1

            print(f"{count} rows, OK")

    return rows_out


def iter_chembl_source(dry_run: bool, source_filter: str | None) -> list[dict]:
    if source_filter and source_filter != "chembl":
        return []

    chembl_dir = BASE_DIR / "chembl"
    files = list(chembl_dir.glob("*.csv"))
    if not files:
        print("  ERROR: No CSV found in data/raw/chembl/", file=sys.stderr)
        sys.exit(1)

    fp = files[0]
    print(f"  [chembl] {fp.name} ...", end=" ", flush=True)

    rows_out = []
    with open(fp, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter=";")
        for i, row in enumerate(reader):
            if dry_run and i >= 100:
                break
            smiles_raw = row.get("Smiles", "").strip()
            cid_or_id  = row.get("Molecule ChEMBL ID", "").strip()
            smiles_std, inchi_key, n_heavy, mw, status = standardize_mol(smiles_raw)
            rows_out.append({
                "source":        "chembl",
                "aid":           "CHEMBL6111",
                "layer":         "chembl",
                "cid_or_id":     cid_or_id,
                "smiles_raw":    smiles_raw,
                "smiles_std":    smiles_std,
                "inchi_key":     inchi_key,
                "std_status":    status,
                "n_heavy_atoms": n_heavy,
                "mol_weight":    mw,
                "aid_preferred": "CHEMBL6111",
            })

    print(f"{len(rows_out)} rows, OK")
    return rows_out


def iter_bindingdb_source(dry_run: bool, source_filter: str | None) -> list[dict]:
    if source_filter and source_filter != "bindingdb":
        return []

    bdb_dir = BASE_DIR / "bindingdb"
    files = list(bdb_dir.glob("*.tsv")) + list(bdb_dir.glob("*.csv"))
    if not files:
        print("  ERROR: No TSV/CSV found in data/raw/bindingdb/", file=sys.stderr)
        sys.exit(1)

    fp = files[0]
    sep = "\t" if fp.suffix == ".tsv" else ","
    print(f"  [bindingdb] {fp.name} ...", end=" ", flush=True)

    rows_out = []
    with open(fp, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter=sep)
        for i, row in enumerate(reader):
            if dry_run and i >= 100:
                break
            smiles_raw = row.get("Ligand SMILES", "").strip()
            cid_or_id  = row.get("BindingDB MonomerID", "").strip()
            smiles_std, inchi_key, n_heavy, mw, status = standardize_mol(smiles_raw)
            rows_out.append({
                "source":        "bindingdb",
                "aid":           "Q9UM07",
                "layer":         "bindingdb",
                "cid_or_id":     cid_or_id,
                "smiles_raw":    smiles_raw,
                "smiles_std":    smiles_std,
                "inchi_key":     inchi_key,
                "std_status":    status,
                "n_heavy_atoms": n_heavy,
                "mol_weight":    mw,
                "aid_preferred": "Q9UM07",
            })

    print(f"{len(rows_out)} rows, OK")
    return rows_out


# =============================================================================
# QC REPORT
# =============================================================================

def build_report(df: pd.DataFrame) -> pd.DataFrame:
    records = []
    for (source, aid, layer), grp in df.groupby(["source", "aid", "layer"], sort=False):
        n_input         = len(grp)
        n_ok            = (grp["std_status"] == "OK").sum()
        n_parse_fail    = (grp["std_status"] == "PARSE_FAIL").sum()
        n_sanitize_fail = (grp["std_status"] == "SANITIZE_FAIL").sum()
        n_no_smiles     = (grp["std_status"] == "NO_SMILES").sum()
        pct_ok          = round(100 * n_ok / n_input, 1) if n_input else 0.0
        ok_rows         = grp[grp["std_status"] == "OK"]["n_heavy_atoms"].dropna()
        median_heavy    = round(ok_rows.median(), 1) if len(ok_rows) else None
        records.append({
            "source":           source,
            "aid":              aid,
            "layer":            layer,
            "n_input":          n_input,
            "n_ok":             n_ok,
            "n_parse_fail":     n_parse_fail,
            "n_sanitize_fail":  n_sanitize_fail,
            "n_no_smiles":      n_no_smiles,
            "pct_ok":           pct_ok,
            "median_heavy_atoms": median_heavy,
        })
    return pd.DataFrame(records)


# =============================================================================
# MAIN
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="PAD4-DB v2 SMILES standardization")
    p.add_argument("--dry-run", action="store_true",
                   help="Process first 100 rows per source; suffix output files with _dryrun")
    p.add_argument("--source", default=None,
                   help="Filter to one source: pubchem, pubchem_hts, pubchem_confirmatory, "
                        "pubchem_literature_derived, pubchem_secondary, chembl, bindingdb")
    return p.parse_args()


def main():
    args = parse_args()
    dry_run       = args.dry_run
    source_filter = args.source

    suffix = "_dryrun" if dry_run else ""

    parquet_out = INTERIM_DIR / f"standardized_compounds{suffix}.parquet"
    report_out  = OUT_DIR / f"01_standardization_report{suffix}.csv"

    INTERIM_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("PAD4-DB v2 — Step 01: SMILES Standardization")
    if dry_run:
        print("  [DRY-RUN: first 100 rows per source]")
    if source_filter:
        print(f"  [SOURCE FILTER: {source_filter}]")
    print("=" * 72)

    all_rows = []

    print("\n── PubChem sources ───────────────────────────────────────────────────")
    all_rows.extend(iter_pubchem_sources(dry_run, source_filter))

    print("\n── ChEMBL ────────────────────────────────────────────────────────────")
    all_rows.extend(iter_chembl_source(dry_run, source_filter))

    print("\n── BindingDB ─────────────────────────────────────────────────────────")
    all_rows.extend(iter_bindingdb_source(dry_run, source_filter))

    if not all_rows:
        print("\nNo rows processed. Check --source filter.")
        sys.exit(1)

    df = pd.DataFrame(all_rows)
    df["n_heavy_atoms"] = pd.to_numeric(df["n_heavy_atoms"], errors="coerce")
    df["mol_weight"]    = pd.to_numeric(df["mol_weight"],    errors="coerce")

    print(f"\n── Writing parquet: {parquet_out} ───────────────────────────────────")
    df.to_parquet(parquet_out, index=False)
    print(f"  {len(df):,} rows written")

    report_df = build_report(df)
    report_df.to_csv(report_out, index=False)
    print(f"  Report: {report_out}")

    # Summary — numbers derived from same df as report
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    total       = len(df)
    n_ok        = (df["std_status"] == "OK").sum()
    n_parse     = (df["std_status"] == "PARSE_FAIL").sum()
    n_sanitize  = (df["std_status"] == "SANITIZE_FAIL").sum()
    n_no_smiles = (df["std_status"] == "NO_SMILES").sum()
    n_uniq_ik   = df["inchi_key"].dropna().nunique()

    print(f"  Total input rows      : {total:>8,}")
    print(f"  OK                    : {n_ok:>8,}  ({100*n_ok/total:.1f}%)")
    print(f"  PARSE_FAIL            : {n_parse:>8,}  ({100*n_parse/total:.1f}%)")
    print(f"  SANITIZE_FAIL         : {n_sanitize:>8,}  ({100*n_sanitize/total:.1f}%)")
    print(f"  NO_SMILES             : {n_no_smiles:>8,}  ({100*n_no_smiles/total:.1f}%)")
    print(f"  Unique InChIKeys      : {n_uniq_ik:>8,}")

    print(f"\n  Per-source pct_ok:")
    for _, row in report_df.iterrows():
        flag = "  ⚠ " if row["pct_ok"] < 85 else "    "
        print(f"  {flag}{row['source']:30} AID {str(row['aid']):>10}  "
              f"{row['n_input']:>6} in  {row['pct_ok']:>5.1f}% OK")

    low_ok = report_df[report_df["pct_ok"] < 85]
    if not low_ok.empty:
        print(f"\n  ⚠  {len(low_ok)} source(s) below 85% OK — review before proceeding:")
        for _, r in low_ok.iterrows():
            print(f"    {r['source']} AID {r['aid']}: {r['pct_ok']}%  "
                  f"(parse={r['n_parse_fail']} sanitize={r['n_sanitize_fail']} no_smiles={r['n_no_smiles']})")
    else:
        print(f"\n  ✅ All sources ≥ 85% OK")

    print(f"\n  Output : {parquet_out}")
    print(f"  Report : {report_out}")


if __name__ == "__main__":
    main()
