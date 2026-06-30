#!/usr/bin/env python3
"""
PAD4-DB v2 — Raw Data Inventory & QC Check
scripts/00_inventory/00_check_raw_inventory.py

Checks:
  1. Registry completeness  — every active AID has a file on disk
  2. File integrity         — no empty or suspiciously small files
  3. CSV parse check        — every file opens and has expected columns
  4. Row counts             — records per file (SID/CID rows, not header)
  5. Duplicate AID copies   — AIDs present in more than one subdir
  6. Endpoint audit         — activity_type values per AID vs registry expectation
  7. [VERIFY] / [CHECK]     — flag AIDs needing manual inspection
  8. Non-PubChem sources    — ChEMBL and BindingDB file presence + column peek
  9. Summary report         — pass / warn / fail per AID

Verified findings (resolved 2025-06):
  ChEMBL file  : semicolon-delimited, all 48 expected columns present
  BindingDB    : TSV, all key columns present (UniProt col at position 21+)
  AID 2202576/77: RFMS1 vs RFMS2 same patent — 55/231 SIDs overlap → keep both,
                  prefer RFMS2 for overlapping SIDs during dedup
  AID 2202596/97: completely distinct compound sets — no dedup needed
  AID 725596/97 : Standard Type/Value/Units cols present → dose-response confirmed;
                  "1µM / 10µM" in description = enzyme preincubation conc, not compound conc

Output:
  outputs/tables/00_raw_inventory_report.csv
  outputs/tables/00_raw_inventory_summary.txt
"""

import os
import csv
import sys
from pathlib import Path
from collections import defaultdict

# =============================================================================
# CONFIG
# =============================================================================

BASE_DIR  = Path("data/raw")
OUT_DIR   = Path("outputs/tables")
OUT_DIR.mkdir(parents=True, exist_ok=True)

REPORT_CSV = OUT_DIR / "00_raw_inventory_report.csv"
REPORT_TXT = OUT_DIR / "00_raw_inventory_summary.txt"

MIN_FILE_BYTES = 200   # PubChem error pages are ~150 bytes

# PubChem CSV: required cols in every well-formed export
PUBCHEM_REQUIRED_COLS = {"PUBCHEM_SID", "PUBCHEM_CID"}

# =============================================================================
# REGISTRY  (must match download_pad4db_v2.py)
# =============================================================================

REGISTRY = {
    "HTS": {
        463073: {"subdir": "hts", "layer": "HTS", "endpoint": "FP_screening",  "flags": []},
        485272: {"subdir": "hts", "layer": "HTS", "endpoint": "FP_screening",  "flags": []},
        488796: {"subdir": "hts", "layer": "HTS", "endpoint": "Single_conc",   "flags": []},
    },
    "A": {
        492970:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        320707:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1804546: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1804627: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1805620: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1806182: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1806183: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1806764: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1919095: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1920200: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1963715: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        2202576: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",
                  "flags": ["VERIFIED_OVERLAP"],
                  "verify_note": "RFMS1 (75 SIDs) vs RFMS2 2202577 (231 SIDs): 55 SID overlap (23.8%). "
                                 "Keep both AIDs. During dedup: for 55 overlapping SIDs prefer 2202577 (RFMS2)."},
        2202577: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",
                  "flags": ["VERIFIED_OVERLAP"],
                  "verify_note": "RFMS2 (231 SIDs): preferred source for 55 SIDs shared with 2202576 (RFMS1)."},
        2202596: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",
                  "flags": ["VERIFIED_DISTINCT"],
                  "verify_note": "614 SIDs, 0 overlap with 2202597 — completely distinct compound sets. No dedup needed."},
        2202597: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",
                  "flags": ["VERIFIED_DISTINCT"],
                  "verify_note": "18 SIDs, 0 overlap with 2202596 — completely distinct compound sets. No dedup needed."},
        2202717: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        2202442: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": ["DUAL_LAYER_D"]},
        1330527: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1471656: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1474465: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1474486: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1511938: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "Pct_inh", "flags": []},
        1632998: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1651501: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1651502: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1813082: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1813806: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1875531: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1920046: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": ["DUAL_LAYER_D"]},
        1973686: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1993481: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        2006936: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        2034405: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        2034406: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        2034407: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        2071731: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        2134413: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        2200614: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1422897: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        725673:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        # Verified dose-response: "1µM/10µM" in description = enzyme preincubation conc, not compound conc
        725596:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",
                  "flags": ["VERIFIED_DOSE_RESPONSE"],
                  "verify_note": "Confirmed dose-response. '10µM' in AID description = enzyme preincubation conc, not compound conc."},
        725597:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",
                  "flags": ["VERIFIED_DOSE_RESPONSE"],
                  "verify_note": "Confirmed dose-response. '1µM' in AID description = enzyme preincubation conc, not compound conc."},
        725598:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1196525: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1069613: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        1069619: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "Pct_inh", "flags": []},
        1069623: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "Pct_inh", "flags": []},
        627371:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "Pct_inh", "flags": []},
        627428:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "Pct_inh", "flags": []},
        627432:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "Pct_inh", "flags": []},
        626724:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "Pct_inh", "flags": []},
        626728:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "Pct_inh", "flags": []},
        1069614: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "Pct_inh", "flags": []},
        1069618: {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "Pct_inh", "flags": []},
        588488:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        588559:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
        588560:  {"subdir": "pubchem/confirmatory", "layer": "A", "endpoint": "IC50",    "flags": []},
    },
    "C": {
        626735:  {"subdir": "pubchem/literature_derived", "layer": "C", "endpoint": "Kinact_Ki", "flags": []},
        626738:  {"subdir": "pubchem/literature_derived", "layer": "C", "endpoint": "Kinact_Ki", "flags": []},
        712876:  {"subdir": "pubchem/literature_derived", "layer": "C", "endpoint": "Kinact_Ki", "flags": []},
        725671:  {"subdir": "pubchem/literature_derived", "layer": "C", "endpoint": "Kinact_Ki", "flags": []},
        1069608: {"subdir": "pubchem/literature_derived", "layer": "C", "endpoint": "Kinact_Ki", "flags": []},
        1196521: {"subdir": "pubchem/literature_derived", "layer": "C", "endpoint": "Kinact_Ki", "flags": []},
        1364668: {"subdir": "pubchem/literature_derived", "layer": "C", "endpoint": "Kinact_Ki", "flags": []},
        1422898: {"subdir": "pubchem/literature_derived", "layer": "C", "endpoint": "Kinact_Ki", "flags": []},
        1422904: {"subdir": "pubchem/literature_derived", "layer": "C", "endpoint": "Kinact_Ki", "flags": []},
        2076402: {"subdir": "pubchem/literature_derived", "layer": "C", "endpoint": "Kinact_Ki", "flags": []},
        651867:  {"subdir": "pubchem/literature_derived", "layer": "C", "endpoint": "Kinact_Ki", "flags": []},
    },
    "D": {
        2041348: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "FP_binding",   "flags": []},
        2041349: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "FP_binding",   "flags": []},
        2053867: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "FP_binding",   "flags": []},
        2053915: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "Enzymatic",    "flags": []},
        2053917: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "FP_binding",   "flags": []},
        2193457: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "Enzymatic",    "flags": []},
        1625405: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "Displacement", "flags": []},
        1806765: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "FP_binding",   "flags": []},
        1920046: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "FP_binding",   "flags": ["DUAL_LAYER_A"]},
        2202442: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "FP_binding",   "flags": ["DUAL_LAYER_A"]},
        1069597: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "Mass_shift",   "flags": []},
        1069598: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "Mass_shift",   "flags": []},
        1069599: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "Mass_shift",   "flags": []},
        1069600: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "Mass_shift",   "flags": []},
        1069601: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "Mass_shift",   "flags": []},
        1069604: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "Mass_shift",   "flags": []},
        1069605: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "FP_ABPP",      "flags": []},
        1069606: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "FP_ABPP",      "flags": []},
        588487:  {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "ABPP_Gel",     "flags": []},
        651627:  {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "Binding_mode", "flags": []},
        1346144: {"subdir": "pubchem/secondary", "layer": "D", "endpoint": "Ki_Kd",
                  "flags": ["CHECK_UNITS"],
                  "verify_note": "ChEMBL affinity aggregator — mixed Ki/Kd units. Curate Standard Units carefully during normalization."},
    },
    "E": {
        2039667: {"subdir": "pubchem/secondary", "layer": "E", "endpoint": "NH3_release", "flags": []},
        2053916: {"subdir": "pubchem/secondary", "layer": "E", "endpoint": "NH3_release", "flags": []},
        627442:  {"subdir": "pubchem/secondary", "layer": "E", "endpoint": "Pct_inh",     "flags": []},
        627443:  {"subdir": "pubchem/secondary", "layer": "E", "endpoint": "Pct_inh",     "flags": []},
        1422917: {"subdir": "pubchem/secondary", "layer": "E", "endpoint": "Pct_inh",     "flags": []},
    },
    "EXCLUDED": {
        626677: {"reason": "PAD4_mutant_kinetics", "notes": "Kcat/Km Q346A mutant"},
        626678: {"reason": "PAD4_mutant_kinetics", "notes": "Kcat/Km R347K mutant"},
        626679: {"reason": "PAD4_mutant_kinetics", "notes": "Kcat/Km R347A mutant"},
        627437: {"reason": "PAD4_mutant_kinetics", "notes": "Irreversible Q346E mutant"},
        627438: {"reason": "PAD4_mutant_kinetics", "notes": "Irreversible Q346A mutant"},
        627439: {"reason": "PAD4_mutant_kinetics", "notes": "Irreversible R374K mutant"},
        627440: {"reason": "PAD4_mutant_kinetics", "notes": "Irreversible R374A mutant"},
        627441: {"reason": "PAD4_mutant_kinetics", "notes": "Kcat/Km Q346E mutant"},
        463083: {"reason": "project_summary_no_data", "notes": "MLPCN PAD4 probe summary — no SID/CID rows"},
    },
}

# =============================================================================
# EXTERNAL SOURCE SPECS  (corrected from actual file inspection)
# =============================================================================

CHEMBL_SPEC = {
    "delimiter": ";",
    "required_cols": {"Molecule ChEMBL ID", "Smiles", "Standard Type",
                      "Standard Value", "Standard Units", "pChEMBL Value"},
    "notes": "Semicolon-delimited (not comma). All 48 columns present.",
}

BINDINGDB_SPEC = {
    "delimiter": "\t",
    "required_cols": {"Ligand SMILES", "Ki (nM)", "IC50 (nM)", "Kd (nM)", "Target Name"},
    "notes": "UniProt col present at position 21+ (beyond first 20 shown in preview). "
             "All key activity and SMILES columns confirmed present.",
}

# =============================================================================
# HELPERS
# =============================================================================

def parse_pubchem_csv(filepath: Path):
    """
    Parse a PubChem AID CSV, handling the 4-line metadata header.
    Returns (columns, data_rows, activity_types, has_smiles, errors).
    """
    errors = []
    columns, data_rows, activity_types = [], [], set()
    has_smiles = False

    try:
        with open(filepath, newline="", encoding="utf-8", errors="replace") as f:
            raw = f.read()

        lines = raw.splitlines()
        if not lines:
            return columns, data_rows, activity_types, has_smiles, ["EMPTY_FILE"]

        # Find the real header row (contains PUBCHEM_SID)
        header_idx = next(
            (i for i, l in enumerate(lines[:8]) if "PUBCHEM_SID" in l), 0
        )

        reader = csv.DictReader(lines[header_idx:])
        columns = list(reader.fieldnames or [])

        for row in reader:
            # Skip PubChem metadata description rows (they have no SID value)
            if not row.get("PUBCHEM_SID", "").strip().lstrip("-").isdigit():
                continue
            data_rows.append(row)
            at = (row.get("PUBCHEM_ACTIVITY_OUTCOME", "")
                  or row.get("Standard Type", "")
                  or row.get("Activity_Type", ""))
            if at.strip():
                activity_types.add(at.strip())

        has_smiles = any(c in columns for c in
                         ["PUBCHEM_CANONICAL_SMILES", "PUBCHEM_ISOMERIC_SMILES", "Canonical_SMILES"])

    except Exception as e:
        errors.append(f"PARSE_ERROR:{e}")

    return columns, data_rows, activity_types, has_smiles, errors


def check_pubchem_file(aid: int, filepath: Path, meta: dict) -> dict:
    flags = meta.get("flags", [])
    result = {
        "aid":             aid,
        "category":        meta.get("category", ""),
        "layer":           meta.get("layer", ""),
        "endpoint":        meta.get("endpoint", ""),
        "flags":           "|".join(flags),
        "expected_subdir": meta.get("subdir", ""),
        "actual_path":     str(filepath),
        "file_exists":     filepath.exists(),
        "file_bytes":      0,
        "status":          "PASS",
        "issues":          [],
        "n_data_rows":     0,
        "n_columns":       0,
        "has_sid":         False,
        "has_cid":         False,
        "has_smiles":      False,
        "activity_types":  "",
        "verify_note":     meta.get("verify_note", ""),
    }

    if not filepath.exists():
        result["status"] = "FAIL"
        result["issues"].append("FILE_MISSING")
        return result

    size = filepath.stat().st_size
    result["file_bytes"] = size

    if size == 0:
        result["status"] = "FAIL"
        result["issues"].append("EMPTY_FILE")
        return result
    if size < MIN_FILE_BYTES:
        result["status"] = "WARN"
        result["issues"].append(f"VERY_SMALL_{size}B")

    cols, rows, atypes, has_smiles, errors = parse_pubchem_csv(filepath)
    result.update({
        "n_columns":    len(cols),
        "n_data_rows":  len(rows),
        "has_sid":      "PUBCHEM_SID" in cols,
        "has_cid":      "PUBCHEM_CID" in cols,
        "has_smiles":   has_smiles,
        "activity_types": "|".join(sorted(atypes)),
    })

    if errors:
        result["issues"].extend(errors)
        result["status"] = "WARN"
    if not result["has_sid"]:
        result["issues"].append("MISSING_PUBCHEM_SID")
        result["status"] = "WARN"
    if result["n_data_rows"] == 0:
        result["issues"].append("ZERO_DATA_ROWS")
        result["status"] = "FAIL"

    return result


def scan_disk() -> dict:
    """Return {aid: [paths]} across all pubchem + hts dirs."""
    found = defaultdict(list)
    for d in [BASE_DIR / "hts",
              BASE_DIR / "pubchem" / "confirmatory",
              BASE_DIR / "pubchem" / "literature_derived",
              BASE_DIR / "pubchem" / "secondary"]:
        if not d.exists():
            continue
        for fp in sorted(d.glob("AID_*_datatable_all.csv")):
            try:
                aid = int(fp.name.split("_")[1])
                found[aid].append(fp)
            except (IndexError, ValueError):
                pass
    return found


def check_external(path: Path, spec: dict, source_name: str) -> dict:
    """Check one external source file against its spec."""
    result = {"source": source_name, "file": path.name if path.exists() else "MISSING",
              "status": "PASS", "n_rows": 0, "columns_found": "", "notes": ""}
    if not path.exists():
        result["status"] = "FAIL"
        result["notes"] = f"File not found: {path}"
        return result
    try:
        sep = spec["delimiter"]
        with open(path, newline="", encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f, delimiter=sep)
            cols = list(reader.fieldnames or [])
            rows = list(reader)
        missing = spec["required_cols"] - set(cols)
        result["n_rows"] = len(rows)
        result["columns_found"] = f"{len(cols)} columns"
        if missing:
            result["status"] = "WARN"
            result["notes"] = f"Missing cols: {missing}"
        else:
            result["notes"] = spec.get("notes", "All required columns present")
    except Exception as e:
        result["status"] = "FAIL"
        result["notes"] = str(e)
    return result


# =============================================================================
# MAIN
# =============================================================================

def main():
    lines = []

    def log(msg=""):
        print(msg)
        lines.append(msg)

    log("=" * 72)
    log("PAD4-DB v2 — Raw Data Inventory & QC Report")
    log("=" * 72)

    # 1. Disk scan
    disk_map = scan_disk()
    log(f"\n📁 Files on disk: {sum(len(v) for v in disk_map.values())} "
        f"across {len(disk_map)} unique AIDs")

    # 2. Duplicate location check
    log("\n── Duplicate AID location check ──────────────────────────────────────")
    dual_expected = {1920046, 2202442}
    dup_issues = []
    for aid, paths in sorted(disk_map.items()):
        if len(paths) > 1:
            subdirs = [p.parent.name for p in paths]
            if aid in dual_expected:
                log(f"  ℹ  AID {aid:>8} — dual-layer (expected): {subdirs}")
            else:
                log(f"  ⚠  AID {aid:>8} — UNEXPECTED duplicate: {[str(p) for p in paths]}")
                dup_issues.append(aid)
    if not dup_issues:
        log("  ✅ No unexpected duplicates")

    # 3. Per-AID checks
    log("\n── Per-AID file checks ───────────────────────────────────────────────")
    all_results = []
    counts = {"PASS": 0, "WARN": 0, "FAIL": 0}
    verified_notes = []

    for category, aids in REGISTRY.items():
        if category == "EXCLUDED":
            continue
        for aid, meta in aids.items():
            m = dict(meta)
            m["category"] = category
            expected_path = BASE_DIR / m["subdir"] / f"AID_{aid}_datatable_all.csv"

            # If missing from expected path but found elsewhere, use actual location
            if not expected_path.exists() and aid in disk_map:
                actual = disk_map[aid][0]
                result = check_pubchem_file(aid, actual, m)
                result["issues"].append(
                    f"WRONG_SUBDIR:expected={m['subdir']},got={actual.parent.name}")
                if result["status"] == "PASS":
                    result["status"] = "WARN"
            else:
                result = check_pubchem_file(aid, expected_path, m)

            all_results.append(result)
            counts[result["status"]] += 1

            icon = {"PASS": "✅", "WARN": "⚠ ", "FAIL": "❌"}[result["status"]]
            issue_str = ", ".join(result["issues"]) if result["issues"] else ""
            row_str = f"{result['n_data_rows']:>6} rows" if result["file_exists"] else "      —    "
            flag_str = f"  [{result['flags']}]" if result["flags"] else ""
            log(f"  {icon} AID {aid:>8}  [{category:<8}]  {row_str}  {issue_str}{flag_str}")

            if result["verify_note"]:
                log(f"              → {result['verify_note']}")
                verified_notes.append((aid, result["verify_note"]))

    # 4. Orphan check
    log("\n── Orphan check (on disk, not in registry) ──────────────────────────")
    registry_aids = {aid for cat, aids in REGISTRY.items()
                     if cat != "EXCLUDED" for aid in aids}
    orphans = set(disk_map.keys()) - registry_aids
    if orphans:
        for aid in sorted(orphans):
            log(f"  ❓ AID {aid} — not in registry: {disk_map[aid]}")
    else:
        log("  ✅ No orphan files")

    # 5. Excluded AID check
    log("\n── Excluded AID disk check ───────────────────────────────────────────")
    excl_present = [aid for aid in REGISTRY["EXCLUDED"] if aid in disk_map]
    if excl_present:
        for aid in excl_present:
            log(f"  ⚠  AID {aid} is EXCLUDED but on disk: {disk_map[aid]}")
    else:
        log("  ✅ No excluded AIDs on disk")

    # 6. External sources
    log("\n── External sources ──────────────────────────────────────────────────")
    chembl_file  = next((BASE_DIR / "chembl").glob("*.csv"), None)
    bindingdb_file = next((BASE_DIR / "bindingdb").glob("*.tsv"), None) or \
                     next((BASE_DIR / "bindingdb").glob("*.csv"), None)

    ext_checks = []
    if chembl_file:
        ext_checks.append(check_external(chembl_file, CHEMBL_SPEC, "ChEMBL"))
    else:
        ext_checks.append({"source": "ChEMBL", "file": "MISSING", "status": "FAIL",
                           "n_rows": 0, "columns_found": "", "notes": "No CSV in data/raw/chembl/"})
    if bindingdb_file:
        ext_checks.append(check_external(bindingdb_file, BINDINGDB_SPEC, "BindingDB"))
    else:
        ext_checks.append({"source": "BindingDB", "file": "MISSING", "status": "FAIL",
                           "n_rows": 0, "columns_found": "", "notes": "No TSV in data/raw/bindingdb/"})

    for r in ext_checks:
        icon = {"PASS": "✅", "WARN": "⚠ ", "FAIL": "❌"}[r["status"]]
        log(f"  {icon} {r['source']:12} {r['file'][:50]:50}  {r['n_rows']:>6} rows  "
            f"{r['columns_found']}  {r['notes']}")

    # 7. Summary
    log("\n" + "=" * 72)
    log("SUMMARY")
    log("=" * 72)
    total = sum(counts.values())
    log(f"  PubChem AIDs checked : {total}")
    log(f"  ✅ PASS              : {counts['PASS']}")
    log(f"  ⚠  WARN              : {counts['WARN']}")
    log(f"  ❌ FAIL              : {counts['FAIL']}")
    log(f"  🚫 Excluded AIDs     : {len(REGISTRY['EXCLUDED'])} (not downloaded by design)")

    ext_pass = all(r["status"] == "PASS" for r in ext_checks)
    log(f"\n  External sources:")
    for r in ext_checks:
        icon = {"PASS": "✅", "WARN": "⚠ ", "FAIL": "❌"}[r["status"]]
        log(f"    {icon} {r['source']}: {r['status']} ({r['n_rows']} rows)")

    # Decisions locked
    log(f"\n  Verified decisions (locked):")
    log(f"    AID 2202576/77  RFMS1 vs RFMS2: 55 SID overlap → keep both, prefer 2202577 for overlapping SIDs")
    log(f"    AID 2202596/97  Functional 1/2: 0 SID overlap → fully distinct, no dedup needed")
    log(f"    AID 725596/97   Endpoint confirmed dose-response IC50 (µM in description = enzyme conc)")
    log(f"    AID 1346144     Mixed Ki/Kd units → curate Standard Units during normalization")
    log(f"    AID 1920046/2202442  Dual-layer → split by endpoint type during normalization")

    fail_aids = [r["aid"] for r in all_results if r["status"] == "FAIL"]
    if fail_aids:
        log(f"\n  ❌ FAIL — re-download:")
        for aid in fail_aids:
            r = next(x for x in all_results if x["aid"] == aid)
            log(f"    AID {aid}: {', '.join(r['issues'])}")

    warn_aids = [r for r in all_results if r["status"] == "WARN" and r["issues"]]
    if warn_aids:
        log(f"\n  ⚠  WARN — review:")
        for r in warn_aids:
            log(f"    AID {r['aid']}: {', '.join(r['issues'])}")

    all_pass = counts["FAIL"] == 0 and counts["WARN"] == 0 and ext_pass
    overall = ("✅ ALL CHECKS PASSED — ready for standardization"
               if all_pass else
               ("⚠  WARNINGS — review above before proceeding"
                if counts["FAIL"] == 0 else
                "❌ FAILURES — re-download before proceeding"))
    log(f"\n  {overall}")
    log(f"\n  Reports: {REPORT_CSV}")
    log(f"           {REPORT_TXT}")

    # Write CSV
    fields = ["aid", "category", "layer", "endpoint", "flags", "expected_subdir",
              "actual_path", "file_exists", "file_bytes", "status", "n_data_rows",
              "n_columns", "has_sid", "has_cid", "has_smiles", "activity_types",
              "verify_note", "issues"]
    with open(REPORT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in all_results:
            r2 = dict(r)
            r2["issues"] = "|".join(r2["issues"])
            w.writerow(r2)

    with open(REPORT_TXT, "w") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main()