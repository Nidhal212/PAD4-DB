#!/usr/bin/env python3
"""
PAD4-DB v2 Full Reproducibility Run
All outputs → data/interim_repro/ and data/processed_repro/
Gold files in data/interim/ and data/processed/ are NOT touched.
"""

import sys
import os
import time
import traceback
from pathlib import Path
import pandas as pd
import numpy as np

# ── must run from project root ────────────────────────────────────────────────
ROOT = Path("/home/nidhal/PAD4-db_V2")
os.chdir(ROOT)
sys.path.insert(0, str(ROOT))

LOG = ROOT / "outputs/audit/REPRODUCIBILITY_RUN.log"
LOG.parent.mkdir(parents=True, exist_ok=True)

def log(msg):
    ts = time.strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")

def fail(msg):
    log(f"STOP: {msg}")
    sys.exit(1)

def patch_and_exec(script_path, patches):
    """Read script source, apply path string patches, exec as __main__."""
    src = Path(script_path).read_text()
    for old, new in patches.items():
        src = src.replace(old, new)
    ns = {"__name__": "__main__", "__file__": str(script_path)}
    try:
        exec(compile(src, str(script_path), "exec"), ns)
    except SystemExit as e:
        if e.code and e.code != 0:
            raise RuntimeError(f"Script exited with code {e.code}") from e

# Path patches applied to ALL scripts.
# Strategy: replace the quoted path strings inside Path() calls.
PATCHES = {
    '"data/interim/standardized'  : '"data/interim_repro/standardized',
    '"data/interim/normalized'    : '"data/interim_repro/normalized',
    '"data/processed/'            : '"data/processed_repro/',
    # scaffold/cliff output tables → tables_repro (inputs to Phase 7)
    '"outputs/tables/05_scaffold' : '"outputs/tables_repro/05_scaffold',
    '"outputs/tables/05_cliff'    : '"outputs/tables_repro/05_cliff',
    '"outputs/tables/05_patent'   : '"outputs/tables_repro/05_patent',
    # QC report outputs (01, 02, 03b, 03c) → tables_repro (non-critical but clean)
    '"outputs/tables/0'           : '"outputs/tables_repro/0',
}

# ════════════════════════════════════════════════════════════════════════════
# PHASE 0 — Raw Data Inventory
# ════════════════════════════════════════════════════════════════════════════
log("=" * 65)
log("PHASE 0 — RAW DATA INVENTORY")
log("=" * 65)

import glob

# ChEMBL
chembl_files = list(ROOT.glob("data/raw/chembl/CHEMBL6111*.csv"))
if not chembl_files:
    fail("MISSING: data/raw/chembl/CHEMBL6111*.csv")
chembl_rows = sum(len(pd.read_csv(f, sep=";", low_memory=False)) for f in chembl_files)
status = "✅" if abs(chembl_rows - 4925) <= 5 else "❌"
log(f"{status} ChEMBL: {chembl_rows} rows (expected 4925)")
if "❌" in status:
    fail(f"ChEMBL row count mismatch: {chembl_rows}")

# BindingDB
bdb_file = ROOT / "data/raw/bindingdb/bindingdb_Q9UM07.tsv"
if not bdb_file.exists():
    fail(f"MISSING: {bdb_file}")
bdb_rows = len(pd.read_csv(bdb_file, sep="\t", low_memory=False))
status = "✅" if abs(bdb_rows - 3087) <= 5 else "❌"
log(f"{status} BindingDB: {bdb_rows} rows (expected 3087)")
if "❌" in status:
    fail(f"BindingDB row count mismatch: {bdb_rows}")

# PubChem dirs
pc_dirs = {
    "pubchem/confirmatory":       57,
    "pubchem/literature_derived": 11,
    "pubchem/secondary":          26,
    "hts":                         3,
}
for subdir, expected in pc_dirs.items():
    p = ROOT / "data/raw" / subdir
    if not p.exists():
        fail(f"MISSING dir: {p}")
    aids = list(p.glob("*.csv")) + list(p.glob("*.txt"))
    status = "✅" if len(aids) == expected else "❌"
    log(f"{status} PubChem/{subdir}: {len(aids)} AIDs (expected {expected})")
    if "❌" in status:
        fail(f"AID count mismatch in {subdir}: {len(aids)}")

log("✅ PHASE 0 PASS — all raw files present")
log("")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 1 — SMILES Standardization
# ════════════════════════════════════════════════════════════════════════════
log("=" * 65)
log("PHASE 1 — SMILES STANDARDIZATION")
log("=" * 65)
t0 = time.time()

patch_and_exec(
    ROOT / "scripts/01_standardize/01_standardize_smiles.py",
    PATCHES
)

elapsed = time.time() - t0
log(f"Phase 1 runtime: {elapsed:.0f}s")

# Verify
std_repro = pd.read_parquet(ROOT / "data/interim_repro/standardized/standardized_compounds.parquet")
std_gold  = pd.read_parquet(ROOT / "data/interim/standardized/standardized_compounds.parquet")

log(f"Repro rows: {len(std_repro)}")
log(f"Gold rows:  {len(std_gold)}")
if len(std_repro) != 341282:
    fail(f"Row count mismatch: {len(std_repro)} ≠ 341282")

repro_iks = set(std_repro["inchi_key"].dropna())
gold_iks  = set(std_gold["inchi_key"].dropna())
only_repro = repro_iks - gold_iks
only_gold  = gold_iks - repro_iks
log(f"InChIKeys only in repro: {len(only_repro)}")
log(f"InChIKeys only in gold:  {len(only_gold)}")
if only_repro or only_gold:
    if only_repro:
        log(f"  New IKs: {list(only_repro)[:3]}")
    if only_gold:
        log(f"  Lost IKs: {list(only_gold)[:3]}")
    fail("InChIKey set mismatch in Phase 1")

log("✅ PHASE 1 PASS")
log("")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 2 — Activity Normalization
# ════════════════════════════════════════════════════════════════════════════
log("=" * 65)
log("PHASE 2 — ACTIVITY NORMALIZATION")
log("=" * 65)
t0 = time.time()

patch_and_exec(
    ROOT / "scripts/02_normalize/02_normalize_activities.py",
    PATCHES
)

elapsed = time.time() - t0
log(f"Phase 2 runtime: {elapsed:.0f}s")

norm_repro = pd.read_parquet(ROOT / "data/interim_repro/normalized/normalized_activities.parquet")
log(f"Repro rows: {len(norm_repro)}")
if len(norm_repro) != 341282:
    fail(f"Normalized row count mismatch: {len(norm_repro)} ≠ 341282")

# pIC50 arithmetic check
if "pic50" in norm_repro.columns and "value_nm" in norm_repro.columns:
    valid = norm_repro.dropna(subset=["pic50", "value_nm"])
    expected = -np.log10(valid["value_nm"] * 1e-9)
    max_err = (expected - valid["pic50"]).abs().max()
    log(f"pIC50 arithmetic max error: {max_err:.6f}")
    if max_err > 0.001:
        fail(f"pIC50 arithmetic error too large: {max_err}")

log("✅ PHASE 2 PASS")
log("")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 3 — Replicate Aggregation + Space Split
# ════════════════════════════════════════════════════════════════════════════
log("=" * 65)
log("PHASE 3 — REPLICATE AGGREGATION + SPACE SPLIT")
log("=" * 65)
t0 = time.time()

# Step 3: replicate aggregate (reads normalized, writes replicate_aggregated)
log("Running 03_replicate_aggregate.py ...")
patch_and_exec(ROOT / "scripts/03_aggregate/03_replicate_aggregate.py", PATCHES)

# Step 3a: split spaces (reads replicate_aggregated, writes potency + hts space)
log("Running 03a_split_spaces.py ...")
patch_and_exec(ROOT / "scripts/03_aggregate/03a_split_spaces.py", PATCHES)

# Step 3b: logspace QC (read-only validation)
log("Running 03b_logspace_qc.py ...")
patch_and_exec(ROOT / "scripts/03_aggregate/03b_logspace_qc.py", PATCHES)

# Step 3c: SMILES integrity (read-only validation)
log("Running 03c_smiles_integrity.py ...")
patch_and_exec(ROOT / "scripts/03_aggregate/03c_smiles_integrity.py", PATCHES)

elapsed = time.time() - t0
log(f"Phase 3 runtime: {elapsed:.0f}s")

pot_repro = pd.read_parquet(ROOT / "data/interim_repro/normalized/potency_space.parquet")
log(f"Potency space rows: {len(pot_repro)}")
if len(pot_repro) != 7319:
    fail(f"Potency space mismatch: {len(pot_repro)} ≠ 7319")

if "inchi_key" in pot_repro.columns:
    n_unique = pot_repro["inchi_key"].nunique()
    log(f"Unique compounds in potency space: {n_unique}  (expected 3093)")
    if n_unique != 3093:
        fail(f"Potency space unique compound mismatch: {n_unique}")

log("✅ PHASE 3 PASS")
log("")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 4 — Deduplication and Assembly
# ════════════════════════════════════════════════════════════════════════════
log("=" * 65)
log("PHASE 4 — DEDUPLICATION AND ASSEMBLY")
log("=" * 65)
t0 = time.time()

log("Running 04_dedup_and_assemble.py ...")
patch_and_exec(ROOT / "scripts/04_dedup/04_dedup_and_assemble.py", PATCHES)

log("Running 04b_add_independence_scores.py ...")
patch_and_exec(ROOT / "scripts/04_dedup/04b_add_independence_scores.py", PATCHES)

elapsed = time.time() - t0
log(f"Phase 4 runtime: {elapsed:.0f}s")

repro = pd.read_parquet(ROOT / "data/processed_repro/pad4_compounds.parquet")
gold  = pd.read_parquet(ROOT / "data/processed/pad4_compounds.parquet")

log(f"Repro compounds: {len(repro)}")
log(f"Gold  compounds: {len(gold)}")
if len(repro) != 3093:
    fail(f"Compound count mismatch: {len(repro)} ≠ 3093")

repro_iks = set(repro["inchi_key"])
gold_iks  = set(gold["inchi_key"])
only_repro = repro_iks - gold_iks
only_gold  = gold_iks - repro_iks
log(f"InChIKeys only in repro: {len(only_repro)}")
log(f"InChIKeys only in gold:  {len(only_gold)}")
if only_repro:
    log(f"  New: {list(only_repro)[:3]}")
if only_gold:
    log(f"  Lost: {list(only_gold)[:3]}")
if only_repro or only_gold:
    fail("InChIKey set mismatch in Phase 4")

# pIC50 concordance
pic50_col = "pic50_consensus" if "pic50_consensus" in repro.columns else "pic50"
merged = repro[["inchi_key", pic50_col]].merge(
    gold[["inchi_key", pic50_col]], on="inchi_key", suffixes=("_repro", "_gold"))
delta = (merged[f"{pic50_col}_repro"] - merged[f"{pic50_col}_gold"]).abs()
log(f"pIC50 concordance: max|Δ|={delta.max():.6f}  mean|Δ|={delta.mean():.6f}")
log(f"  Compounds with |Δ|>0.001: {(delta > 0.001).sum()}")
if delta.max() > 0.001:
    bad = merged[delta > 0.001]
    log(f"  First discrepant: {bad.iloc[0]['inchi_key']}  "
        f"repro={bad.iloc[0][f'{pic50_col}_repro']:.4f}  "
        f"gold={bad.iloc[0][f'{pic50_col}_gold']:.4f}")
    fail(f"pIC50 discrepancy: max|Δ|={delta.max():.6f}")

if "source_independence_score" in repro.columns:
    n_multi = (repro["source_independence_score"] >= 0.6).sum()
    log(f"is_true_multi_source (score>=0.6): {n_multi}  (expected 528)")
    if n_multi != 528:
        fail(f"Multi-source count mismatch: {n_multi} ≠ 528")

log("✅ PHASE 4 PASS")
log("")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 5 — Scaffold and Cliff Analysis
# ════════════════════════════════════════════════════════════════════════════
log("=" * 65)
log("PHASE 5 — SCAFFOLD AND CLIFF ANALYSIS")
log("=" * 65)
t0 = time.time()

log("Running 05_scaffold_and_cliffs.py  (this is the slow step)...")
patch_and_exec(ROOT / "scripts/05_cliffs/05_scaffold_and_cliffs.py", PATCHES)

elapsed = time.time() - t0
log(f"Phase 5 runtime: {elapsed:.0f}s")

cliffs_repro = pd.read_parquet(ROOT / "data/processed_repro/activity_cliffs.parquet")
cliffs_gold  = pd.read_parquet(ROOT / "data/processed/activity_cliffs.parquet")

log(f"Repro cliff pairs: {len(cliffs_repro)}")
log(f"Gold  cliff pairs: {len(cliffs_gold)}")
if len(cliffs_repro) != 867:
    fail(f"Cliff count mismatch: {len(cliffs_repro)} ≠ 867")

# Use cliff_tier column if available, else delta_pic50 thresholds
if "cliff_tier" in cliffs_repro.columns:
    severe_repro   = cliffs_repro[cliffs_repro["cliff_tier"] == "severe"]
    moderate_repro = cliffs_repro[cliffs_repro["cliff_tier"] == "moderate"]
    broad_repro    = cliffs_repro[cliffs_repro["cliff_tier"] == "broad"]
else:
    dpic = cliffs_repro["delta_pic50"].abs()
    severe_repro   = cliffs_repro[dpic >= 2.0]
    moderate_repro = cliffs_repro[(dpic >= 1.5) & (dpic < 2.0)]
    broad_repro    = cliffs_repro[(dpic >= 1.0) & (dpic < 1.5)]

log(f"Severe:   {len(severe_repro)}  (expected 94)")
log(f"Moderate: {len(moderate_repro)}  (expected 193)")
log(f"Broad:    {len(broad_repro)}  (expected 580)")

for label, repro_n, expected in [
    ("Severe",   len(severe_repro),   94),
    ("Moderate", len(moderate_repro), 193),
    ("Broad",    len(broad_repro),    580),
]:
    if repro_n != expected:
        fail(f"{label} cliff count mismatch: {repro_n} ≠ {expected}")

max_delta = cliffs_repro["delta_pic50"].abs().max()
log(f"Max |ΔpIC50|: {max_delta:.4f}  (expected 3.045 ± 0.01)")
if abs(max_delta - 3.045) > 0.01:
    fail(f"Max delta mismatch: {max_delta:.4f}")

all_sev_iks = pd.concat([severe_repro["inchi_key_a"], severe_repro["inchi_key_b"]]).unique()
log(f"Unique compounds in severe cliffs: {len(all_sev_iks)}  (expected 99)")
if len(all_sev_iks) != 99:
    fail(f"Severe compound count mismatch: {len(all_sev_iks)}")

degree = pd.concat([severe_repro["inchi_key_a"], severe_repro["inchi_key_b"]]).value_counts()
log("Top 10 hub degrees:")
for ik, cnt in degree.head(10).items():
    log(f"  {ik}  degree={cnt}")

EXPECTED_HUBS = {
    "SMADULGDNOCLOP-GISFHXKWSA-N": 15,
    "RAVBZQAQTVGKIV-XBPDSQQVSA-N": 12,
    "UDCDEKJNAMHBFH-HSZRJFAPSA-N": 12,
    "DVCKJOQIVOGXEI-XMMPIXPASA-N": 11,
}
for ik, expected_deg in EXPECTED_HUBS.items():
    actual = int(degree.get(ik, 0))
    status = "✅" if actual == expected_deg else "❌"
    log(f"  {status} {ik[:20]}: degree={actual} (expected {expected_deg})")
    if actual != expected_deg:
        fail(f"Hub degree mismatch for {ik}: {actual} ≠ {expected_deg}")

sc_csv = ROOT / "outputs/tables_repro/05_scaffold_summary.csv"
if sc_csv.exists():
    sc = pd.read_csv(sc_csv)
    log(f"Scaffolds: {len(sc)}  (expected 1244)")
    log(f"Series (>=2): {(sc['n_compounds'] >= 2).sum()}  (expected 375)")
    log(f"Largest series: {sc['n_compounds'].max()}  (expected 174)")
    if len(sc) != 1244:
        fail(f"Scaffold count mismatch: {len(sc)}")
    if (sc["n_compounds"] >= 2).sum() != 375:
        fail(f"Series count mismatch: {(sc['n_compounds'] >= 2).sum()}")
    if sc["n_compounds"].max() != 174:
        fail(f"Largest series mismatch: {sc['n_compounds'].max()}")
else:
    log(f"WARNING: scaffold summary CSV not found at {sc_csv}")
    sc = None

log("✅ PHASE 5 PASS")
log("")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 6 — Tanimoto Recomputation
# ════════════════════════════════════════════════════════════════════════════
log("=" * 65)
log("PHASE 6 — TANIMOTO RECOMPUTATION (94 severe pairs)")
log("=" * 65)

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, DataStructs

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
smiles_map = dict(zip(repro["inchi_key"], repro["smiles_std"]))

tani_errors = []
for _, row in severe_repro.iterrows():
    smi_a = smiles_map.get(row["inchi_key_a"])
    smi_b = smiles_map.get(row["inchi_key_b"])
    if not smi_a or not smi_b:
        tani_errors.append(f"Missing SMILES: {row['inchi_key_a'][:14]}")
        continue
    mol_a = Chem.MolFromSmiles(smi_a)
    mol_b = Chem.MolFromSmiles(smi_b)
    if not mol_a or not mol_b:
        tani_errors.append(f"Parse fail: {row['inchi_key_a'][:14]}")
        continue
    fp_a = gen.GetFingerprint(mol_a)
    fp_b = gen.GetFingerprint(mol_b)
    tani_repro = DataStructs.TanimotoSimilarity(fp_a, fp_b)
    tani_gold  = abs(row["tanimoto"]) if "tanimoto" in row.index else None
    if tani_gold is not None and abs(tani_repro - tani_gold) > 0.002:
        tani_errors.append(
            f"Mismatch: {row['inchi_key_a'][:14]} repro={tani_repro:.4f} gold={tani_gold:.4f}")

log(f"Tanimoto recomputation errors: {len(tani_errors)}")
if tani_errors:
    for e in tani_errors[:5]:
        log(f"  {e}")
    fail("Tanimoto recomputation mismatch")
log("✅ PHASE 6 PASS — all 94 pairs agree within 0.002")
log("")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 7 — Full Comparison Table
# ════════════════════════════════════════════════════════════════════════════
log("=" * 65)
log("PHASE 7 — FULL REPRODUCIBILITY COMPARISON TABLE")
log("=" * 65)

def _n(series_or_val):
    try:
        return float(series_or_val)
    except Exception:
        return float("nan")

checks = [
    ("SAR compounds",             3093,   len(repro)),
    ("pIC50 min",                 2.00,   repro[pic50_col].min()),
    ("pIC50 max",                 8.52,   repro[pic50_col].max()),
    ("pIC50 mean",                6.550,  repro[pic50_col].mean()),
    ("pIC50 median",              6.845,  repro[pic50_col].median()),
    ("pIC50 std",                 0.992,  repro[pic50_col].std()),
    ("Patent-exclusive (PubChem only)", 233,
        (repro["source_list"] == "pubchem_confirmatory").sum()
        if "source_list" in repro.columns else -1),
    ("BindingDB-only compounds",  95,
        (repro["source_list"] == "bindingdb").sum()
        if "source_list" in repro.columns else -1),
    ("ChEMBL-only compounds",     10,
        (repro["source_list"] == "chembl").sum()
        if "source_list" in repro.columns else -1),
    ("is_true_multi_source (≥0.6)", 528,
        (repro["source_independence_score"] >= 0.6).sum()
        if "source_independence_score" in repro.columns else -1),
    ("Total cliff pairs",          867,   len(cliffs_repro)),
    ("Severe cliff pairs",          94,   len(severe_repro)),
    ("Moderate cliff pairs",       193,   len(moderate_repro)),
    ("Broad cliff pairs",          580,   len(broad_repro)),
    ("Severe cliff compounds",      99,   len(all_sev_iks)),
    ("Max |ΔpIC50|",             3.045,   round(cliffs_repro["delta_pic50"].abs().max(), 3)),
    ("Unique scaffolds",          1244,   len(sc) if sc is not None else -1),
    ("Series scaffolds (≥2)",      375,   int((sc["n_compounds"] >= 2).sum()) if sc is not None else -1),
    ("Largest series",             174,   int(sc["n_compounds"].max()) if sc is not None else -1),
    ("Hub A1 degree",               15,   int(degree.get("SMADULGDNOCLOP-GISFHXKWSA-N", -1))),
    ("Hub A2 degree",               12,   int(degree.get("RAVBZQAQTVGKIV-XBPDSQQVSA-N", -1))),
    ("Hub B1 degree",               12,   int(degree.get("UDCDEKJNAMHBFH-HSZRJFAPSA-N", -1))),
    ("Hub B2 degree",               11,   int(degree.get("DVCKJOQIVOGXEI-XMMPIXPASA-N", -1))),
]

hdr = f"{'Metric':<45} {'Expected':>10} {'Repro':>10}  Status"
log(hdr)
log("-" * 70)

all_pass = True
report_lines = [hdr, "-" * 70]

for metric, gold_val, repro_val in checks:
    try:
        if isinstance(gold_val, float):
            match = abs(float(repro_val) - gold_val) < 0.01
            rv_str = f"{float(repro_val):.3f}"
        else:
            match = int(repro_val) == int(gold_val)
            rv_str = str(int(repro_val))
    except Exception:
        match = False
        rv_str = "ERROR"
    status = "✅ PASS" if match else "❌ FAIL"
    if not match:
        all_pass = False
    line = f"{metric:<45} {str(gold_val):>10} {rv_str:>10}  {status}"
    log(line)
    report_lines.append(line)

log("=" * 70)
if all_pass:
    log("REPRODUCIBILITY CHECK: ALL PASS")
    log("Pipeline is fully reproducible from raw data.")
else:
    log("REPRODUCIBILITY CHECK: FAILURES DETECTED")

# Write report
report_path = ROOT / "outputs/audit/REPRODUCIBILITY_REPORT.txt"
with open(report_path, "w") as f:
    import datetime
    f.write(f"PAD4-DB v2 Reproducibility Report\n")
    f.write(f"Run: {datetime.datetime.now()}\n\n")
    for line in report_lines:
        f.write(line + "\n")
    f.write("\n")
    f.write("RESULT: ALL PASS\n" if all_pass else "RESULT: FAILURES DETECTED\n")

log(f"\nReport written to: {report_path}")
log("=" * 65)
