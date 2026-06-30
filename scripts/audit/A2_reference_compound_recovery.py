#!/usr/bin/env python3
"""
PAD4-DB v2 — Reference Compound Recovery Audit
scripts/audit/A2_reference_compound_recovery.py

Searches 14 manually curated PAD4 inhibitors across all raw source files
and traces each through the pipeline. VERIFICATION ONLY — no files modified.

Output: outputs/audit/A2_reference_recovery.json
        (full report also printed to terminal)
"""

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.MolStandardize import rdMolStandardize

# ── Paths ──────────────────────────────────────────────────────────────────
RAW_DIRS = [
    Path("data/raw/hts"),
    Path("data/raw/pubchem/confirmatory"),
    Path("data/raw/pubchem/literature_derived"),
    Path("data/raw/pubchem/secondary"),
]
CHEMBL_CSV  = Path("data/raw/chembl/CHEMBL6111_Protein-arginine deiminase type_4.csv")
BDB_TSV     = Path("data/raw/bindingdb/bindingdb_Q9UM07.tsv")
STD_PARQ    = Path("data/interim/standardized/standardized_compounds.parquet")
NORM_PARQ   = Path("data/interim/normalized/normalized_activities.parquet")
COMP_PARQ   = Path("data/processed/pad4_compounds.parquet")
OUT_JSON    = Path("outputs/audit/A2_reference_recovery.json")

for p in [CHEMBL_CSV, BDB_TSV, STD_PARQ, NORM_PARQ, COMP_PARQ]:
    if not p.exists():
        sys.exit(f"ERROR: required input not found: {p}")

OUT_JSON.parent.mkdir(parents=True, exist_ok=True)

# =============================================================================
# Reference compound definitions
# =============================================================================

REF_COMPOUNDS = [
    ("Streptonigrin",  "PVYJZLYGTZKPJE-UHFFFAOYSA-N",  5298,
     "CC1=C(C(=C(N=C1C(=O)O)C2=NC3=C(C=C2)C(=O)C(=C(C3=O)N)OC)N)C4=C(C(=C(C=C4)OC)OC)O",
     "HIGH", "Pan-PAD IC50=2.5µM; streptonigrin AID cluster"),

    ("Cl-amidine",     "BPWATVWOHQZVRP-NSHDSACASA-N",  24970878,
     "C1=CC=C(C=C1)C(=O)N[C@@H](CCCN=C(CCl)N)C(=O)N",
     "HIGH", "Pan-PAD IC50=5.9µM; Cl-amidine AID cluster"),

    ("o-F-Amidine",    "HBEIARVCIYYMOR-UHFFFAOYSA-N",  54669783,
     "C1=CC=C(C(=C1)C(=O)NC(CCCN=C(CF)N)C(=O)N)C(=O)O",
     "HIGH", "PAD4 inhibitor; crystal structure 3B1U"),

    ("F-Amidine",      "OLFDULIIJWCYCK-NSHDSACASA-N",  11589584,
     "C1=CC=C(C=C1)C(=O)N[C@@H](CCCN=C(CF)N)C(=O)N",
     "HIGH", "PAD1/4 IC50=21.6µM; F-amidine AID cluster"),

    ("GSK484",         "MULKOGJHUZTANI-ADMBKAPUSA-N",  168896895,
     "O=C(N1C[C@H](N)[C@H](O)CC1)C2=CC(OC)=C3C(N=C(C(N4CC5CC5)=CC6=C4C=CC=C6)N3C)=C2.[H]Cl",
     "MEDIUM", "Selective PAD4 IC50=50nM; ChEMBL CHEMBL3545376"),

    ("GSK199",         "KRGMIOKDGHBYQE-UNTBIKODSA-N",  86340155,
     "CCN1C(=CC2=C1N=CC=C2)C3=NC4=C(N3C)C(=CC(=C4)C(=O)N5CCCC(C5)N)OC.Cl",
     "MEDIUM", "Selective PAD4 IC50=200nM; ChEMBL CHEMBL3545375"),

    ("Amodiaquine",    "OVCDSSHSILBFBN-UHFFFAOYSA-N",  2165,
     "CCN(CC)CC1=C(C=CC(=C1)NC2=C3C=CC(=CC3=NC=C2)Cl)O",
     "MEDIUM", "Repositioning candidate; old compound"),

    ("Pyroxamide",     "PTJGLFIIZFVFJV-UHFFFAOYSA-N",  4996,
     "C1=CC(=CN=C1)NC(=O)CCCCCCC(=O)NO",
     "MEDIUM", "Repositioning candidate"),

    ("BB-Cl-Amidine",  "YDOAWJHYHGBQFI-QHCPKHFHSA-N",  129021946,
     "C1=CC=C(C=C1)C2=CC=C(C=C2)C(=O)NC(CCCN=C(CCl)N)C3=NC4=CC=CC=C4N3",
     "LOW", "Irreversible pan-PAD; 2016+ compound"),

    ("TDFA",           "SOZMHIJABUOUSN-ORMVGFHCSA-N",  121513865,
     "[H]N([C@@H]([C@@H](C)O)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CCCNC(=N)CF)C(N)=O)C(C)=O",
     "LOW", "Selective PAD4 IC50=2.3µM; peptide inhibitor"),

    ("BMS-P5",         None,                            118158953,
     "C[C@H]1CC[C@H](CN1C(=O)C2=CC3=C(C(=C2)OC)N(C(=N3)C4=CC5=C(N4CC6CC6)N=CC=C5)C)N",
     "LOW", "Selective PAD4 IC50=98nM; patent compound"),

    ("JBI-589",        "DUVCPNSLXBKGOK-XMMPIXPASA-N",  138578805,
     None,
     "LOW", "2023 compound, likely too recent"),

    ("PAD-PF1",        None,                            None,
     "O=C(N1CCCC1)C(N)C2CCN(C3=NC=NC4=C3C(C)=C(C)N4C5=CC=C(Br)C=C5)CC2",
     "LOW", "Allosteric PAD4 IC50=15.9µM; no PubChem CID"),

    ("AFM-30a",        None,                            None,
     "COC1=C2C(N(C([C@H](CCCNC(CF)=N)NC(C3=C4C(CNC4=O)=CC=C3)=O)=N2)C)=CC=C1",
     "EXCLUDE", "PAD2-specific inhibitor; should NOT be in PAD4 potency model"),
]

# ChEMBL molecule IDs for targeted search
CHEMBL_IDS = {
    "GSK484":        ["CHEMBL3545376", "CHEMBL4539512", "CHEMBL5081624"],
    "GSK199":        ["CHEMBL3545375"],
    "Cl-amidine":    ["CHEMBL1213462"],
    "o-F-Amidine":   ["CHEMBL1213463"],
    "Streptonigrin": ["CHEMBL14432"],
    "Amodiaquine":   ["CHEMBL633"],
}

NAME_VARIANTS = [
    "streptonigrin", "gsk484", "gsk199", "gsk-484", "gsk-199",
    "cl-amidine", "clamidine", "cl_amidine", "o-f-amidine",
    "f-amidine", "bb-cl", "amodiaquine", "pyroxamide",
    "bms-p5", "jbi-589", "jbi589", "tdfa", "pad-pf1",
]

# =============================================================================
# RDKit standardization (identical to Step 01)
# =============================================================================

_fc  = rdMolStandardize.LargestFragmentChooser()
_unc = rdMolStandardize.Uncharger()


def std_ik(smiles: str | None) -> str | None:
    """Standardize SMILES → InChIKey using same pipeline as Step 01."""
    if not smiles:
        return None
    # Strip Daylight extended annotations
    idx = smiles.find(" |")
    if idx != -1:
        smiles = smiles[:idx]
    smiles = smiles.strip()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        mol = _fc.choose(mol)
        mol = _unc.uncharge(mol)
        return inchi.MolToInchiKey(mol)
    except Exception:
        return None


def std_ik_safe(smiles: str | None) -> str | None:
    try:
        return std_ik(smiles)
    except Exception:
        return None


# =============================================================================
# Step 0 — Pre-compute standardized InChIKeys for all ref compounds
# =============================================================================

def precompute_ref_inchikeys(refs: list) -> dict:
    """Compute InChIKey from SMILES; compare to manual InChIKey."""
    print("\n=== STEP 0: PRE-COMPUTE REF INCHIKEYS ===")
    result = {}
    for name, manual_ik, cid, smiles, priority, notes in refs:
        computed_ik = std_ik_safe(smiles) if smiles else None
        match = None
        if manual_ik and computed_ik:
            match = (manual_ik == computed_ik)
        elif manual_ik and not computed_ik:
            match = None   # can't compare
        # For lookup: prefer computed (salt-stripped) if available, else manual
        effective_ik = computed_ik if computed_ik else manual_ik
        result[name] = {
            "manual_ik":   manual_ik,
            "computed_ik": computed_ik,
            "effective_ik": effective_ik,
            "ik_match":    match,
        }
        if match is False:
            print(f"  MISMATCH {name}:")
            print(f"    manual:   {manual_ik}")
            print(f"    computed: {computed_ik}")
            print(f"    → using computed (salt-stripped) for pipeline search")
        elif computed_ik and match:
            print(f"  MATCH     {name}: {computed_ik}")
        elif computed_ik and match is None:
            print(f"  COMPUTED  {name}: {computed_ik} (no manual IK)")
        else:
            print(f"  NO SMILES {name}: manual={manual_ik}")
    return result


# =============================================================================
# Step 1 — PubChem raw CSV search
# =============================================================================

def _parse_pubchem_csv(path: Path) -> tuple[list[str], list[dict]]:
    """Return (headers, rows) parsing PubChem 4-line metadata header."""
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = []
        for i, line in enumerate(f):
            lines.append(line.rstrip("\n"))
            if i > 10:
                break
    # Find header row (contains PUBCHEM_SID)
    hdr_idx = None
    for i, line in enumerate(lines[:8]):
        if "PUBCHEM_SID" in line:
            hdr_idx = i
            break
    if hdr_idx is None:
        return [], []

    df = pd.read_csv(path, skiprows=hdr_idx, dtype=str, low_memory=False)
    # Drop metadata rows (SID not a digit string)
    if "PUBCHEM_SID" in df.columns:
        df = df[df["PUBCHEM_SID"].str.match(r"^\d+$", na=False)]
    return list(df.columns), df.to_dict("records")


def search_pubchem_raw(refs: list, ref_iks: dict) -> dict[str, list]:
    """Search all PubChem CSV files by CID and by pre-computed InChIKey."""
    print("\n=== SEARCH 1: PubChem raw CSV files ===")

    # Build fast lookup sets
    cid_to_name:  dict[str, str] = {}
    ik_to_name:   dict[str, str] = {}
    for name, manual_ik, cid, smiles, priority, notes in refs:
        if cid:
            cid_to_name[str(cid)] = name
        info = ref_iks[name]
        if info["effective_ik"]:
            ik_to_name[info["effective_ik"]] = name
        # also index manual if different from computed
        if info["manual_ik"] and info["manual_ik"] != info["effective_ik"]:
            ik_to_name[info["manual_ik"]] = name

    # Load CID→InChIKey from standardized_compounds (avoid re-standardizing all SMILES)
    std_df = pd.read_parquet(
        STD_PARQ,
        columns=["source", "aid", "cid_or_id", "inchi_key"],
    )
    cid_ik_map: dict[str, str] = {}
    for _, row in std_df[std_df["inchi_key"].notna()].iterrows():
        cid_ik_map[str(row["cid_or_id"])] = str(row["inchi_key"])

    # Invert: ik → [(source, aid, cid_or_id)]
    ik_loc: dict[str, list] = {}
    for _, row in std_df[std_df["inchi_key"].notna()].iterrows():
        ik = str(row["inchi_key"])
        ik_loc.setdefault(ik, []).append((str(row["source"]), str(row["aid"]), str(row["cid_or_id"])))

    hits: dict[str, list] = {name: [] for name, *_ in refs}
    value_cols = ["Standard Type", "Standard Relation", "Standard Value", "Standard Units",
                  "PUBCHEM_ACTIVITY_OUTCOME", "PUBCHEM_ACTIVITY_SCORE"]
    smiles_cols = ["PUBCHEM_EXT_DATASOURCE_SMILES", "PUBCHEM_OPENEYE_ISO_SMILES", "PUBCHEM_CANONICAL_SMILES"]

    for raw_dir in RAW_DIRS:
        for csv_path in sorted(raw_dir.glob("AID_*_datatable_all.csv")):
            aid_str = csv_path.name.split("_")[1]
            try:
                df = pd.read_csv(csv_path, skiprows=lambda i: i in (1, 2, 3),
                                 dtype=str, low_memory=False)
                # Remove metadata rows
                if "PUBCHEM_SID" in df.columns:
                    df = df[df["PUBCHEM_SID"].str.match(r"^\d+$", na=False)]
            except Exception as e:
                print(f"  WARN: could not load {csv_path.name}: {e}")
                continue

            if df.empty:
                continue

            cid_col = "PUBCHEM_CID" if "PUBCHEM_CID" in df.columns else None

            for _, row in df.iterrows():
                # Priority (a): CID match
                name_hit = None
                method = None
                if cid_col:
                    cid_val = str(row.get(cid_col, "")).strip()
                    if cid_val in cid_to_name:
                        name_hit = cid_to_name[cid_val]
                        method = "HIT_CID"
                    elif cid_val in cid_ik_map:
                        ik = cid_ik_map[cid_val]
                        if ik in ik_to_name:
                            name_hit = ik_to_name[ik]
                            method = "HIT_CID_VIA_IK"

                # Priority (b): SMILES match via standardization
                if name_hit is None:
                    for sc in smiles_cols:
                        if sc not in df.columns:
                            continue
                        smi = str(row.get(sc, "")).strip()
                        if not smi or smi == "nan":
                            continue
                        ik = std_ik_safe(smi)
                        if ik and ik in ik_to_name:
                            name_hit = ik_to_name[ik]
                            method = f"HIT_SMILES({sc})"
                            break

                if name_hit is None:
                    continue

                rec = {
                    "aid":              aid_str,
                    "sid":              str(row.get("PUBCHEM_SID", "")),
                    "cid":              str(row.get(cid_col, "")) if cid_col else "",
                    "match_method":     method,
                    "activity_outcome": str(row.get("PUBCHEM_ACTIVITY_OUTCOME", "")),
                    "activity_score":   str(row.get("PUBCHEM_ACTIVITY_SCORE", "")),
                    "standard_type":    str(row.get("Standard Type", "")),
                    "standard_relation":str(row.get("Standard Relation", "")),
                    "standard_value":   str(row.get("Standard Value", "")),
                    "standard_units":   str(row.get("Standard Units", "")),
                }
                hits[name_hit].append(rec)

    for name, h in hits.items():
        if h:
            aids = sorted({x["aid"] for x in h})
            print(f"  {name}: {len(h)} raw rows across AIDs {aids}")
        else:
            print(f"  {name}: NOT FOUND in PubChem raw files")

    return hits


# =============================================================================
# Step 2 — ChEMBL search
# =============================================================================

def search_chembl(refs: list, ref_iks: dict) -> dict[str, dict]:
    print("\n=== SEARCH 2: ChEMBL ===")
    chembl = pd.read_csv(CHEMBL_CSV, sep=";", dtype=str)
    ik_to_name: dict[str, str] = {
        info["effective_ik"]: name
        for name in (r[0] for r in refs)
        for info in [ref_iks[name]]
        if info["effective_ik"]
    }

    results: dict[str, dict] = {name: {} for name, *_ in refs}

    for _, row in chembl.iterrows():
        mol_id    = str(row.get("Molecule ChEMBL ID", "")).strip()
        mol_name  = str(row.get("Molecule Name", "")).strip().lower()
        smiles    = str(row.get("Smiles", "")).strip()

        matched_name = None
        method = None

        # (a) ChEMBL ID match
        for name, chembl_ids in CHEMBL_IDS.items():
            if mol_id in chembl_ids:
                matched_name = name
                method = f"HIT_CHEMBLID({mol_id})"
                break

        # (b) Name variant match
        if matched_name is None:
            for variant in NAME_VARIANTS:
                if variant in mol_name:
                    # Map variant to compound name
                    vmap = {
                        "streptonigrin": "Streptonigrin",
                        "gsk484": "GSK484", "gsk-484": "GSK484",
                        "gsk199": "GSK199", "gsk-199": "GSK199",
                        "cl-amidine": "Cl-amidine", "clamidine": "Cl-amidine",
                        "cl_amidine": "Cl-amidine",
                        "o-f-amidine": "o-F-Amidine",
                        "f-amidine": "F-Amidine",
                        "bb-cl": "BB-Cl-Amidine",
                        "amodiaquine": "Amodiaquine",
                        "pyroxamide": "Pyroxamide",
                        "bms-p5": "BMS-P5",
                        "jbi-589": "JBI-589", "jbi589": "JBI-589",
                        "tdfa": "TDFA", "pad-pf1": "PAD-PF1",
                    }
                    matched_name = vmap.get(variant)
                    if matched_name:
                        method = f"HIT_NAME({mol_name[:30]})"
                        break

        # (c) SMILES InChIKey match
        if matched_name is None and smiles and smiles != "nan":
            ik = std_ik_safe(smiles)
            if ik and ik in ik_to_name:
                matched_name = ik_to_name[ik]
                method = "HIT_IK"

        if matched_name is None:
            continue

        # Record first meaningful hit (prefer quantitative data)
        std_type  = str(row.get("Standard Type", ""))
        std_val   = str(row.get("Standard Value", ""))
        std_units = str(row.get("Standard Units", ""))
        pchembl   = str(row.get("pChEMBL Value", ""))

        if not results[matched_name]:
            results[matched_name] = {
                "hit": True,
                "molecule_chembl_id": mol_id,
                "molecule_name":      str(row.get("Molecule Name", "")).strip(),
                "standard_type":      std_type,
                "standard_value":     std_val,
                "standard_units":     std_units,
                "pchembl_value":      None if pchembl in ("nan", "", "None") else pchembl,
                "match_method":       method,
            }

    for name, h in results.items():
        if h:
            print(f"  {name}: FOUND — {h['molecule_chembl_id']} {h['standard_type']}={h['standard_value']}{h['standard_units']} method={h['match_method']}")
        else:
            print(f"  {name}: NOT FOUND in ChEMBL")
    return results


# =============================================================================
# Step 3 — BindingDB search
# =============================================================================

def search_bindingdb(refs: list, ref_iks: dict) -> dict[str, dict]:
    print("\n=== SEARCH 3: BindingDB ===")
    bdb = pd.read_csv(BDB_TSV, sep="\t", dtype=str, low_memory=False)

    ik_to_name: dict[str, str] = {
        info["effective_ik"]: name
        for name in (r[0] for r in refs)
        for info in [ref_iks[name]]
        if info["effective_ik"]
    }

    smi_col  = "Ligand SMILES"
    name_col = "Ligand HET ID (if a PDB ligand)" if "Ligand HET ID (if a PDB ligand)" in bdb.columns else None
    # Prefer a proper name column
    for c in ["BindingDB Ligand Name", "Ligand Aliases", "Ligand Name"]:
        if c in bdb.columns:
            name_col = c
            break

    results: dict[str, dict] = {name: {} for name, *_ in refs}

    for _, row in bdb.iterrows():
        matched_name = None
        method = None
        smi = str(row.get(smi_col, "")).strip()

        # (a) SMILES → InChIKey
        if smi and smi != "nan":
            ik = std_ik_safe(smi)
            if ik and ik in ik_to_name:
                matched_name = ik_to_name[ik]
                method = "HIT_IK"

        # (b) Name variant match
        if matched_name is None and name_col:
            lig_name = str(row.get(name_col, "")).strip().lower()
            for variant in NAME_VARIANTS:
                if variant in lig_name:
                    vmap = {
                        "streptonigrin": "Streptonigrin",
                        "gsk484": "GSK484", "gsk-484": "GSK484",
                        "gsk199": "GSK199", "gsk-199": "GSK199",
                        "cl-amidine": "Cl-amidine", "clamidine": "Cl-amidine",
                        "o-f-amidine": "o-F-Amidine",
                        "f-amidine": "F-Amidine",
                        "bb-cl": "BB-Cl-Amidine",
                        "amodiaquine": "Amodiaquine",
                        "pyroxamide": "Pyroxamide",
                        "bms-p5": "BMS-P5",
                        "jbi-589": "JBI-589", "jbi589": "JBI-589",
                        "tdfa": "TDFA", "pad-pf1": "PAD-PF1",
                    }
                    matched_name = vmap.get(variant)
                    if matched_name:
                        method = f"HIT_NAME({lig_name[:30]})"
                        break

        if matched_name is None:
            continue

        if not results[matched_name]:
            results[matched_name] = {
                "hit": True,
                "monomer_id":   str(row.get("BindingDB Reactant_set_id", "")),
                "ligand_name":  str(row.get(name_col, "")) if name_col else "",
                "IC50_nM":      str(row.get("IC50 (nM)", "")),
                "Ki_nM":        str(row.get("Ki (nM)", "")),
                "Kd_nM":        str(row.get("Kd (nM)", "")),
                "match_method": method,
            }

    for name, h in results.items():
        if h:
            print(f"  {name}: FOUND — IC50={h['IC50_nM']} Ki={h['Ki_nM']} method={h['match_method']}")
        else:
            print(f"  {name}: NOT FOUND in BindingDB")
    return results


# =============================================================================
# Step 4 — Pipeline trace
# =============================================================================

def trace_pipeline(refs: list, ref_iks: dict) -> dict[str, dict]:
    print("\n=== SEARCH 4: Pipeline trace ===")
    std_df  = pd.read_parquet(STD_PARQ,  columns=["inchi_key", "std_status", "source", "aid"])
    norm_df = pd.read_parquet(NORM_PARQ, columns=[
        "inchi_key", "source", "aid", "layer", "endpoint_type",
        "value_nM", "norm_status", "use_in_potency_model",
    ])
    comp_df = pd.read_parquet(COMP_PARQ, columns=["inchi_key", "pic50_consensus", "source_list"])

    std_iks  = set(std_df["inchi_key"].dropna())
    norm_iks = set(norm_df["inchi_key"].dropna())
    comp_iks = set(comp_df["inchi_key"].dropna())

    results = {}
    for name, manual_ik, cid, smiles, priority, notes in refs:
        info = ref_iks[name]
        ik   = info["effective_ik"]  # salt-stripped computed IK (best for lookup)
        also_manual = (manual_ik and manual_ik != ik)

        in_std  = (ik in std_iks)  or (also_manual and manual_ik in std_iks)
        in_norm = (ik in norm_iks) or (also_manual and manual_ik in norm_iks)
        in_comp = (ik in comp_iks) or (also_manual and manual_ik in comp_iks)

        # Effective IK that hit
        eff = ik
        if not in_std and also_manual and manual_ik in std_iks:
            eff = manual_ik

        pic50 = None
        if in_comp:
            row = comp_df[comp_df["inchi_key"].isin([ik, manual_ik] if also_manual else [ik])].iloc[0]
            pic50 = float(row["pic50_consensus"])

        exclusion = None
        if in_norm and not in_comp:
            sub = norm_df[norm_df["inchi_key"] == eff]
            # Trace why it didn't make it to pad4_compounds
            if (sub["use_in_potency_model"] == False).all():
                ep_types   = sub["endpoint_type"].unique().tolist()
                nst_types  = sub["norm_status"].unique().tolist()
                layers     = sub["layer"].unique().tolist()
                if all(ep in ("Kinact_Ki", "kon", "k_off", "Kcat", "Km", "Kcat/Km", "Unknown")
                       for ep in ep_types):
                    exclusion = f"endpoint_type not IC50 ({ep_types}) → use_in_potency_model=False"
                elif all(n in ("UNCONVERTIBLE_UNITS", "NO_VALUE") for n in nst_types):
                    exclusion = f"norm_status={nst_types} → value_nM=null → use_in_potency_model=False"
                elif all(l in ("D", "E", "HTS") for l in layers if l):
                    exclusion = f"layer={layers} (secondary/HTS) → use_in_potency_model=False"
                else:
                    exclusion = f"use_in_potency_model=False (ep={ep_types}, norm={nst_types}, layer={layers})"

        results[name] = {
            "found_in_standardized":  in_std,
            "found_in_normalized":    in_norm,
            "found_in_pad4_compounds": in_comp,
            "pic50_in_db":            pic50,
            "exclusion_reason":       exclusion,
        }

        marker = "✓ pad4_compounds" if in_comp else ("↳ std+norm only" if in_norm else ("↳ std only" if in_std else "✗ not in pipeline"))
        print(f"  {name}: {marker}", end="")
        if pic50:
            print(f"  pIC50={pic50:.3f}", end="")
        if exclusion:
            print(f"  [{exclusion}]", end="")
        print()

    return results


# =============================================================================
# Classification
# =============================================================================

def classify(name: str, priority: str,
             pubchem_hits: list, chembl_hit: dict, bdb_hit: dict,
             trace: dict) -> str:
    if priority == "EXCLUDE":
        return "excluded_correct"
    if trace["found_in_pad4_compounds"]:
        return "already_present"
    found_in_any_raw = bool(pubchem_hits) or bool(chembl_hit) or bool(bdb_hit)
    if not found_in_any_raw and not trace["found_in_standardized"]:
        return "absent_by_design"
    if trace["found_in_standardized"] and not trace["found_in_pad4_compounds"]:
        if trace["exclusion_reason"]:
            return "present_but_not_mapped"
        return "present_but_not_mapped"
    if found_in_any_raw and not trace["found_in_standardized"]:
        return "recoverable"
    return "present_but_not_mapped"


# =============================================================================
# Main
# =============================================================================

def main():
    print("PAD4-DB v2 — A2 Reference Compound Recovery Audit")
    print("=" * 60)

    # Step 0: standardize ref InChIKeys
    ref_iks = precompute_ref_inchikeys(REF_COMPOUNDS)

    # Step 1: PubChem
    pubchem_hits = search_pubchem_raw(REF_COMPOUNDS, ref_iks)

    # Step 2: ChEMBL
    chembl_hits = search_chembl(REF_COMPOUNDS, ref_iks)

    # Step 3: BindingDB
    bdb_hits = search_bindingdb(REF_COMPOUNDS, ref_iks)

    # Step 4: Pipeline trace
    traces = trace_pipeline(REF_COMPOUNDS, ref_iks)

    # ── Build results ─────────────────────────────────────────────────────
    all_results = []
    for name, manual_ik, cid, smiles, priority, notes in REF_COMPOUNDS:
        info = ref_iks[name]
        ph   = pubchem_hits.get(name, [])
        ch   = chembl_hits.get(name, {})
        bh   = bdb_hits.get(name, {})
        tr   = traces.get(name, {})
        rec  = classify(name, priority, ph, ch, bh, tr)

        all_results.append({
            "name":                   name,
            "priority":               priority,
            "notes":                  notes,
            "pubchem_cid":            cid,
            "inchikey_manual":        manual_ik,
            "inchikey_computed":      info["computed_ik"],
            "inchikey_match":         info["ik_match"],
            "effective_inchikey":     info["effective_ik"],
            "pubchem_hits":           ph[:5],  # cap at 5 for JSON readability
            "pubchem_n_raw_rows":     len(ph),
            "pubchem_aids":           sorted({x["aid"] for x in ph}),
            "chembl_hit":             bool(ch),
            "chembl_value":           f"{ch.get('standard_value','')} {ch.get('standard_units','')}" if ch else None,
            "chembl_pchembl":         ch.get("pchembl_value") if ch else None,
            "chembl_match_method":    ch.get("match_method") if ch else None,
            "bindingdb_hit":          bool(bh),
            "bindingdb_value":        f"IC50={bh.get('IC50_nM','')} Ki={bh.get('Ki_nM','')}" if bh else None,
            "bindingdb_match_method": bh.get("match_method") if bh else None,
            "found_in_standardized":  tr.get("found_in_standardized", False),
            "found_in_normalized":    tr.get("found_in_normalized", False),
            "found_in_pad4_compounds": tr.get("found_in_pad4_compounds", False),
            "pic50_in_db":            tr.get("pic50_in_db"),
            "exclusion_reason":       tr.get("exclusion_reason"),
            "recommendation":         rec,
        })

    # ── Print grouped report ──────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("RESULTS GROUPED BY PRIORITY")
    print("=" * 60)
    for prio in ("HIGH", "MEDIUM", "LOW", "EXCLUDE"):
        group = [r for r in all_results if r["priority"] == prio]
        if not group:
            continue
        found = sum(1 for r in group if r["found_in_pad4_compounds"])
        print(f"\n── {prio} priority ({found}/{len(group)} in pad4_compounds) ──")
        for r in group:
            status_icon = {
                "already_present":        "✓",
                "present_but_not_mapped": "↳",
                "recoverable":            "~",
                "absent_by_design":       "✗",
                "excluded_correct":       "⊘",
            }.get(r["recommendation"], "?")
            print(f"  [{status_icon}] {r['name']}")
            print(f"      IK computed: {r['inchikey_computed']}  (manual={'MATCH' if r['inchikey_match'] else ('MISMATCH' if r['inchikey_match'] is False else 'n/a')})")
            if r["found_in_pad4_compounds"]:
                print(f"      pIC50={r['pic50_in_db']:.3f} in pad4_compounds")
            if r["exclusion_reason"]:
                print(f"      EXCLUDED: {r['exclusion_reason']}")
            src_hits = []
            if r["pubchem_n_raw_rows"]: src_hits.append(f"PubChem ({r['pubchem_n_raw_rows']} rows, AIDs={r['pubchem_aids']})")
            if r["chembl_hit"]:         src_hits.append(f"ChEMBL ({r['chembl_value']}, {r['chembl_match_method']})")
            if r["bindingdb_hit"]:      src_hits.append(f"BindingDB ({r['bindingdb_value']}, {r['bindingdb_match_method']})")
            if src_hits:
                print(f"      Raw hits: {' | '.join(src_hits)}")
            else:
                print(f"      Raw hits: NONE")
            print(f"      Recommendation: {r['recommendation']}")

    # ── Summary counts ────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for prio in ("HIGH", "MEDIUM", "LOW"):
        group = [r for r in all_results if r["priority"] == prio]
        found = sum(1 for r in group if r["found_in_pad4_compounds"])
        in_std = sum(1 for r in group if r["found_in_standardized"])
        print(f"  {prio}: {found}/{len(group)} in pad4_compounds  {in_std}/{len(group)} in standardized")

    by_rec = {}
    for r in all_results:
        by_rec.setdefault(r["recommendation"], []).append(r["name"])
    print("\n  By recommendation:")
    for rec, names in sorted(by_rec.items()):
        print(f"    {rec}: {names}")

    # ── Write JSON ────────────────────────────────────────────────────────
    OUT_JSON.write_text(json.dumps(all_results, indent=2, default=str))
    print(f"\nWritten → {OUT_JSON}")


if __name__ == "__main__":
    main()
