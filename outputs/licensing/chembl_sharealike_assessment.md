# ChEMBL ShareAlike Assessment — PAD4-DB v1.0
**Date:** 2026-06-22  
**Auditor role:** Licensing compliance — no science modifications

---

## 1. ChEMBL License Summary

ChEMBL is distributed by EMBL-EBI under **Creative Commons Attribution-ShareAlike 4.0 International (CC-BY-SA 4.0)** as of ChEMBL 29 (2021-02-01). This applies to all subsequent releases including the version downloaded for PAD4-DB (2026-06-14).

**CC-BY-SA 4.0 ShareAlike clause (relevant section):**  
> "If You Share Adapted Material You produce, the following conditions also apply:
> The Adapter's License You apply must be a Creative Commons license with the same License Elements..."

**UNCERTAINTY:** Whether "Adapted Material" under CC-BY-SA 4.0 applies to (a) individual chemical structures extracted from ChEMBL, or (b) only to databases derived from ChEMBL as a whole, is not definitively established in case law. This assessment applies a conservative reading.

---

## 2. ChEMBL-Sourced Records in the Deposition Package

| File | Type of ChEMBL Content | Assessment |
|------|------------------------|------------|
| `standardized_structures.sdf` | 1,566 structures (out of 3,093) with source_list containing "chembl" | **TRANSFORMED** — RDKit-standardized SMILES, new RDKit-generated 2D coordinates, canonicalized InChIKey as identifier. No ChEMBL IDs retained. |
| `source_independence_scores.csv` | "chembl" appears in source_list column as provenance string | **ATTRIBUTION ONLY** — not a ChEMBL record; records traced back to ChEMBL for provenance |
| `activity_cliffs.csv` | "chembl" in source_combination string | **ATTRIBUTION ONLY** — derived output; no ChEMBL records |
| `cliff_hubs.csv` | 2 of 4 hub compounds have "chembl" in source_list | **TRANSFORMED** — standardized SMILES for 2 compounds |
| All other CSVs | Derived statistics only | **NOT APPLICABLE** — no ChEMBL structures or records |

---

## 3. Assessment by Category

### 3a. Direct Redistribution
**Finding: NO direct redistribution detected.**  
No raw ChEMBL records, SDF exports, or bulk activity tables are included. No ChEMBL IDs (CHEMBL\d+), assay IDs, or SID/CID identifiers appear in any released file.

### 3b. Transformed Redistribution (⚠ UNCERTAIN)
**Finding: PRESENT — conservatively flagged as potentially triggering ShareAlike.**

The `standardized_structures.sdf` contains 3,093 molecular structures. Approximately 1,566 of these have `source_list` containing "chembl" in `source_independence_scores.csv`, meaning the standardized structure traces back in part to a ChEMBL-sourced SMILES. The structures have been:
- SMILES-standardized by RDKit (desalting, canonicalization, charge normalization)
- Assigned new RDKit-generated 2D coordinates (not from ChEMBL)
- Named by InChIKey (computed, not a ChEMBL identifier)

**Conservative position:** A European Database Directive reading could treat the 1,566 ChEMBL-traced structures as a substantial extraction from ChEMBL (≈ all 1,566 records from CHEMBL6111, which had 4,925 rows before dedup). The transformation (RDKit standardization) is unlikely to be considered "sufficient" under database right doctrine to break the derivative chain.

**Permissive position:** In US copyright law, factual data (chemical structures) generally lack copyright protection unless there is sufficient originality in selection/arrangement. A curated subset of 3,093 compounds from three sources is an original selection; the individual structures are not copyrightable. Under this reading, no ShareAlike obligation arises.

**Recommendation:** Apply the **conservative position** and license the dataset as CC-BY-SA 4.0 to avoid risk.

### 3c. Derived Analytical Output (not subject to ShareAlike)
All other files (null_model_comparison.csv, scaffold_ruggedness_table.csv, transformation_impact_table.csv, hub metrics, permutation statistics) contain **only derived statistics** — aggregated counts, means, ranks, p-values. These do not constitute "Adapted Material" under CC-BY-SA 4.0.

---

## 4. Attribution Requirement (ChEMBL CC-BY-SA 4.0)

Even under the permissive position, ChEMBL requires attribution. Current README_zenodo.md contains:

> "Source bioactivity data derive from PubChem BioAssay, ChEMBL (assay CHEMBL6111), and BindingDB (UniProt Q9UM07)"

**Gap:** This statement does not:
- Name the ChEMBL release version (e.g., ChEMBL 34)
- Include the ChEMBL license statement ("CC-BY-SA 4.0")
- Provide the EMBL-EBI URL or DOI for ChEMBL

**Required addition (see README_release_revision.md for full text):**
> "ChEMBL data were obtained from ChEMBL release [version] (EMBL-EBI, https://www.ebi.ac.uk/chembl/), licensed under CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/)."

---

## 5. Summary

| Question | Finding |
|----------|---------|
| Raw ChEMBL records redistributed? | NO |
| Transformed ChEMBL structures redistributed? | YES — 1,566 structures in SDF, 2 in cliff_hubs.csv |
| ShareAlike triggered (conservative)? | **UNCERTAIN — YES under conservative/EU reading** |
| ShareAlike triggered (permissive)? | NO under US copyright reading |
| Attribution present? | PARTIAL — missing version, license name, URL |
| Recommended action | Change license to CC-BY-SA 4.0; add ChEMBL attribution with version and license URL |

**FLAG:** This assessment is not legal advice. The ShareAlike question for chemical structure databases is not settled. If the authors have institutional legal counsel, this specific question should be reviewed before deposition.
