# BindingDB Compliance Assessment — PAD4-DB v1.0
**Date:** 2026-06-22  
**Auditor role:** Licensing compliance — no science modifications

---

## 1. BindingDB License

BindingDB distributes its data under **Creative Commons Attribution 4.0 International (CC-BY 4.0)**. This license:
- Permits redistribution (in original or adapted form)
- Permits commercial use
- **Does NOT impose ShareAlike**
- **Requires attribution** (credit, copyright notice, license link, indication of changes)

Source: https://www.bindingdb.org/bind/info.jsp (BindingDB Terms of Use, as of 2024)

---

## 2. Raw Bulk Dump Check

**Finding: NO raw BindingDB bulk dump is redistributed.**

The BindingDB download (UniProt Q9UM07, downloaded 2026-06-10, 3,087 rows) is **not included in the deposition package**. Confirmed by examining all files in `deposition/zenodo/`:

| File | Raw BindingDB Records? | Finding |
|------|------------------------|---------|
| `standardized_structures.sdf` | No | Structures transformed; no BindingDB IDs; no raw Ki/Kd/IC50 values from BindingDB |
| `source_independence_scores.csv` | No | "bindingdb" in source_list is a provenance label, not a BindingDB record |
| `activity_cliffs.csv` | No | Derived pairs; no raw measurements |
| `cliff_hubs.csv` | No | 4 standardized SMILES; source_list references BindingDB as provenance |
| All other CSVs | No | Derived statistics only |

---

## 3. Attribution Requirement

CC-BY 4.0 requires attribution with:
1. Credit to BindingDB creators/source
2. URI or hyperlink to license
3. Indication that material has been modified (if applicable)

**Current README_zenodo.md attribution:**
> "Source bioactivity data derive from PubChem BioAssay, ChEMBL (assay CHEMBL6111), and BindingDB (UniProt Q9UM07)"

**Gaps identified:**

| Requirement | Current state | Gap |
|-------------|--------------|-----|
| Credit BindingDB | ✅ Present | — |
| UniProt accession cited | ✅ Present | — |
| Download date | ❌ Missing | "downloaded 2026-06-10" not in README_zenodo.md (present in manuscript Methods only) |
| BindingDB license statement | ❌ Missing | CC-BY 4.0 not named in README |
| BindingDB URL | ❌ Missing | https://www.bindingdb.org not included |
| Indication of modification | ❌ Missing | Structures were standardized (RDKit); this should be stated |

**Required addition (see README_release_revision.md):**
> "BindingDB data (UniProt Q9UM07) were downloaded 2026-06-10 from BindingDB (https://www.bindingdb.org), licensed under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/). Structures were standardized using RDKit 2025.09.5."

---

## 4. Summary

| Check | Status |
|-------|--------|
| Raw BindingDB bulk dump in package? | ✅ NO |
| BindingDB credited? | ✅ YES (partial) |
| CC-BY 4.0 license named? | ❌ MISSING |
| Download date present? | ❌ MISSING from README_zenodo.md |
| BindingDB URL present? | ❌ MISSING |
| Modification noted (standardization)? | ❌ MISSING |
| ShareAlike triggered? | ✅ NO (CC-BY 4.0 has no SA clause) |

**Verdict:** Attribution is partially present but incomplete. No ShareAlike obligation. Fix README_zenodo.md before deposition. See README_release_revision.md for corrected text.
