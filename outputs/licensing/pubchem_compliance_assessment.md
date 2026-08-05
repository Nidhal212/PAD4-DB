# PubChem Attribution Assessment — PAD4-DB v1.0
**Date:** 2026-06-22  
**Auditor role:** Licensing compliance — no science modifications

---

## 1. PubChem License Status

PubChem is a product of the **National Center for Biotechnology Information (NCBI), National Library of Medicine (NLM), National Institutes of Health (NIH)**, a US federal government agency. As a US government work, PubChem data is in the **public domain** under 17 U.S.C. § 105:

> "Copyright protection under this title is not available for any work of the United States Government."

**There are no copyright restrictions on PubChem data. No license is required. Attribution is recommended but not legally mandated.**

---

## 2. Raw PubChem Bulk Dump Check

**Finding: NO raw PubChem bioassay data is redistributed.**

The raw PubChem data (97 files, 341,282 rows, downloaded 2026-06-10 to 2026-06-14) is **not included in the deposition package**. Confirmed by examining all files in `deposition/zenodo/`:

| File | Raw PubChem Records? | Finding |
|------|----------------------|---------|
| `standardized_structures.sdf` | No | Structures standardized; no SID/CID/AID identifiers present |
| `source_independence_scores.csv` | No | "pubchem_confirmatory" in source_list is a provenance label |
| All other CSVs | No | Derived statistics only |

No AID numbers, SID numbers, CID numbers, or raw bioassay outcome values appear in any released file.

---

## 3. Attribution Assessment

PubChem attribution is recommended practice even though not legally required. NLM/NCBI requests acknowledgment in publications.

**Current README_zenodo.md statement:**
> "Source bioactivity data derive from PubChem BioAssay, ChEMBL (assay CHEMBL6111), and BindingDB (UniProt Q9UM07)"

**Assessment by requirement:**

| Requirement | Status | Note |
|-------------|--------|------|
| PubChem named | ✅ YES | "PubChem BioAssay" present |
| Download date | ❌ MISSING from README | In manuscript Methods ("2026-06-10 to 2026-06-14"); not in README_zenodo.md |
| PubChem URL | ❌ MISSING | https://pubchem.ncbi.nlm.nih.gov not included |
| AID-level provenance | ✅ ADEQUATE | Not required in README; full AID list is in manuscript/audit files |
| Copyright statement | N/A | Public domain — no copyright notice required |

---

## 4. Provenance Statement Accuracy

**Verified accurate:**
- Data source: "PubChem BioAssay" ✅ (57 confirmatory AIDs + 11 Layer C + 26 Layer D/E + 3 HTS = 97 files)
- Assay categories: confirmatory, secondary, literature-derived, HTS ✅
- Download period: 2026-06-10 to 2026-06-14 ✅ (per DOWNLOAD_DATES.txt)

No inaccurate provenance statements detected.

---

## 5. Summary

| Check | Status |
|-------|--------|
| Raw PubChem bulk data redistributed? | ✅ NO |
| PubChem credited? | ✅ YES |
| Legal restriction on redistribution? | ✅ NONE (public domain) |
| Download date in README_zenodo.md? | ❌ MISSING |
| PubChem URL in README? | ❌ MISSING |
| Provenance statements accurate? | ✅ YES |

**Verdict:** No compliance risk. Public domain source. Attribution should be improved (download date + URL) for scholarly completeness. See README_release_revision.md for corrected text.
