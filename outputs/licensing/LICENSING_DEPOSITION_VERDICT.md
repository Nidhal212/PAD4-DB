# Licensing and Deposition Compliance Audit — Final Verdict
**PAD4-DB v1.0**  
**Date:** 2026-06-22  
**Scope:** Licensing, attribution, redistribution rights, Zenodo compliance only. No science changes.

---

## Verdict

```
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║   READY AFTER README / LICENSE UPDATE                                        ║
║                                                                              ║
║   DO NOT DEPOSIT WITH CURRENT CC-BY-4.0 LICENSE.                            ║
║   Required changes are administrative (no data changes, no science changes)  ║
║   and can be completed in approximately 1–2 hours.                           ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
```

---

## Rationale

### What is correct ✅

- No raw source records from any database are redistributed
- All source databases are credited in README_zenodo.md
- The Zenodo package contains no manuscript files (correctly held back)
- All 9 data files and 7 code files are appropriate for open deposition
- The SDF and CSVs contain only derived/transformed content; no raw activity tables
- No source database identifiers (ChEMBL IDs, PubChem SIDs/CIDs, BindingDB IDs) are present
- Package structure is coherent and complete (README_data.md, README_code.md, environment.yml)

### What must be fixed before deposition

#### BLOCKER 1 — License field is incorrect
**Current:** `"license": "CC-BY-4.0"` in DEPOSITION_METADATA.json  
**Required:** `"license": "CC-BY-SA-4.0"`  
**Reason:** ChEMBL (source of ~1,566 structures in the SDF) is licensed CC-BY-SA 4.0. Including ChEMBL-derived structures in a released database and licensing it as CC-BY-4.0 conflicts with ChEMBL's ShareAlike clause under a conservative reading. This is the primary license risk.  
**Uncertainty:** Whether individual chemical structures trigger database ShareAlike is not definitively settled in law. The conservative position (CC-BY-SA) is recommended. If institutional legal counsel is available, confirm before deposition.  
**File to change:** `PAD4_DB_FREEZE_v1.0/deposition/metadata/DEPOSITION_METADATA.json`

#### BLOCKER 2 — README_zenodo.md license statement is incorrect
**Current:** "CC-BY-4.0"  
**Required:** "CC-BY-SA 4.0" + per-source attribution block  
**File to change:** `PAD4_DB_FREEZE_v1.0/deposition/zenodo/README_zenodo.md`  
**Full corrected text:** see `outputs/licensing/README_release_revision.md`

#### BLOCKER 3 — Author names are placeholders
**Current:** "[Author 1, Last, First]" × 3 in DEPOSITION_METADATA.json  
**Required:** Real names, affiliations, ORCIDs  
**This is a metadata completeness requirement, not a licensing issue, but blocks deposition.**

### What should be fixed (HIGH, not blocking deposition)

- ChEMBL release version not named (add to README)
- ChEMBL CC-BY-SA 4.0 not cited by name or URL (add to README)
- BindingDB CC-BY 4.0 not cited by name or URL (add to README)
- Download dates not in README_zenodo.md (add — already in manuscript Methods)
- RDKit standardization not disclosed in README (add — important for modification notice)
- GitHub URL missing from README_zenodo.md and DEPOSITION_METADATA.json related_identifiers
- Funding field missing from DEPOSITION_METADATA.json
- Publication_date should be updated to actual upload date (currently 2026-06-21, a placeholder)

---

## Ordered Unblocking Sequence

**Step 1 (~10 min):** Update `DEPOSITION_METADATA.json`:
- Change `"license"` to `"CC-BY-SA-4.0"`
- Add GitHub to `related_identifiers`
- Add funding field
- Update publication_date at actual upload time

**Step 2 (~30 min):** Replace `README_zenodo.md` with corrected version from `outputs/licensing/README_release_revision.md`
- Fill ChEMBL release version (check download metadata)
- Fill contact name and email

**Step 3 (~30 min):** Fill author metadata (names, affiliations, ORCIDs) in DEPOSITION_METADATA.json

**Step 4:** Proceed to Zenodo deposition once author metadata + DOI is available

---

## Files Produced by This Audit

| File | Purpose |
|------|---------|
| `outputs/licensing/licensing_inventory.csv` | Per-file provenance and content classification |
| `outputs/licensing/chembl_sharealike_assessment.md` | ChEMBL CC-BY-SA analysis |
| `outputs/licensing/bindingdb_compliance_assessment.md` | BindingDB CC-BY 4.0 compliance |
| `outputs/licensing/pubchem_compliance_assessment.md` | PubChem public domain verification |
| `outputs/licensing/recommended_license_decision.md` | Option A/B/C analysis; recommendation = CC-BY-SA 4.0 |
| `outputs/licensing/README_release_revision.md` | Full corrected README_zenodo.md text |
| `outputs/licensing/zenodo_metadata_checklist.md` | Field-by-field metadata audit + corrected JSON |
| `outputs/licensing/LICENSING_DEPOSITION_VERDICT.md` | This file |

---

## Science Freeze Confirmed

No scientific results, numbers, figures, or manuscript content were modified or reviewed in this audit. All canonical statistics remain locked as per CLAUDE.md.
