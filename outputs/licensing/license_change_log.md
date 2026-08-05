# License Change Log — PAD4-DB v1.0
**Date:** 2026-06-22  
**Change:** CC-BY-4.0 → CC-BY-SA-4.0  
**Reason:** ChEMBL (source of ~1,566 structures in standardized_structures.sdf) is licensed CC-BY-SA 4.0. Conservative compliance with ChEMBL's ShareAlike clause requires the derivative dataset to carry the same license.

---

## Files Modified

### 1. `PAD4_DB_FREEZE_v1.0/deposition/metadata/DEPOSITION_METADATA.json`

| Field | Before | After |
|-------|--------|-------|
| `license` | `"CC-BY-4.0"` | `"CC-BY-SA-4.0"` |
| `title` | "PAD4-DB: A Curated..." | "PAD4-DB v1.0: A Curated..." (version added) |
| `description` | 150 words, no license info | Updated to include source license summary |
| `keywords` | 9 keywords | 15 keywords (added: PAD4 inhibitors, bioactivity database, RDKit, PubChem, ChEMBL, BindingDB) |
| `related_identifiers` | 1 entry (manuscript DOI placeholder) | 2 entries (+ GitHub repo) |
| `funding` | Not present | Added as placeholder field |
| `publication_date` | "2026-06-21" (stale) | Replaced with instruction to set at upload time |

### 2. `PAD4_DB_FREEZE_v1.0/deposition/zenodo/README_zenodo.md`

| Section | Before | After |
|---------|--------|-------|
| License statement | "CC-BY-4.0" | "CC BY-SA 4.0" with full URL |
| Header | Title only | Added: version, license, contact placeholder, GitHub URL, manuscript DOI placeholder |
| Provenance | 1-line sentence | Full per-source block: PubChem (public domain + URL + date), ChEMBL (CC-BY-SA + URL + date + version placeholder + attribution DOI + ShareAlike notice), BindingDB (CC-BY + URL + date + attribution DOI) |
| Software provenance | Not present | Added: Python 3.10.19, RDKit 2025.09.5, pandas 2.3.3, NumPy 2.2.5, SciPy 1.15.3, matplotlib 3.10.8, networkx 3.4.2 |
| "Package contents" note | "release candidate" language | Cleaned to describe final package; held-back note retained |

---

## Files NOT Modified

| File | Reason |
|------|--------|
| All CSV data files | No license changes touch data |
| standardized_structures.sdf | No changes to data |
| All Python scripts | Code is author copyright, not subject to source database licenses |
| environment.yml | No changes required |
| README_data.md | No license statements present; data attribution in README_zenodo.md |
| README_code.md | No license statements present |

---

## LICENSE File — Status

No standalone `LICENSE` file was present in the Zenodo deposition package (`deposition/zenodo/`). The license is declared in:
- `DEPOSITION_METADATA.json` (machine-readable, Zenodo API field)
- `README_zenodo.md` (human-readable, with full URL)

**Recommendation:** Consider adding a `LICENSE` file containing the CC-BY-SA 4.0 full text or a standard SPDX reference. This is not required by Zenodo but is good practice. Standard content:

```
Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)
https://creativecommons.org/licenses/by-sa/4.0/

This dataset is licensed under CC BY-SA 4.0.
See https://creativecommons.org/licenses/by-sa/4.0/legalcode for the full license text.

Source data attributions:
- PubChem BioAssay (NCBI/NLM/NIH): public domain
- ChEMBL (EMBL-EBI): CC BY-SA 4.0
- BindingDB: CC BY 4.0
```

---

## Verification

```
JSON validity:        ✅ PASS (python3 json.load confirmed)
license field:        ✅ CC-BY-SA-4.0
GitHub in metadata:   ✅ PRESENT
README license line:  ✅ CC BY-SA 4.0 with URL
ChEMBL attribution:   ✅ PRESENT (version placeholder remains — author action required)
BindingDB attribution:✅ PRESENT
PubChem attribution:  ✅ PRESENT
Software versions:    ✅ PRESENT (all 7 packages)
```
