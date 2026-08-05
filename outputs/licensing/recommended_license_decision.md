# License Recommendation — PAD4-DB v1.0
**Date:** 2026-06-22  
**Auditor role:** Licensing compliance — no science modifications  
**Disclaimer:** This is not legal advice. Flag all uncertainties to institutional legal counsel if available.

---

## Options Evaluated

| Option | License | ShareAlike | Summary |
|--------|---------|------------|---------|
| A | CC-BY 4.0 | No | Current DEPOSITION_METADATA.json selection |
| B | CC-BY-SA 4.0 | Yes | Inherited from ChEMBL (conservative reading) |
| C | Mixed-license (structure data CC-BY-SA; derived stats CC-BY) | Partial | Maximum legal clarity; highest complexity |

---

## Analysis

### Option A — CC-BY 4.0

**Justification for:** PubChem is public domain. BindingDB is CC-BY 4.0. Molecular structures may not be individually copyrightable under US law. The transformations applied (RDKit standardization, new 2D coordinates, InChIKey naming) could be argued to create a sufficiently distinct work.

**Justification against:**
- ChEMBL is CC-BY-SA 4.0. Under a conservative (European Database Directive) reading, extracting and redistributing ~1,566 structures derived from CHEMBL6111 — even standardized — may constitute "adapted material" requiring ShareAlike.
- If challenged, the defense would require demonstrating that individual chemical structures are uncopyrightable facts — a position that is legally untested for databases in EU jurisdictions.
- EMBL-EBI has explicitly licensed ChEMBL under CC-BY-SA (not CC-BY), which signals intent to require ShareAlike.
- Using CC-BY when the source requires CC-BY-SA would be a license violation under the conservative reading.

**Risk level: MEDIUM-HIGH in EU jurisdictions. LOW in US-only context.**

### Option B — CC-BY-SA 4.0

**Justification for:**
- Inherits ChEMBL's own license — legally unambiguous.
- Compatible with BindingDB (CC-BY) and PubChem (public domain) in both directions.
- Downstream users receive clear notice that derivative databases must also be CC-BY-SA.
- Lowest legal risk.
- ChEMBL is itself CC-BY-SA; aligning with the most restrictive source is standard practice.
- No journal or Zenodo policy prevents CC-BY-SA for datasets.

**Justification against:**
- Requires downstream commercial database builders to open-source their derivatives, which may limit some industrial applications.
- Slightly more restrictive than CC-BY.

**Risk level: LOW.**

### Option C — Mixed license

**Justification for:**
- Theoretically maximizes freedom: derived-statistics-only files (null_model_comparison.csv, scaffold_ruggedness_table.csv, etc.) are clearly original work with no database-right entanglement and could be CC-BY.
- Structure files (standardized_structures.sdf) carry the ShareAlike obligation.

**Justification against:**
- Mixed licensing on a single Zenodo record requires splitting files into separate deposits with separate license fields, or clear per-file license headers — neither of which Zenodo easily supports in a single record.
- Significantly increases documentation burden.
- Most readers/users will not track per-file licenses.
- Scientific Data and similar journals prefer a single dataset license.

**Risk level: LOW (if done correctly). Implementation complexity: HIGH.**

---

## Recommendation

**OPTION B — CC-BY-SA 4.0**

This is the conservative and legally sound choice given that ChEMBL structures are included in the SDF. It:
1. Satisfies the ChEMBL license requirement under all readings
2. Is compatible with BindingDB (CC-BY ⊆ CC-BY-SA compatibility: both allow ShareAlike downstream)
3. Is compatible with PubChem (public domain ⊆ any CC license)
4. Does not conflict with Zenodo, Scientific Data, or JCIM deposition requirements
5. Is the license already used by ChEMBL itself — citing it avoids any impression of license laundering

**Required change:** Update `DEPOSITION_METADATA.json` field `"license"` from `"CC-BY-4.0"` to `"CC-BY-SA-4.0"`. Update `README_zenodo.md` license statement accordingly.

---

## Compatibility Matrix

| Source | Source License | Compatible with CC-BY-SA? |
|--------|---------------|--------------------------|
| PubChem | Public domain | ✅ Yes |
| BindingDB | CC-BY 4.0 | ✅ Yes |
| ChEMBL | CC-BY-SA 4.0 | ✅ Yes (same license) |
| Original code/scripts | Author copyright | ✅ Yes (authors can choose) |
| Derived statistics | Author copyright | ✅ Yes (authors can choose) |

---

## Required Follow-up Actions

1. Update `DEPOSITION_METADATA.json`: `"license": "CC-BY-SA-4.0"`
2. Update `README_zenodo.md`: change "CC-BY-4.0" to "CC-BY-SA 4.0" and add source attribution block
3. Add per-source attribution statements (see README_release_revision.md)
4. If institutional legal counsel is available: confirm this reading before public deposition
