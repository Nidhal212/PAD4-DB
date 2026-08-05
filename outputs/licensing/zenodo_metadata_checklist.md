# Zenodo Metadata Checklist — PAD4-DB v1.0
**Date:** 2026-06-22  
**Source file audited:** `PAD4_DB_FREEZE_v1.0/deposition/metadata/DEPOSITION_METADATA.json`

---

## Field-by-Field Audit

### Creators / Authors
| Field | Current value | Status | Required action |
|-------|--------------|--------|-----------------|
| authors[0].name | "[Author 1, Last, First]" | ❌ PLACEHOLDER | Fill with real name |
| authors[0].affiliation | "[Affiliation 1]" | ❌ PLACEHOLDER | Fill with institution |
| authors[0].orcid | "[ORCID]" | ❌ PLACEHOLDER | Fill or remove if no ORCID |
| authors[1].name | "[Author 2, Last, First]" | ❌ PLACEHOLDER | Fill with real name |
| authors[1].affiliation | "[Affiliation 2]" | ❌ PLACEHOLDER | Fill with institution |
| authors[1].orcid | "[ORCID]" | ❌ PLACEHOLDER | Fill or remove if no ORCID |
| authors[2].name | "[Corresponding Author, Last, First]" | ❌ PLACEHOLDER | Fill with real name |
| authors[2].affiliation | "[Affiliation 1]" | ❌ PLACEHOLDER | Fill with institution |
| authors[2].orcid | "[ORCID]" | ❌ PLACEHOLDER | Fill or remove if no ORCID |

### Title
| Field | Current value | Status |
|-------|--------------|--------|
| title | "PAD4-DB: A Curated Structure-Activity Resource for PAD4 Inhibitors" | ✅ PRESENT |

Note: Paper title includes "hub-organized activity cliffs and scaffold-dependent SAR ruggedness" — consider whether Zenodo title should match manuscript title exactly for discoverability.

### Description
| Field | Current value | Status |
|-------|--------------|--------|
| description | Present, ~170 words | ✅ ADEQUATE |

Gap: description does not mention CC-BY-SA 4.0 or ChEMBL ShareAlike. Minor — not required in description field.

### Keywords
| Current keywords | Status |
|-----------------|--------|
| PAD4 | ✅ |
| PADI4 | ✅ |
| protein arginine deiminase | ✅ |
| activity cliffs | ✅ |
| structure-activity relationships | ✅ |
| cheminformatics | ✅ |
| matched molecular pairs | ✅ |
| SAR ruggedness | ✅ |
| drug discovery | ✅ |

**Missing suggested keywords:**
- `bioactivity database`
- `PAD4 inhibitors`
- `RDKit`
- `PubChem`
- `ChEMBL`
- `BindingDB`

Zenodo allows up to ~30 keywords. Adding source database names improves discoverability.

### License
| Field | Current value | Status | Required action |
|-------|--------------|--------|-----------------|
| license | "CC-BY-4.0" | ❌ INCORRECT | Change to "CC-BY-SA-4.0" (see recommended_license_decision.md) |

Zenodo's license identifier for CC-BY-SA 4.0: `CC-BY-SA-4.0`

### Access Right
| Field | Current value | Status |
|-------|--------------|--------|
| access_right | "open" | ✅ CORRECT |

### Version
| Field | Current value | Status |
|-------|--------------|--------|
| version | "1.0" | ✅ PRESENT |

### Upload Type
| Field | Current value | Status |
|-------|--------------|--------|
| upload_type | "dataset" | ✅ CORRECT |

### Publication Date
| Field | Current value | Status | Note |
|-------|--------------|--------|------|
| publication_date | "2026-06-21" | ⚠ VERIFY | Should match actual Zenodo upload date; update if deposition is later |

### Funding
| Field | Current value | Status | Required action |
|-------|--------------|--------|-----------------|
| funding | **NOT PRESENT** | ❌ MISSING | Add funding field with grant agency and number once author metadata is provided |

Zenodo supports: `{"funder": {"name": "..."}, "award": {"number": "...", "title": "..."}}`

### Related Identifiers
| Field | Current value | Status | Required action |
|-------|--------------|--------|-----------------|
| related_identifiers[0].identifier | "[manuscript DOI - to be assigned]" | ❌ PLACEHOLDER | Fill once manuscript DOI is assigned |
| related_identifiers[0].relation | "isSupplementTo" | ✅ CORRECT relation type | — |
| GitHub repository | **NOT PRESENT** | ❌ MISSING | Add `{"identifier": "https://github.com/Nidhal212/PAD4-DB", "relation": "isSupplementTo", "resource_type": "software"}` |

### GitHub Linkage
| Check | Status |
|-------|--------|
| GitHub repo exists | ✅ YES — https://github.com/Nidhal212/PAD4-DB |
| GitHub URL in related_identifiers | ❌ MISSING — must be added |
| GitHub URL in README_zenodo.md | ❌ MISSING — must be added |

---

## Corrected DEPOSITION_METADATA.json

Replace the current file with:

```json
{
  "title": "PAD4-DB v1.0: A Curated Structure–Activity Resource for PAD4 Inhibitors",
  "description": "PAD4-DB is a standardized, deduplicated knowledge base of 3,093 structurally resolved PAD4 (PADI4) inhibitors with consensus pIC50 values (range 2.00–8.52, median 6.84), assembled from PubChem BioAssay, ChEMBL (assay CHEMBL6111), and BindingDB (UniProt Q9UM07) through a fully scripted curation pipeline. Includes Murcko scaffold annotations, per-compound source-independence scores, the severe activity-cliff network (94 severe pairs), four hub-compound labels (53.2% hub concentration), a per-scaffold ruggedness ranking (96.8% smooth series), and matched-molecular-pair transformation typology (80/94 severe cliffs MMP-confirmed). Structures standardized using RDKit 2025.09.5. Source databases: PubChem (public domain), ChEMBL CC BY-SA 4.0, BindingDB CC BY 4.0.",
  "upload_type": "dataset",
  "authors": [
    {"name": "[Author 1, Last, First]", "affiliation": "[Affiliation 1]", "orcid": "[ORCID or remove field]"},
    {"name": "[Author 2, Last, First]", "affiliation": "[Affiliation 2]", "orcid": "[ORCID or remove field]"},
    {"name": "[Corresponding Author, Last, First]", "affiliation": "[Affiliation 1]", "orcid": "[ORCID or remove field]"}
  ],
  "keywords": [
    "PAD4", "PADI4", "protein arginine deiminase", "PAD4 inhibitors",
    "activity cliffs", "structure-activity relationships", "SAR ruggedness",
    "cheminformatics", "matched molecular pairs", "bioactivity database",
    "RDKit", "PubChem", "ChEMBL", "BindingDB", "drug discovery"
  ],
  "license": "CC-BY-SA-4.0",
  "access_right": "open",
  "version": "1.0",
  "publication_date": "[date of actual Zenodo upload]",
  "funding": [
    {"funder": {"name": "[Grant agency]"}, "award": {"number": "[grant number]", "title": "[grant title if applicable]"}}
  ],
  "related_identifiers": [
    {
      "identifier": "[manuscript DOI - to be assigned]",
      "relation": "isSupplementTo",
      "resource_type": "publication-article"
    },
    {
      "identifier": "https://github.com/Nidhal212/PAD4-DB",
      "relation": "isSupplementTo",
      "resource_type": "software"
    }
  ],
  "_placeholders_to_fill_before_deposition": [
    "authors[].name, affiliation, orcid (or remove orcid field if not available)",
    "funding[].funder.name and award.number",
    "related_identifiers[0].identifier (manuscript DOI)",
    "publication_date (use actual upload date)"
  ]
}
```

---

## Checklist Summary

| Field | Status |
|-------|--------|
| Title | ✅ PRESENT |
| Description | ✅ ADEQUATE |
| Authors | ❌ 3× PLACEHOLDER |
| Keywords | ⚠ PARTIAL (missing 6 recommended) |
| License | ❌ INCORRECT — must change to CC-BY-SA-4.0 |
| Access right | ✅ CORRECT |
| Version | ✅ PRESENT |
| Upload type | ✅ CORRECT |
| Publication date | ⚠ VERIFY at upload time |
| Funding | ❌ MISSING |
| Related identifier: manuscript DOI | ❌ PLACEHOLDER |
| Related identifier: GitHub | ❌ MISSING |
