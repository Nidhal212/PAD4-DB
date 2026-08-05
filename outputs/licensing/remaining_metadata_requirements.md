# Remaining Metadata Requirements — PAD4-DB v1.0
**Date:** 2026-06-22  
**Status after licensing corrections:** All licensing fields resolved. Remaining gaps are author-supplied metadata only.

---

## Author Metadata (BLOCKER for deposition)

These fields are in `PAD4_DB_FREEZE_v1.0/deposition/metadata/DEPOSITION_METADATA.json` and must be filled before Zenodo upload.

| Field | Location | Current value | Required |
|-------|----------|--------------|---------|
| Author 1 name | `authors[0].name` | "[Author 1, Last, First]" | Real name, format "Last, First" |
| Author 1 affiliation | `authors[0].affiliation` | "[Affiliation 1]" | Full institution name |
| Author 1 ORCID | `authors[0].orcid` | "[ORCID...]" | 16-digit ORCID or **remove field** |
| Author 2 name | `authors[1].name` | "[Author 2, Last, First]" | Real name |
| Author 2 affiliation | `authors[1].affiliation` | "[Affiliation 2]" | Full institution name |
| Author 2 ORCID | `authors[1].orcid` | "[ORCID...]" | 16-digit ORCID or **remove field** |
| Corresponding author name | `authors[2].name` | "[Corresponding Author, Last, First]" | Real name |
| Corresponding author affiliation | `authors[2].affiliation` | "[Affiliation 1]" | Full institution name |
| Corresponding author ORCID | `authors[2].orcid` | "[ORCID...]" | 16-digit ORCID or **remove field** |

**ORCID note:** If an author does not have an ORCID, remove the `"orcid"` key from that author object entirely. Do not leave a placeholder string — Zenodo will reject it.

---

## README Contact Block (BLOCKER for deposition)

In `PAD4_DB_FREEZE_v1.0/deposition/zenodo/README_zenodo.md`, the contact section reads:

```
[Corresponding author full name]
[Title, Department, Institution]
[email]
```

Replace with actual name, title, and email before upload.

---

## ChEMBL Release Version (HIGH — required for reproducibility)

In `README_zenodo.md`, the ChEMBL block contains:

```
Release version: [INSERT ChEMBL release number — check ebi.ac.uk/chembl/db_info/]
```

**How to find:** Go to https://www.ebi.ac.uk/chembl/ → "Release Notes" or check the download page. The release number appears as "ChEMBL XX" (e.g., ChEMBL 34 or ChEMBL 35). Cross-reference with download date of 2026-06-14.

The CHEMBL6111 CSV file (the only ChEMBL file downloaded) does not embed a release version number; it must be retrieved from the ChEMBL website or from the original download browser tab history.

---

## Funding (MEDIUM — Zenodo field; not blocking)

In `DEPOSITION_METADATA.json`:

```json
"funding": [
  {
    "funder": {"name": "[Grant agency]"},
    "award": {"number": "[grant number]"}
  }
]
```

Replace with actual grant agency and number, or remove the `funding` field entirely if no external funding applies.

---

## Manuscript DOI (MEDIUM — cannot be assigned until after journal acceptance)

In `DEPOSITION_METADATA.json` `related_identifiers[0]` and in `README_zenodo.md` Citation section:

```
"[manuscript DOI — to be assigned upon publication]"
```

This is expected to remain a placeholder until after journal acceptance. Zenodo allows updating metadata after deposition, so this can be filled post-acceptance.

---

## Publication Date (LOW — update at upload time)

In `DEPOSITION_METADATA.json`:

```json
"publication_date": "[UPDATE TO ACTUAL UPLOAD DATE — format: YYYY-MM-DD]"
```

Set this to the actual date when uploading to Zenodo (e.g., "2026-07-15"). The field accepts YYYY-MM-DD only.

---

## Manuscript Citation Block (LOW — author decision)

In `README_zenodo.md`:

```
> [Authors]. PAD4-DB: a curated structure–activity resource reveals...
```

Replace [Authors] and [Journal], [Year] when manuscript is accepted.

---

## Summary Table

| Item | Blocking deposition? | Author action required? |
|------|---------------------|------------------------|
| Author names (×3) | ✅ YES | ✅ YES |
| Affiliations (×3) | ✅ YES | ✅ YES |
| ORCIDs (×3) | ✅ YES (or remove) | ✅ YES |
| Contact name + email in README | ✅ YES | ✅ YES |
| ChEMBL release version | ⚠ HIGH | ✅ YES — check ebi.ac.uk |
| Funding | NO | ✅ YES (or remove field) |
| Manuscript DOI | NO (post-acceptance) | ✅ YES — fill after journal acceptance |
| Publication date | NO (set at upload) | ✅ YES — set at upload time |
| Citation block in README | NO | ✅ YES — fill after acceptance |

**Estimated time to complete all blocking items:** 15–30 minutes once author metadata are in hand.

---

## Items Fully Resolved (no further action needed)

| Item | Status |
|------|--------|
| License field | ✅ CC-BY-SA-4.0 |
| ChEMBL attribution text | ✅ Present (version placeholder remains) |
| BindingDB attribution text | ✅ Present |
| PubChem attribution text | ✅ Present |
| Software versions in README | ✅ All 7 packages listed |
| GitHub URL in metadata | ✅ Present |
| GitHub URL in README | ✅ Present |
| ShareAlike notice | ✅ Present |
| Keywords | ✅ 15 keywords |
| Zenodo JSON validity | ✅ Validated |
