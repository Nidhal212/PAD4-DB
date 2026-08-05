# README_zenodo.md — Release Revision
**Date:** 2026-06-22  
**Auditor role:** Identify gaps; provide corrected text. Do not change scientific content.

---

## Gaps Identified in Current README_zenodo.md

| Gap | Severity | Details |
|-----|----------|---------|
| License field incorrect | **BLOCKER** | States "CC-BY-4.0"; should be "CC-BY-SA 4.0" (see recommended_license_decision.md) |
| ChEMBL license not named | **HIGH** | CC-BY-SA 4.0 must be cited per ChEMBL attribution requirements |
| ChEMBL release version missing | **HIGH** | Required for reproducibility and attribution |
| ChEMBL URL missing | **MEDIUM** | Should link to ebi.ac.uk/chembl |
| BindingDB license not named | **MEDIUM** | CC-BY 4.0 should be cited |
| BindingDB URL missing | **MEDIUM** | Should link to bindingdb.org |
| PubChem URL missing | **LOW** | Recommended for attribution completeness |
| Download dates missing | **MEDIUM** | In manuscript Methods but not in dataset README |
| RDKit standardization not mentioned | **MEDIUM** | Attribution requires noting modifications to source structures |
| Contact info is placeholder | **HIGH** | "[Corresponding author — name and email to be inserted]" must be filled |
| Repository version not stated | **LOW** | "v1.0" should appear explicitly |
| GitHub URL missing | **MEDIUM** | Related repository should be linked |
| Redistribution caveat missing | **MEDIUM** | Downstream users should be notified of ChEMBL ShareAlike obligation |

---

## Corrected README_zenodo.md (full replacement text)

Replace the entire current `README_zenodo.md` with the following:

---

```markdown
# PAD4-DB v1.0: A Curated Structure–Activity Resource for PAD4 Inhibitors

**Version:** 1.0  
**License:** CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/)  
**Corresponding author:** [Name] ([email])  
**Repository:** https://github.com/Nidhal212/PAD4-DB  
**Manuscript DOI:** [to be assigned upon publication]

---

## What is PAD4-DB?

PAD4-DB is a standardized, deduplicated knowledge base of **3,093 structurally resolved PAD4 
(PADI4) inhibitors** with consensus pIC50 values, assembled from three public bioactivity 
databases (PubChem, ChEMBL, BindingDB) through a fully scripted, fail-loud curation pipeline.

The resource includes:
- RDKit-standardized structures (3,093 compounds, SDF format)
- Consensus pIC50 values (range 2.00–8.52, median 6.84)
- Murcko scaffold annotations and per-scaffold SAR ruggedness metrics
- Source-independence scores (17.1% of compounds independently replicated across databases)
- Activity cliff network (94 severe pairs, 4 hub compounds, 53.2% hub share)
- Matched-molecular-pair transformation typology (80/94 severe cliffs MMP-confirmed)
- Three-null permutation framework output

---

## Package Contents

| Directory | Contents |
|-----------|----------|
| `data/` | Machine-readable data files (CSV + SDF) — see `data/README_data.md` |
| `code/` | Reproducibility/analysis scripts + `environment.yml` — see `code/README_code.md` |

---

## Canonical Numbers (validated)

3,093 compounds · 94 severe activity cliffs · 4 cliff-hub compounds · 1,244 unique Murcko scaffolds

---

## Data Provenance and Source Licenses

Structures in this dataset derive from three public bioactivity databases. The standardized 
structures (standardized_structures.sdf) have been processed using RDKit 2025.09.5 
(SMILES canonicalization, desalting, charge normalization, new 2D coordinate generation). 
Original source identifiers have not been retained in the released files.

### PubChem BioAssay
- Source: PubChem BioAssay (https://pubchem.ncbi.nlm.nih.gov), NCBI/NLM/NIH
- Downloaded: 2026-06-10 to 2026-06-14
- 95 assay IDs covering confirmatory, secondary, HTS, and literature-derived assay layers
- License: **Public domain** (US government work, 17 U.S.C. § 105)
- No restrictions on redistribution or adaptation

### ChEMBL
- Source: ChEMBL release [INSERT VERSION] (EMBL-EBI, https://www.ebi.ac.uk/chembl/)
- Assay: CHEMBL6111
- Downloaded: 2026-06-14
- License: **CC BY-SA 4.0** (https://creativecommons.org/licenses/by-sa/4.0/)
- ⚠ ShareAlike notice: This dataset includes structures derived from ChEMBL. In accordance 
  with ChEMBL's CC BY-SA 4.0 license, derivative databases that include these structures 
  must also be licensed under CC BY-SA 4.0 or a compatible license.
- Attribution: Mendez D, et al. ChEMBL: towards direct deposition of bioassay data. 
  Nucleic Acids Res. 2019;47(D1):D930-D940. doi:10.1093/nar/gky1075

### BindingDB
- Source: BindingDB (https://www.bindingdb.org), UniProt accession Q9UM07
- Downloaded: 2026-06-10
- License: **CC BY 4.0** (https://creativecommons.org/licenses/by/4.0/)
- Modifications: structures were standardized using RDKit 2025.09.5
- Attribution: Gilson MK, et al. BindingDB in 2015: A public database for medicinal 
  chemistry, computational chemistry and systems pharmacology. 
  Nucleic Acids Res. 2016;44(D1):D1045-D1053. doi:10.1093/nar/gkv1072

---

## License

This dataset is licensed under **Creative Commons Attribution-ShareAlike 4.0 International 
(CC BY-SA 4.0)**.

You are free to:
- Share — copy and redistribute the material in any medium or format
- Adapt — remix, transform, and build upon the material for any purpose, even commercially

Under the following terms:
- **Attribution** — give appropriate credit, link to the license, indicate changes made
- **ShareAlike** — if you remix, transform, or build upon this material, distribute your 
  contributions under the same CC BY-SA 4.0 license

Full license text: https://creativecommons.org/licenses/by-sa/4.0/

---

## Citation

If you use PAD4-DB, please cite:

> [Authors]. PAD4-DB: a curated structure–activity resource reveals hub-organized activity 
> cliffs and scaffold-dependent SAR ruggedness in PAD4 inhibitors. [Journal], [Year]. 
> DOI: [to be assigned].

Also cite the source databases as listed under Data Provenance above.

---

## Software

Structures were standardized using RDKit 2025.09.5 (https://www.rdkit.org).  
Pipeline environment: Python 3.10.19, conda (see `code/environment.yml`).

---

## Contact

[Corresponding author full name]  
[Title, Department, Institution]  
[email]
```

---

## Change Summary

| Section | Change |
|---------|--------|
| License | CC-BY-4.0 → CC-BY-SA 4.0 |
| Header | Added version, license, contact, GitHub URL |
| Provenance block | Added per-source license names, download dates, URLs, modification notice |
| ChEMBL | Added ShareAlike notice, ChEMBL DOI attribution |
| BindingDB | Added CC-BY 4.0 statement, Gilson 2016 citation |
| PubChem | Added public domain statement, URL, date |
| Software | Added RDKit version and URL |
| Contact | Placeholder to fill |
| Manuscript DOI | Placeholder to fill |
| GitHub URL | Added |
