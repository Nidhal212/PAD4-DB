# PAD4-DB v1.0: A Curated Structure–Activity Resource for PAD4 Inhibitors

**Version:** 1.0  
**License:** CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/)  
**Corresponding author:** [Name] — [email]  
**GitHub repository:** https://github.com/Nidhal212/PAD4-DB  
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

**Note:** The manuscript and rendered figures are held back until publication. All figures are
regenerable from the released data and scripts; the manuscript will be linked by DOI upon
publication.

---

## Canonical Numbers (validated)

3,093 compounds · 94 severe activity cliffs · 4 cliff-hub compounds · 1,244 unique Murcko scaffolds

---

## Data Provenance and Source Attribution

All structures have been standardized using RDKit 2025.09.5 (SMILES canonicalization,
desalting, charge normalization, new 2D coordinate generation via RDKit's 2D layout engine).
Original source identifiers (SIDs, CIDs, ChEMBL IDs, BindingDB IDs) are not retained
in the released files. Provenance is recorded as a source-label string per compound
(e.g., "bindingdb|chembl|pubchem_confirmatory") in source_independence_scores.csv.

### PubChem BioAssay

- **Source:** PubChem BioAssay, National Center for Biotechnology Information (NCBI/NLM/NIH)
  https://pubchem.ncbi.nlm.nih.gov
- **Data retrieved:** 2026-06-10 to 2026-06-14
- **Coverage:** 95 assay IDs (57 confirmatory, 11 literature-derived, 26 secondary/HTS, 3 HTS)
- **License:** Public domain (US government work, 17 U.S.C. § 105 — no copyright restrictions)
- **Attribution:** Wang Y, et al. PubChem 2019 update: improved access/downloads, new
  compound comparison, periodic table page, and updated compound, substance and assay pages.
  *Nucleic Acids Res.* 2019;47(D1):D1102–D1109. doi:10.1093/nar/gky1033

### ChEMBL

- **Source:** ChEMBL database, EMBL-EBI, https://www.ebi.ac.uk/chembl/
- **Assay:** CHEMBL6111 (Protein-arginine deiminase type-4, Homo sapiens)
- **Release version:** [INSERT ChEMBL release number — check ebi.ac.uk/chembl/db_info/]
- **Data retrieved:** 2026-06-14
- **License:** CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/)
- **⚠ ShareAlike notice:** This dataset includes structures derived from ChEMBL. In
  accordance with ChEMBL's CC BY-SA 4.0 license, any derivative database that includes
  these structures must be distributed under CC BY-SA 4.0 or a compatible license.
- **Attribution:** Mendez D, et al. ChEMBL: towards direct deposition of bioassay data.
  *Nucleic Acids Res.* 2019;47(D1):D930–D940. doi:10.1093/nar/gky1075

### BindingDB

- **Source:** BindingDB, https://www.bindingdb.org
- **Query:** UniProt accession Q9UM07 (Protein-arginine deiminase type-4, Homo sapiens)
- **Data retrieved:** 2026-06-10
- **License:** CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
- **Modifications:** structures were standardized using RDKit 2025.09.5 (see above)
- **Attribution:** Gilson MK, et al. BindingDB in 2015: A public database for medicinal
  chemistry, computational chemistry and systems pharmacology.
  *Nucleic Acids Res.* 2016;44(D1):D1045–D1053. doi:10.1093/nar/gkv1072

---

## Software Provenance

All data processing, standardization, and analysis were performed with the following
software versions (pinned in `code/environment.yml`):

| Software | Version |
|----------|---------|
| Python | 3.10.19 |
| RDKit | 2025.09.5 |
| pandas | 2.3.3 |
| NumPy | 2.2.5 |
| SciPy | 1.15.3 |
| matplotlib | 3.10.8 |
| networkx | 3.4.2 |

Reproduce the environment:
```bash
conda env create -f code/environment.yml
conda activate pad4bench
```

---

## Citation

If you use PAD4-DB, please cite:

> [Authors]. PAD4-DB: a curated structure–activity resource reveals hub-organized activity
> cliffs and scaffold-dependent SAR ruggedness in PAD4 inhibitors. *[Journal]*, [Year].
> DOI: [to be assigned].

Also cite the primary source databases (PubChem, ChEMBL, BindingDB) as listed above.

---

## License

This dataset is licensed under **Creative Commons Attribution-ShareAlike 4.0 International
(CC BY-SA 4.0)** — https://creativecommons.org/licenses/by-sa/4.0/

**You are free to:**
- Share — copy and redistribute in any medium or format
- Adapt — remix, transform, and build upon for any purpose, including commercially

**Under the following terms:**
- **Attribution** — give appropriate credit, provide a link to the license, indicate
  if changes were made
- **ShareAlike** — if you remix, transform, or build upon this material, distribute
  your contributions under CC BY-SA 4.0

License compatibility note: PubChem (public domain) and BindingDB (CC-BY 4.0) are
compatible with CC BY-SA 4.0. The ShareAlike requirement originates from ChEMBL.

---

## Contact

[Corresponding author full name]  
[Title, Department, Institution]  
[email]
