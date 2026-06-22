# PAD4-DB v1.0: A Curated Structure–Activity Resource for PAD4 Inhibitors

**Version:** 1.0  
**License:** CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/)  
**Corresponding author:** Nidhal Tarhouni — Nidhal1t@gmail.com  
                          https://orcid.org/0000-0001-5042-2947  
**GitHub repository:** https://github.com/Nidhal212/PAD4-DB  
**Manuscript DOI:** [to be assigned upon publication]

---

## What is PAD4-DB?

PAD4-DB is a standardized, deduplicated knowledge base of **3,093 structurally resolved PAD4
(PADI4) inhibitors** with consensus pIC50 values, assembled from three public bioactivity
databases (PubChem, ChEMBL, BindingDB) through a fully scripted, fail-loud curation pipeline.
A parallel high-throughput screening (HTS) layer covers 327,336 screened compounds.

The resource includes:
- RDKit-standardized structures (3,093 compounds, SDF format)
- Consensus pIC50 values (range 2.00–8.52, median 6.84)
- Murcko scaffold annotations and per-scaffold SAR ruggedness metrics
- Source-independence scores (17.1% of compounds independently replicated across databases)
- Activity cliff network (94 severe pairs, 4 hub compounds, 53.2% hub share)
- Matched-molecular-pair transformation typology (80/94 severe cliffs MMP-confirmed)
- Three-null permutation framework output

---

## What Is Inside

### Compound counts
| Layer | Compounds |
|-------|-----------|
| Dose-response (pIC50) database | **3,093** |
| HTS screening layer | **327,336** |

### pIC50 definition
pIC50 = −log₁₀(IC50 in molar units). IC50 values were converted to nM, then
pIC50 = −log₁₀(IC50_nM × 10⁻⁹). Where multiple measurements existed for the same
compound × assay × source combination, replicates were aggregated by geometric mean
(log-space arithmetic mean). Consensus pIC50 across sources was computed as the
unweighted mean of per-source pIC50 values. Range: 2.00–8.52, median 6.84.

### Structural fingerprints
Activity cliffs and Tanimoto similarities were computed using:
- **Fingerprint type:** ECFP4 (Extended Connectivity Fingerprint, circular)
- **Algorithm:** Morgan algorithm (RDKit implementation)
- **Radius:** 2 (equivalent to ECFP4)
- **Bit vector length:** 2048 bits
- **Software:** RDKit 2025.09.5

Activity cliff definition: Tanimoto similarity ≥ 0.8 AND |ΔpIC50| ≥ 2.0 log units (severe);
≥ 1.5 (moderate); ≥ 1.0 (broad).

---

## Package Contents

```
PAD4-DB v1.0/
├── README.md                           This file
├── environment.yml                     Pinned conda environment (pad4bench)
├── requirements.txt                    Pip-compatible package list
│
├── data/
│   ├── PAD4_DB_pIC50.csv              3,093 PAD4 inhibitors with consensus pIC50
│   ├── PAD4_HTS_327k.csv              327,336 HTS-screened compounds
│   └── standardized_structures.sdf    3,093 structures in SDF format (RDKit V2000)
│
├── analysis/
│   ├── activity_cliffs.csv            867 cliff pairs (94 severe / 193 moderate / 580 broad)
│   ├── cliff_hubs.csv                 4 hub compounds with class labels and SMILES
│   ├── scaffold_ruggedness.csv        375 scaffold series with ruggedness metrics
│   ├── mmp_transformations.csv        MMP transformation impact by category
│   └── source_independence_scores.csv Per-compound source provenance and independence scores
│
└── code/
    ├── pipeline_step1_standardization.py   SMILES standardization (RDKit)
    ├── pipeline_step2_normalization.py     Activity normalization (unit conversion → pIC50)
    ├── cliff_analysis.py                   Activity cliff detection + scaffold analysis
    ├── scaffold_analysis.py                Per-scaffold ruggedness metrics
    └── null_models.py                      Three-null permutation framework
```

### PAD4_DB_pIC50.csv column descriptions
| Column | Description |
|--------|-------------|
| inchi_key | Standard InChIKey (computed by RDKit, 27-char) |
| smiles_std | RDKit-standardized canonical SMILES |
| pIC50 | Consensus pIC50 (unweighted mean across sources) |
| pic50_min | Minimum pIC50 across all measurements |
| pic50_max | Maximum pIC50 across all measurements |
| n_sources | Number of databases contributing measurements |
| source_list | Pipe-separated source labels |
| source_independence_score | Independence score (0.3–1.0; see paper) |
| mol_weight | Molecular weight (Da) |
| n_heavy_atoms | Heavy atom count |
| mechanism_class | Assay mechanism classification |
| is_covalent | SMARTS-based covalent warhead flag |
| warhead_class | Warhead type if covalent |
| hub_class | Activity cliff hub class (A / B / none) |
| fragment_flag | Fragment-like compound flag (MW<200 AND pIC50<4.0) |
| high_confidence | High-confidence flag (multi-source, concordant) |
| multi_source | Present in ≥2 source databases |

### PAD4_HTS_327k.csv column descriptions
| Column | Description |
|--------|-------------|
| inchi_key | Standard InChIKey |
| smiles_std | RDKit-standardized SMILES |
| n_hts_assays | Number of HTS assays in which compound was screened |
| max_pct_inh | Maximum % inhibition across all HTS assays |
| any_active | True if max_pct_inh ≥ 50% in any assay |
| hts_activity_score | Aggregate HTS activity score |
| hts_consensus_confidence | Confidence in HTS consensus outcome |
| hts_outcome | Consensus HTS outcome (active/inactive/inconclusive) |
| confirmed_in_potency_space | True if compound also has pIC50 in PAD4_DB_pIC50.csv |

---

## How to Reproduce

### Option A — Conda (recommended)

```bash
# 1. Create the environment
conda env create -f environment.yml
conda activate pad4bench

# 2. Run pipeline (from project root, with raw data in data/raw/)
python code/pipeline_step1_standardization.py
python code/pipeline_step2_normalization.py

# 3. Run analysis
python code/cliff_analysis.py
python code/scaffold_analysis.py
python code/null_models.py
```

### Option B — Pip

```bash
pip install -r requirements.txt
```

Note: RDKit via pip requires Python ≥ 3.8. If installation fails, use Option A (conda).

### Notes on running from this Zenodo package
The pipeline scripts reference input paths relative to the full project directory.
To run from a relocated copy, update the `ROOT` variable at the top of each script
to point to the directory containing `data/`. Permutation tests use a fixed random
seed (42) for reproducibility.

---

## Software Versions (locked)

| Software | Version |
|----------|---------|
| Python | 3.10.19 |
| RDKit | 2025.09.5 |
| pandas | 2.3.3 |
| NumPy | 2.2.5 |
| SciPy | 1.15.3 |
| matplotlib | 3.10.8 |
| networkx | 3.4.2 |

---

## Data Provenance and Source Attribution

All structures have been standardized using RDKit 2025.09.5 (SMILES canonicalization,
desalting, charge normalization, new 2D coordinate generation via RDKit's 2D layout engine).
Original source identifiers (SIDs, CIDs, ChEMBL IDs, BindingDB IDs) are not retained
in the released files. Provenance is recorded as a source-label string per compound
(e.g., "bindingdb|chembl|pubchem_confirmatory") in analysis/source_independence_scores.csv.

### PubChem BioAssay
- **Source:** PubChem BioAssay, NCBI/NLM/NIH — https://pubchem.ncbi.nlm.nih.gov
- **Data retrieved:** 2026-06-10 to 2026-06-14
- **Coverage:** 95 assay IDs (57 confirmatory, 11 literature-derived, 26 secondary/HTS, 3 HTS)
- **License:** Public domain (US government work, 17 U.S.C. § 105)
- **Attribution:** Wang Y, et al. *Nucleic Acids Res.* 2019;47(D1):D1102–D1109.
  doi:10.1093/nar/gky1033

### ChEMBL
- **Source:** ChEMBL 34, EMBL-EBI — https://www.ebi.ac.uk/chembl/
- **Assay:** CHEMBL6111 (Protein-arginine deiminase type-4, Homo sapiens)
- **Data retrieved:** 2026-06-14
- **License:** CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/)
- **⚠ ShareAlike notice:** Derivative databases including these structures must use CC BY-SA 4.0.
- **Attribution:** Mendez D, et al. *Nucleic Acids Res.* 2019;47(D1):D930–D940.
  doi:10.1093/nar/gky1075

### BindingDB
- **Source:** BindingDB — https://www.bindingdb.org (UniProt Q9UM07)
- **Data retrieved:** 2026-06-10
- **License:** CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
- **Attribution:** Gilson MK, et al. *Nucleic Acids Res.* 2016;44(D1):D1045–D1053.
  doi:10.1093/nar/gkv1072

---

## Canonical Numbers (validated)

3,093 compounds · 94 severe activity cliffs · 4 cliff-hub compounds · 1,244 unique Murcko scaffolds · pIC50 median 6.84

---

## Citation

> Tarhouni N. PAD4-DB: a curated structure–activity resource reveals hub-organized activity
> cliffs and scaffold-dependent SAR ruggedness in PAD4 inhibitors. *[Journal]*, [Year].
> DOI: [to be assigned].

Also cite PubChem, ChEMBL, and BindingDB as listed above.

---

## License

**Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)**  
https://creativecommons.org/licenses/by-sa/4.0/

ShareAlike requirement derives from ChEMBL (CC BY-SA 4.0). PubChem (public domain) and
BindingDB (CC BY 4.0) are both compatible with CC BY-SA 4.0.

---

## Contact

Nidhal Tarhouni, PhD (Postdoctoral Researcher)  
Laboratory of Enzyme Engineering and Microbiology  
National Engineering School of Sfax (ENIS), University of Sfax  
P.O. Box 1173, Sfax 3038, Tunisia  
Email: Nidhal1t@gmail.com  
ORCID: https://orcid.org/0000-0001-5042-2947
