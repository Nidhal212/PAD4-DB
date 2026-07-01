# PAD4-DB

**A curated structure–activity resource for PAD4 inhibitors: hub-organized activity cliffs and scaffold-dependent SAR ruggedness**

[![License: CC-BY-SA-4.0](https://img.shields.io/badge/License-CC--BY--SA--4.0-lightgrey.svg)](LICENSE)
[![Pipeline: Snakemake](https://img.shields.io/badge/workflow-Snakemake-green)](Snakefile)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue)](environment.yml)

> Manuscript under review. Dataset and code are released ahead of publication for transparency and reproducibility.

---

## Overview

PAD4-DB is a standardized, deduplicated database of **3,093 human PAD4 (PADI4) inhibitors** assembled from three public bioactivity sources: PubChem bioassay campaigns (95 AIDs), ChEMBL (CHEMBL6111), and BindingDB (UniProt Q9UM07). The dataset provides:

- Consensus pIC50 values (range 2.00–8.52, median 6.84) from dose-response IC50 measurements
- Source provenance and independence scores (distinguishing pipeline re-curation from independent replication)
- Systematic activity cliff characterization (94 severe pairs, Tanimoto ≥ 0.8, ΔpIC50 ≥ 2.0)
- Scaffold annotations (1,244 unique Murcko scaffolds; 375 series ≥ 2 members)
- MMP-confirmed cliff pairs (80/94 severe pairs, 85.1%)
- Cliff-hub annotations for two structural classes covering 53.2% of severe cliff pairs
- A 327,336-compound HTS structural reference

The full pipeline is implemented as a Snakemake workflow and reproduces all results from raw source files in under 10 minutes on 4 cores.

---

## Repository Structure

```
PAD4-DB/
├── Snakefile                    End-to-end Snakemake workflow (11 rules)
├── environment.yml              Conda environment specification
├── config/
│   └── config.yaml              AID lists, QC thresholds, canonical assertions
├── workflow/
│   └── envs/pad4bench.yaml      Pinned environment (Python 3.10, RDKit 2025.09.5)
│
├── scripts/                     Pipeline scripts (standalone + Snakemake)
│   ├── 01_standardize/          SMILES standardization, salt stripping, InChIKey
│   ├── 02_normalize/            Unit conversion → pIC50 (nM scale)
│   ├── 03_aggregate/            Replicate aggregation, potency/HTS split
│   ├── 04_dedup/                Cross-AID deduplication, compound assembly
│   ├── 05_cliffs/               Scaffold analysis, activity cliffs, MMP
│   ├── audit/                   Biological audit (A1) and reference recovery (A2)
│   └── nature/                  Figure generation scripts (Nature style)
│
├── data/
│   └── processed/               Released dataset files
│       ├── pad4_compounds.parquet          3,093 PAD4 inhibitors  ← main file
│       ├── hts_compound_index.parquet      327,336 HTS compounds
│       ├── activity_cliffs.parquet         867 cliff pairs
│       └── activity_pairs_sim_ge06.parquet 358,416 similarity pairs
│
└── outputs/
    └── tables/                  Pipeline QC and analysis tables
        ├── 01_standardization_report.csv
        ├── 02_normalization_report.csv
        ├── 05_scaffold_summary.csv
        ├── 05_cliff_summary.json
        └── supp_*.csv           Supplementary analysis tables
```

---

## Quick Start

### Load the dataset

```python
import pandas as pd

# Main compound dataset
df = pd.read_parquet("data/processed/pad4_compounds.parquet")
print(f"{len(df)} PAD4 inhibitors, pIC50 {df.pic50_consensus.min():.2f}–{df.pic50_consensus.max():.2f}")

# Activity cliffs
cliffs = pd.read_parquet("data/processed/activity_cliffs.parquet")
severe = cliffs[cliffs.cliff_tier == "severe"]
print(f"{len(severe)} severe cliff pairs (Tanimoto ≥ 0.8, ΔpIC50 ≥ 2.0)")
```

### Reproduce the full pipeline

```bash
# 1. Create environment
conda env create -f environment.yml
conda activate pad4bench
pip install "snakemake==7.32.4" "pulp==2.7.0"

# 2. Download raw data into data/raw/
#    - PubChem: https://pubchem.ncbi.nlm.nih.gov/bioassay/
#    - ChEMBL:  https://www.ebi.ac.uk/chembl/  (assay CHEMBL6111)
#    - BindingDB: https://www.bindingdb.org    (UniProt Q9UM07)

# 3. Dry-run to verify DAG
snakemake --dry-run --cores 4

# 4. Run full pipeline (~7 min on 4 cores)
snakemake --cores 4
```

---

## Dataset: `pad4_compounds.parquet`

Primary file — 3,093 rows, one per unique compound (by InChIKey).

| Column | Description |
|--------|-------------|
| `inchi_key` | Standard InChIKey — primary identifier |
| `smiles_std` | Canonical SMILES (RDKit 2025.09.5, salt-stripped) |
| `pic50_consensus` | Consensus pIC50 = −log₁₀(IC50 / M); median 6.84 |
| `source_list` | Pipe-delimited source databases (e.g. `bindingdb\|chembl\|pubchem_confirmatory`) |
| `source_independence_score` | 0.3–1.0; ≥ 0.6 indicates genuinely independent multi-source data |
| `multi_source` | True if present in ≥ 2 databases (89.1%) |
| `high_confidence` | True if concordant AND source_independence_score ≥ 0.3 (88.8%) |
| `concordant` | True if max cross-source ΔpIC50 ≤ 1.0 log unit (99.7%) |
| `patent_exclusive` | True if present only in PubChem patent-deposited data (n = 233) |

For the full data dictionary see [`publication/DATA_DICTIONARY.md`](publication/DATA_DICTIONARY.md).

---

## Key Results

| Metric | Value |
|--------|-------|
| Total compounds | 3,093 |
| HTS reference compounds | 327,336 |
| pIC50 range | 2.00 – 8.52 |
| pIC50 median | 6.84 |
| Severe activity cliffs | 94 pairs (Tanimoto ≥ 0.8, ΔpIC50 ≥ 2.0) |
| MMP-confirmed severe cliffs | 80 / 94 (85.1%) |
| Cliff-hub compounds | 4 (2 classes, 53.2% of severe pairs) |
| Unique Murcko scaffolds | 1,244 |
| Compounds in scaffold series | 71.9% |
| Pipeline validation | 47 / 47 canonical metrics reproduced |

---

## Source Data

Raw data is not included in this repository (total ~44 MB). Download instructions:

| Source | Target | URL |
|--------|--------|-----|
| PubChem | 95 bioassay CSVs | Listed in `config/config.yaml` |
| ChEMBL | CHEMBL6111 activity export | https://www.ebi.ac.uk/chembl/ |
| BindingDB | UniProt Q9UM07 export | https://www.bindingdb.org |

---

## Citation

> [Authors]. PAD4-DB: a curated structure–activity resource reveals hub-organized activity cliffs and scaffold-dependent SAR ruggedness in PAD4 inhibitors. *[Journal]*, [Year]. DOI: [to be assigned after acceptance].

BibTeX and full citation metadata: [`publication/CITATION.cff`](publication/CITATION.cff)

---

## License

Dataset (data/processed/) and outputs (outputs/tables/) are released under [CC-BY-SA-4.0](LICENSE).  
Code (scripts/, Snakefile) is released under [MIT License](LICENSE).
