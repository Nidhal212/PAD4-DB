# PAD4-DB v2

**A provenance-first, curated database of 3,093 human PAD4 inhibitors
with activity cliff characterization and source independence scoring.**

[![DOI](https://zenodo.org/badge/DOI/[TBD].svg)](https://doi.org/[TBD])
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Overview

PAD4-DB v2 integrates bioactivity data from 95 PubChem bioassays,
ChEMBL (CHEMBL6111), and BindingDB (Q9UM07) into a single curated
SAR database for human peptidylarginine deiminase 4 (PAD4, UniProt Q9UM07).

**Key statistics:**
- 3,093 curated SAR compounds with consensus pIC50 values (range: 2.00–8.52)
- 327,336 HTS structural reference compounds from 3 PubChem campaigns
- 1,244 Bemis-Murcko scaffolds; 71.9% compound scaffold coverage; Gini = 0.532
- 94 severe activity cliff pairs (Tanimoto ≥ 0.8, ΔpIC50 ≥ 2.0 log units)
- 85/94 (90.4%) severe cliff pairs validated by matched molecular pair (MMP) analysis
- 528 genuinely independent multi-source compounds (source independence score ≥ 0.6)
- Four cliff-hub compounds in two structural classes account for 53.2% of severe pairs
- Full pipeline reproducibility in < 4 minutes on a standard workstation

## Repository Structure

```
PAD4-db_V2/
├── data/                         Primary data files (also on Zenodo)
│   ├── pad4_compounds.parquet        3,093 SAR compounds (main file)
│   ├── activity_cliffs.parquet         867 cliff pairs (severe/moderate/broad)
│   ├── activity_pairs_with_sali.parquet  358,416 Tanimoto≥0.6 pairs + SALI
│   ├── hts_compound_index.parquet    327,336 HTS compounds
│   ├── mmp_pairs_cliff99.csv              707 MMP pairs among cliff compounds
│   ├── mmp_discontinuity_scores.csv        99 per-compound MMP scores
│   └── fingerprint_sensitivity_94pairs.csv  ECFP4 vs ECFP6 for all 94 severe pairs
├── figures/
│   ├── main/                     Figures 1–6 (PNG 600 dpi + PDF)
│   └── supplementary/            Figures S1–S5 (PNG 600 dpi + PDF)
├── tables/                       Main text and supplementary tables (HTML + CSV)
├── scripts/
│   ├── 01_standardize/           SMILES standardization (Step 01)
│   ├── 02_normalize/             Activity normalization (Step 02)
│   ├── 03_aggregate/             Replicate aggregation + QC (Step 03)
│   ├── 04_dedup/                 Deduplication, assembly, scoring (Step 04)
│   ├── 05_cliffs/                Scaffold + cliff analysis (Step 05)
│   ├── audit/                    Biological, chemical, reference audits
│   └── figures/                  Figure and table generation scripts
├── audit/                        Audit trail and verification reports
├── manuscript/                   Manuscript draft
├── README.md                     This file
├── DATA_DICTIONARY.md            Column descriptions for all data files
├── CHANGELOG.md                  Version history
├── CITATION.cff                  Citation metadata (cff format)
└── environment.yml               Conda environment specification
```

## Quick Start

### 1. Set up the environment

```bash
conda env create -f environment.yml
conda activate pad4bench
```

### 2. Load the main database

```python
import pandas as pd

# Main compound file
df = pd.read_parquet('data/pad4_compounds.parquet')
print(f"{len(df):,} compounds, columns: {list(df.columns)}")

# Filter to high-confidence compounds
hc = df[df['high_confidence'] == True]
print(f"{len(hc):,} high-confidence compounds")

# Filter to genuinely multi-source
ms = df[df['source_independence_score'] >= 0.6]
print(f"{len(ms):,} truly independent multi-source compounds")
```

### 3. Load activity cliffs

```python
cliffs = pd.read_parquet('data/activity_cliffs.parquet')
severe = cliffs[cliffs['cliff_tier'] == 'severe']
print(f"{len(severe)} severe cliff pairs, max ΔpIC50 = {severe['delta_pic50'].max():.3f}")

# Hub compounds
hubs = df[df['hub_class'].isin(['A', 'B'])]
print(f"Hub compounds: {hubs[['hub_class','pic50_consensus']].to_string()}")
```

### 4. Reproduce all figures

```bash
# Pre-compute shared objects (fingerprints, embeddings)
python scripts/figures/precompute_shared.py

# Generate all main text figures
python scripts/figures/fig01_headline.py
python scripts/figures/fig02_source_overlap.py
python scripts/figures/fig03_potency.py
python scripts/figures/fig04_scaffold.py
python scripts/figures/fig05_cliff_network.py
python scripts/figures/fig06_mmp.py

# Supplementary figures
python scripts/figures/supp_s01_pipeline.py
python scripts/figures/supp_s02_sali.py
python scripts/figures/supp_s03_patent.py
python scripts/figures/supp_s04_independence.py
python scripts/figures/supp_s05_scaffold_structures.py
```

### 5. Reproduce the full pipeline from raw data

Raw source files (not included — too large for GitHub) must be downloaded first.
See `audit/DOWNLOAD_DATES.txt` for source URLs and download dates.

```bash
python scripts/01_standardize/01_standardize_smiles.py
python scripts/02_normalize/02_normalize_activities.py
python scripts/03_aggregate/03_replicate_aggregate.py
python scripts/03_aggregate/03a_split_spaces.py
python scripts/03_aggregate/03b_logspace_qc.py
python scripts/03_aggregate/03c_smiles_integrity.py
python scripts/04_dedup/04_dedup_and_assemble.py
python scripts/04_dedup/04b_add_independence_scores.py
python scripts/05_cliffs/05_scaffold_and_cliffs.py
```

## Data Files

| File | Rows | Description |
|------|------|-------------|
| `pad4_compounds.parquet` | 3,093 | Main SAR compound file with consensus pIC50, source annotations, scaffold assignments, cliff-hub labels, source independence scores |
| `activity_cliffs.parquet` | 867 | All cliff pairs (severe/moderate/broad tiers) with ECFP4 Tanimoto, ECFP6 Tanimoto, ΔpIC50, MMP status, ecfp4_only_cliff flag |
| `activity_pairs_with_sali.parquet` | 358,416 | All compound pairs with Tanimoto ≥ 0.6, source labels, SALI values |
| `hts_compound_index.parquet` | 327,336 | HTS structural reference with activity scores and confirmed_in_potency_space flag |
| `mmp_pairs_cliff99.csv` | 707 | MMP pairs among the 99 severe cliff compounds with change-type classification |
| `mmp_discontinuity_scores.csv` | 99 | Per-compound MMP discontinuity scores for cliff hub ranking |
| `fingerprint_sensitivity_94pairs.csv` | 94 | ECFP4 and ECFP6 Tanimoto for all 94 severe pairs; full Supplementary Table S6 |

See `DATA_DICTIONARY.md` for complete column descriptions.

## Source Independence Scoring

89.1% of PAD4-DB v2 compounds appear in two or more source databases, but this
reflects the architecture of ChEMBL and BindingDB as aggregators of PubChem
bioassay data, not independent experimental replication. The `source_independence_score`
column explicitly quantifies this:

| Score | Combination | Compounds | Interpretation |
|-------|-------------|-----------|----------------|
| 0.3 | BindingDB + ChEMBL + PubChem | 1,366 | Pipeline re-curation (all three aggregate PubChem) |
| 0.5 | BindingDB + PubChem | 1,199 | Pipeline re-curation |
| 0.6 | BindingDB + ChEMBL | 167 | Threshold-independent |
| 0.7 | ChEMBL + PubChem | 23 | Genuinely independent |
| 1.0 | Single source | 338 | No redundancy |

Use `source_independence_score >= 0.6` (n=528, 17.1%) to select compounds with
genuinely independent measurements.

## Activity Cliff Hubs

Four compounds collectively participate in 53.2% of severe cliff pairs:

| InChIKey | Class | pIC50 | Hub pairs | Archetype |
|----------|-------|-------|-----------|-----------|
| SMADULGDNOCLOP-GISFHXKWSA-N | A | 5.390 | 15 | Series-embedded potency floor |
| RAVBZQAQTVGKIV-XBPDSQQVSA-N | A | 5.341 | 12 | Series-embedded potency floor |
| UDCDEKJNAMHBFH-HSZRJFAPSA-N | B | 4.301 | 12 | Scaffold-singleton structural attractor |
| DVCKJOQIVOGXEI-XMMPIXPASA-N | B | 4.301 | 11 | Scaffold-singleton structural attractor |

Class A hubs are mid-potency members of the dominant 174-compound azaindole-benzimidazole
scaffold series. Class B hubs are scaffold singletons (cyclobutyl/cyclopentyl sulfonamide)
whose broad ECFP4 similarity to diverse active compounds generates cliff pairs across
multiple chemotypes. These represent two independent and mechanistically interpretable
failure modes for similarity-based ML models.

## Citation

If you use PAD4-DB v2, please cite:

> [Author names TBD]. PAD4-DB v2: A Provenance-First Database of PAD4 Inhibitors
> with Activity Cliff Characterization and Source Independence Scoring.
> *[Journal TBD]*, [Year TBD]. DOI: [TBD]

See `CITATION.cff` for machine-readable citation metadata.

## License

Data: Creative Commons Attribution 4.0 International (CC BY 4.0)
Code: MIT License

## Contact

[Corresponding author TBD]
