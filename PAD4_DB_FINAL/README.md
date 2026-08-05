# PAD4-DB — Ultimate Final Version

**A curated structure–activity resource reveals hub-organized activity cliffs and scaffold-dependent SAR ruggedness in PAD4 inhibitors**

| Item | Value |
|------|-------|
| Compounds | 3,093 PAD4 inhibitors |
| HTS screened | 327,336 compounds |
| Sources | PubChem (95 AIDs) · ChEMBL (CHEMBL6111) · BindingDB (Q9UM07) |
| pIC50 range | 2.00 – 8.52 (median 6.84) |
| Severe activity cliffs | 94 pairs (Tan ≥ 0.8, ΔpIC50 ≥ 2.0) |
| Unique scaffolds | 1,244 (375 series, 869 singletons) |
| Pipeline | Snakemake 7.32.4, fully automated |
| Validated | 47/47 canonical metrics reproduced exactly (2026-07-01) |

---

## Directory Layout

```
PAD4_DB_FINAL/
├── README.md               this file
├── CLAUDE.md               canonical numbers, locked statements, project guide
├── LICENSE                 CC-BY-SA-4.0
│
├── raw_data/               SOURCE DATA AS DOWNLOADED
│   ├── pubchem/
│   │   ├── hts/            3 HTS assays (AIDs 463073, 485272, 488796)
│   │   ├── confirmatory/   57 Layer A confirmatory assays
│   │   ├── literature_derived/  11 Layer C assays
│   │   └── secondary/      26 Layer D+E assays
│   ├── chembl/             CHEMBL6111 (4,925 rows)
│   └── bindingdb/          Q9UM07 (3,087 rows)
│
├── pipeline/               REPRODUCIBLE PIPELINE
│   ├── Snakefile           11-rule end-to-end Snakemake workflow
│   ├── environment.yml
│   ├── config/config.yaml  AID lists, thresholds, canonical assertions
│   ├── workflow/envs/pad4bench.yaml  pinned env (Python 3.10, RDKit 2025.09.5)
│   └── scripts/
│       ├── 01_standardize/   SMILES standardization + InChIKey
│       ├── 02_normalize/     unit conversion to pIC50
│       ├── 03_aggregate/     replicate aggregation + space split
│       ├── 04_dedup/         cross-AID dedup + compound assembly
│       ├── 05_cliffs/        scaffold + activity cliff + MMP analysis
│       ├── audit/            A1 biological audit, A2 reference recovery
│       └── nature/           all figure scripts (Nature style)
│
├── data/                   PROCESSED OUTPUTS
│   ├── interim/standardized/  standardized_compounds.parquet (341,282 rows)
│   ├── interim/normalized/    normalized → aggregated → potency/hts split
│   └── processed/
│       ├── pad4_compounds.parquet        MAIN DATASET (3,093 compounds)
│       ├── hts_compound_index.parquet    (327,336 HTS compounds)
│       ├── activity_cliffs.parquet       (867 cliff pairs)
│       └── activity_pairs_sim_ge06.parquet  (358,416 similarity pairs)
│
├── analysis/               QC REPORTS AND AUDIT
│   ├── tables/             pipeline QC tables + supplementary analysis CSVs
│   └── audit/              A1 biological audit, A2 reference recovery
│
├── figures/                FINAL FIGURES (Nature style, 600 dpi)
│   ├── main/               6 main figures (PNG + PDF)
│   └── supplementary/      5 supplementary figures (PNG + PDF)
│
└── manuscript/             FINAL MANUSCRIPT
    ├── PAD4_DB_manuscript_FINAL.md      Markdown source
    ├── PAD4_DB_manuscript_FINAL.docx    Word (for journal submission)
    └── PAD4_DB_manuscript_FINAL.pdf     PDF
```

---

## Reproduce the Full Pipeline

```bash
# 1. Install environment
conda env create -f pipeline/environment.yml
conda activate pad4bench
pip install "snakemake==7.32.4" "pulp==2.7.0"

# 2. Run from project root (PAD4_DB_FINAL/)
snakemake -s pipeline/Snakefile --cores 4

# Dry-run first
snakemake -s pipeline/Snakefile --dry-run --cores 4
```

---

## Canonical Numbers (locked 2026-07-01, verified by Snakemake re-run)

| Metric | Value |
|--------|-------|
| Total raw rows ingested | 341,282 |
| Rows with valid SMILES | 341,276 (100.0%) |
| Unique InChIKeys (raw) | 328,976 |
| norm_status OK | 338,021 (99.0%) |
| Replicate groups | 339,687 |
| **Final compound count** | **3,093** |
| pIC50 range | 2.00 – 8.52 |
| pIC50 median | 6.84 |
| multi_source | 89.1% |
| Severe activity cliffs | 94 |
| max ΔpIC50 severe | 3.045 |
| Unique scaffolds | 1,244 |
| Largest scaffold series | 174 compounds |
