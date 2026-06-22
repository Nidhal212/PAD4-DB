# PAD4-DB

A standardized, deduplicated knowledge base of **3,093 structurally resolved PAD4 (PADI4) inhibitors** with consensus pIC50 values, assembled from PubChem bioassay campaigns, ChEMBL (CHEMBL6111), and BindingDB (Q9UM07). Includes scaffold annotations, source-independence scores, the severe activity-cliff network, four hub-compound labels, a per-scaffold ruggedness ranking, and matched-molecular-pair transformations.

**Manuscript:** under review. Will be linked here upon publication.

## Repository layout

```
PAD4_DB_FREEZE_v1.0/
├── figures/main/           6 main figures (PNG, 600 dpi)
├── figures/supplementary/  5 supplementary figures
├── supplementary/tables/   11 supplementary tables (CSV)
├── supplementary/datasets/ standardized_structures.sdf + hub metrics
├── reproducibility/code/   full pipeline + figure scripts + environment
└── deposition/zenodo/      self-contained Zenodo deposition package

publication/
├── data/                   pad4_compounds.parquet · activity_cliffs.parquet
│                           hts_compound_index.parquet · mmp_pairs · etc.
├── scripts/                canonical pipeline (stages 01–05) + figure scripts
├── environment.yml         pinned conda environment (pad4bench)
└── audit/                  verification reports and audit notes
```

## Canonical numbers (validated)

3,093 compounds · 94 severe activity cliffs · 4 hub compounds · 1,244 Murcko scaffolds · pIC50 median 6.84

## Reproduce

```bash
conda env create -f publication/environment.yml
conda activate pad4bench
# Run any pipeline stage, e.g.:
python publication/scripts/05_cliffs/05_scaffold_and_cliffs.py
```

## Citation

> [Authors]. PAD4-DB: a curated structure–activity resource reveals hub-organized activity cliffs and scaffold-dependent SAR ruggedness in PAD4 inhibitors. *[Journal]*, [Year]. DOI: [to be assigned].

## License

CC-BY-4.0.
