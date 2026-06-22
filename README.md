# PAD4-DB

A standardized, deduplicated knowledge base of **3,093 structurally resolved PAD4 (PADI4) inhibitors** with consensus pIC50 values, plus a parallel HTS layer of 327,336 screened compounds. Includes scaffold annotations, source-independence scores, the severe activity-cliff network, four hub-compound labels, a per-scaffold ruggedness ranking, and matched-molecular-pair transformations.

## Contents
- `publication/manuscript/` — manuscript (Markdown, DOCX, PDF)
- `publication/figures/` — main figures 1–6 and supplementary figures S1–S5
- `publication/scripts/` — analysis and manuscript-building scripts
- `outputs/supplementary_package/` — machine-readable data package (CSV + SDF)
- `zenodo/` — self-contained Zenodo deposition package
- `environment.yml` — conda environment

## Canonical numbers (validated)
3,093 compounds · 94 severe activity cliffs · 4 hub compounds · 1,244 Murcko scaffolds.

## Reproduce
```bash
conda env create -f environment.yml
conda activate pad4bench
python publication/scripts/analysis/audit_constrained_permutation.py   # example
```

## Citation
> [Authors]. PAD4-DB: a curated structure–activity resource reveals hub-organized activity cliffs and scaffold-dependent SAR ruggedness in PAD4 inhibitors. [Journal], [Year]. DOI: [to be assigned].

## License
CC-BY-4.0.
