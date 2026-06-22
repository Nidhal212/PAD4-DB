# PAD4-DB: A Curated Structure–Activity Resource for PAD4 Inhibitors

**PAD4-DB** is a standardized, deduplicated knowledge base of **3,093 structurally resolved PAD4 (PADI4) inhibitors** with consensus pIC50 values, plus a parallel high-throughput-screening layer of 327,336 screened compounds. The resource includes Murcko scaffold annotations, per-compound source-independence scores, the severe activity-cliff network, four hub-compound labels, a per-scaffold ruggedness ranking, and matched-molecular-pair transformation typology.

## Package contents (release candidate — code + data only)

| Directory | Contents |
|-----------|----------|
| `data/` | All machine-readable data files (CSV + 3,093-structure SDF) — see `data/README_data.md` |
| `code/` | Reproducibility/analysis scripts + `environment.yml` — see `code/README_code.md` |
| `DEPOSITION_METADATA.json` | Zenodo metadata template |

**Note:** the manuscript (MD/DOCX/PDF) and rendered figures are **held back** until the
publication/submission strategy is finalized, and are therefore not part of this release
candidate. All figures are regenerable from the released data and scripts; the manuscript
and figures will be added (or linked by DOI) upon publication.

## Canonical numbers (validated)
3,093 compounds · 94 severe activity cliffs · 4 cliff-hub compounds · 1,244 unique Murcko scaffolds.

## Key findings
- Cross-source agreement (99.7% concordance) largely reflects shared re-curation, not independent replication; only 17.1% of compounds are source-independent.
- The severe-cliff landscape is hub-organized: four compounds in two classes account for 53.2% of severe cliffs, an effect robust to scaffold- and assay-constrained permutation nulls.
- SAR ruggedness is scaffold- and vector-specific: 96.8% of scaffold series are perfectly smooth.

## Citation
If you use PAD4-DB, please cite:

> [Authors]. PAD4-DB: a curated structure–activity resource reveals hub-organized activity cliffs and scaffold-dependent SAR ruggedness in PAD4 inhibitors. [Journal], [Year]. DOI: [to be assigned].

## License
CC-BY-4.0 (see `DEPOSITION_METADATA.json`).

## Contact
[Corresponding author — name and email to be inserted before deposition].

## Provenance
Source bioactivity data derive from PubChem BioAssay, ChEMBL (assay CHEMBL6111), and BindingDB (UniProt Q9UM07), which remain available from the respective repositories.
