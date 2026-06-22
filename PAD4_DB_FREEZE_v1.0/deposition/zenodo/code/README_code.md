# PAD4-DB — Code

Analysis and manuscript-building scripts for PAD4-DB.

## Environment
Reproduced via `environment.yml` (conda):

```bash
conda env create -f environment.yml
conda activate pad4bench
```

Core versions: Python 3.10.19, RDKit 2025.09.5, pandas 2.3.3, numpy 2.2.5, scipy 1.15.3, matplotlib 3.10.8, networkx 3.4.2.

## Scripts

| Script | Description |
|--------|-------------|
| `analysis/audit_scaffold_ruggedness.py` | Per-scaffold ruggedness metrics + scaffold diversity (Gini, Shannon). |
| `analysis/audit_constrained_permutation.py` | Unrestricted / scaffold- / assay-constrained permutation nulls for cliff rarity and hub concentration. |
| `analysis/audit_hub_characterization.py` | Hub physicochemical comparison (FDR) + Class A/B neighborhood validation. |
| `analysis/audit_mmp_typology.py` | Expanded matched-molecular-pair transformation typology of severe cliffs. |
| `analysis/audit_source_and_qc.py` | Source-independence validation, results revalidation, and QC/provenance ledger. |
| `analysis/audit_export_package.py` | Assembles the machine-readable supplementary data package. |

## Usage
Run from the project root with the environment active:

```bash
python analysis/audit_constrained_permutation.py
```

Scripts use absolute project paths in their original form; adjust the `ROOT` variable at the top of each script if running from a relocated copy. Permutation tests use a fixed random seed (42) for reproducibility.
