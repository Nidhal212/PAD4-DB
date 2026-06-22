# Supplementary Materials — PAD4-DB v1.0

## Supplementary Tables (`tables/`)

| File | Label | Description | Rows |
|------|-------|-------------|------|
| `Table_S-cliffs.csv` | Table S-cliffs | All 867 activity cliff pairs (94 severe, 193 moderate, 580 broad); Tanimoto, ΔpIC50, cliff tier, MMP flag, ECFP4/ECFP6 Tanimoto, ecfp4_only_cliff | 867 |
| `Table_S-scaf.csv` | Table S-scaf | 375 multi-member Murcko scaffold series: size, cliff density, smoothness score, mean/std pIC50 | 375 |
| `Table_S-null.csv` | Table S-null | Three-null permutation comparison (unrestricted/scaffold/assay-constrained) for cliff rarity and hub concentration | 3 |
| `Table_S-indep.csv` | Table S-indep | Per-compound source list, count, independence score, concordance flags | 3,093 |
| `Table_S-hubs.csv` | Table S-hubs | The 4 cliff-hub compounds with class labels, pIC50, MMP degree | 4 |
| `Table_S-mmp.csv` | Table S-mmp | Severe-cliff ΔpIC50 by MMP transformation category | — |
| `Table_S-assay.csv` | Table S-assay | Within-assay cliff robustness (88/94 = 94% within genuine PubChem AID) | 3 |
| `Table_S-mw.csv` | Table S-mw | ΔMW sensitivity: 3 bands (<20 / 20-50 / ≥50 Da), mean |ΔpIC50| per band | 3 |
| `Table_S-indep_robust.csv` | Table S-indep_robust | Independence score robustness under ±20% weight perturbation (1,000 permutations) | — |
| `Table_S-threshold.csv` | Table S-threshold | Cliff count sensitivity to Tanimoto and ΔpIC50 thresholds | — |
| `Table_S-qc.csv` | Table S-qc | Dataset QC/provenance ledger: raw → final compound counts with reasons for every exclusion | — |

## Supplementary Datasets (`datasets/`)

| File | Description |
|------|-------------|
| `standardized_structures.sdf` | All 3,093 standardized structures in RDKit V2000 SDF; SD tags: pIC50, Murcko scaffold SMILES, source-independence score, hub_class, mechanism_class |
| `hub_neighborhood_metrics.csv` | Per-hub neighborhood metrics: potency gradient, neighbor count, similarity, scaffold occupancy, local cliff count |
| `hub_physchem_table.csv` | Hub vs other-cliff physicochemical comparison (Mann–Whitney U + Benjamini–Hochberg FDR at 3 degree thresholds) |

## Canonical Counts

- 3,093 compounds · 94 severe activity cliffs · 4 hub compounds · 1,244 unique Murcko scaffolds
- All counts verified against source files (see `manuscript/FREEZE_MANIFEST_v1.0.md`)
