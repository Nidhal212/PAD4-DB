# PAD4-DB — Data files

Complete machine-readable data package for PAD4-DB (3,093 PAD4 inhibitors).

## Files

| File | Rows | Description |
|------|------|-------------|
| `standardized_structures.sdf` | 3,093 | All standardized structures (RDKit), tagged with pIC50, Murcko scaffold, source-independence score, and hub_class |
| `activity_cliffs.csv` | 867 | All activity-cliff pairs (94 severe, 193 moderate, 580 broad) with Tanimoto, ΔpIC50, tier, ECFP6 Tanimoto, and `ecfp4_only_cliff` flag |
| `cliff_hubs.csv` | 4 | The four cliff-hub compounds with class labels (A = series floor, B = singleton attractor) |
| `scaffold_ruggedness_table.csv` | 375 | Every multi-member Murcko scaffold series: size, severe cliffs, cliff density, smoothness score, mean/median pIC50, potency range, mean/max SALI |
| `transformation_impact_table.csv` | — | Severe-cliff \|ΔpIC50\| by MMP transformation category (heteroatom, ring, carbon-only, halogen, aromatic) |
| `hub_physchem_table.csv` | — | Hub vs other-cliff physicochemical comparison at three degree thresholds (Mann–Whitney U + Benjamini–Hochberg FDR) |
| `hub_neighborhood_metrics.csv` | 4 | Per-hub neighborhood metrics (potency gradient, neighbor count, neighbor similarity, scaffold occupancy, local cliff count) |
| `null_model_comparison.csv` | 3 | Three-null permutation comparison (unrestricted / scaffold-constrained / assay-constrained) for cliff rarity and hub concentration |
| `source_independence_scores.csv` | 3,093 | Per-compound source list, source count, independence score, multi-source / concordance flags |

## Formats
- CSV: UTF-8, comma-separated, header row.
- SDF: RDKit V2000; molecule name = InChIKey; SD tags as listed above.

## Canonical counts (validated)
3,093 compounds · 94 severe activity cliffs · 4 hub compounds · 1,244 unique Murcko scaffolds.

See `../README_zenodo.md` for the overall package and `../code/` for the scripts that generate these files.
