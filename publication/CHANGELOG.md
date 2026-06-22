# PAD4-DB Changelog

## v2.0.0 (2026-06-16) — Current release

Initial public release of PAD4-DB v2.

### Database
- 3,093 unique PAD4 inhibitors from 95 PubChem bioassays, ChEMBL (CHEMBL6111),
  and BindingDB (Q9UM07)
- 327,336 HTS structural reference compounds from 3 PubChem screening campaigns
  (AIDs 463073, 485272, 488796)
- Six-layer architecture (A–F) with explicit mechanism class annotation

### Curation pipeline
- RDKit 2025.09.5 SMILES standardization (salt stripping, charge neutralization)
- Percent-inhibition HTS intercept raising normalization success from 89.8% to 99.0%
- InChIKey-level deduplication; log-space mean replicate aggregation
- Source independence scoring distinguishing pipeline re-curation (score 0.3–0.5)
  from genuinely independent measurements (score 0.6–1.0)

### Analysis
- ECFP4 (radius=2, 2048 bits) similarity landscape: 358,416 pairs at Tanimoto ≥ 0.6
- Activity cliff classification: 94 severe / 193 moderate / 580 broad pairs
- MMP validation via rdMMPA (maxCuts=1): 85/94 (90.4%) severe pairs confirmed
- ECFP6 (radius=3) sensitivity analysis with ecfp4_only_cliff flag
- Cliff hub discovery: 4 compounds in 2 structural classes (53.2% of severe pairs)
- Bemis-Murcko scaffold analysis: 1,244 scaffolds, Gini = 0.532

### Quality assurance
- 10-phase biological and chemical audit (86 checks, all passed)
- Reference compound recovery audit (7/13 reference inhibitors recovered)
- Fingerprint sensitivity analysis (Supplementary Table S6)
- Pre-submission stress test (17 checks, all passed)
- Full pipeline reproducibility confirmed from raw source files

### Data sources downloaded
- PubChem bioassay data: 2026-06-10 to 2026-06-14
- ChEMBL (CHEMBL6111): 2026-06-14
- BindingDB (Q9UM07): 2026-06-10

---

## v1.0.0 (internal) — Not publicly released

Internal prototype. Single-source PubChem integration only. Not peer-reviewed.
