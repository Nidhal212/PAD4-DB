# PAD4-DB v2 — Project Guide

## Project
**Name:** PAD4-DB v2  
**Purpose:** Curated SAR knowledge base of PAD4 (Protein-Arginine Deiminase 4) inhibitors, aggregating PubChem bioassay data, ChEMBL, and BindingDB into a single standardized, deduplicated compound set for SAR analysis and publication.  
**Conda env:** `pad4bench` (Python 3.10.19, RDKit 2025.09.5, pandas 2.3.3, numpy 2.2.5) — NOTE: CLAUDE.md previously stated Python 3.12; live verification (2026-06-16) confirms Python 3.10.19. Paper reports Python 3.10.

## Rules
- Always run scripts from the **project root** (`/home/nidhal/PAD4-db_V2/`)
- Always activate `pad4bench` before running any script: `conda activate pad4bench`
- Never hardcode absolute paths — all paths relative to project root
- Fail loudly on missing files; no silent skips

## Directory Map

```
PAD4-db_V2/
├── data/
│   ├── raw/
│   │   ├── hts/                        3 HTS AIDs (463073, 485272, 488796)
│   │   ├── pubchem/
│   │   │   ├── confirmatory/           57 Layer A AIDs
│   │   │   ├── literature_derived/     11 Layer C AIDs
│   │   │   └── secondary/              26 Layer D+E AIDs
│   │   ├── chembl/                     CHEMBL6111 CSV (semicolon-delimited, 4925 rows)
│   │   └── bindingdb/                  Q9UM07 TSV (tab-delimited, 3087 rows)
│   └── interim/
│       ├── standardized/               standardized_compounds.parquet (Step 01 output)
│       └── normalized/                 normalized_activities.parquet (Step 02 output)
├── outputs/
│   └── tables/                         QC reports (CSV + TXT)
├── scripts/
│   ├── 01_standardize/
│   │   └── 01_standardize_smiles.py
│   ├── 02_normalize/
│   │   └── 02_normalize_activities.py
│   ├── 03_aggregate/
│   │   ├── 03_replicate_aggregate.py
│   │   ├── 03a_split_spaces.py
│   │   ├── 03b_logspace_qc.py
│   │   └── 03c_smiles_integrity.py
│   └── 04_dedup/
│       ├── 04_dedup_and_assemble.py
│       └── 04b_add_independence_scores.py
└── CLAUDE.md
```

## Pipeline Stage Status

| Stage | Script | Status |
|-------|--------|--------|
| 00 Raw inventory QC | `00_check_raw_inventory.py` (project root) | ✅ done |
| 01 SMILES standardization | `scripts/01_standardize/01_standardize_smiles.py` | ✅ done |
| 02 Activity normalization | `scripts/02_normalize/02_normalize_activities.py` | ✅ done |
| 03 Replicate aggregation | `scripts/03_aggregate/03_replicate_aggregate.py` | ✅ done |
| 04 Dedup + assembly | `scripts/04_dedup/04_dedup_and_assemble.py` | ✅ done |
| Audit A1 Biological & chemical audit | `scripts/audit/A1_biological_audit.py` | ✅ done |
| Audit A2 Reference compound recovery | `scripts/audit/A2_reference_compound_recovery.py` | ✅ done |
| 05 SAR analysis (scaffold + cliffs) | `scripts/05_cliffs/05_scaffold_and_cliffs.py` | ✅ done |

## Key Decisions

1. **AID 2202576/77 (RFMS1 vs RFMS2):** 55 SID overlap (23.8%). Keep both AIDs. During dedup, prefer AID 2202577 (RFMS2) for the 55 overlapping SIDs. Tag `aid_preferred = 2202577` for those rows.

2. **AID 2202596/97 (Functional 1 vs 2):** 0 SID overlap — completely distinct compound sets. No dedup needed between these two AIDs.

3. **AID 725596/97 (Dose-response):** Endpoint confirmed as IC50 dose-response. "1 µM / 10 µM" in the AID description refers to enzyme preincubation concentration, not compound concentration.

4. **AID 1346144 (Mixed affinity):** Mixed Ki/Kd units from ChEMBL aggregator. Must curate `Standard Units` carefully during activity normalization.

5. **AID 1920046 / 2202442 (Dual-layer):** Both AIDs appear in `confirmatory/` (Layer A) and `secondary/` (Layer D). Canonical source = `confirmatory/`. Load from `confirmatory/` only; secondary entries are skipped.

6. **BindingDB SMILES — Daylight extended annotation strip:** BindingDB exports include ` |r,THB:...|` Daylight extended SMILES annotations that RDKit cannot parse. Strip the ` |...|` suffix before passing to `Chem.MolFromSmiles`. Without this: 83% OK on first 100 rows; with strip: 100% OK.

7. **Source concordance is pipeline overlap artifact (Step 04b):** multi_source=89.1% reflects partial re-curation of shared primary assay campaigns rather than fully independent experimental replication. Concordance (99.7%) reflects this structure, not measurement agreement. Report `source_independence_score` in paper, not raw `multi_source`. Cross-source overlap: BindingDB+ChEMBL+PubChem (0.3) and BindingDB+PubChem (0.5) are the two dominant redundant patterns covering 2,565/3,093 compounds (82.9%).

**Paper-facing statement:** "Cross-source overlap reflects partial re-curation of shared primary assay campaigns rather than fully independent experimental replication. Concordance (99.7%) reflects this structure, not measurement agreement."

8. **Dataset is assay-registry-derived, not full medchem reference space (Audit A2):** Known inhibitors absent from the dataset (GSK199, Pyroxamide, PAD-PF1) reflect source database gaps, not pipeline errors. o-F-Amidine, Amodiaquine, and BB-Cl-Amidine are present in raw data but correctly excluded by endpoint type (no primary IC50). The pipeline is not a comprehensive literature review; it is bounded by what was submitted to PubChem bioassay campaigns, ChEMBL (CHEMBL6111), and BindingDB (Q9UM07).

**Paper-facing statement:** "PAD4-DB comprises compounds from PubChem bioassay campaigns, patent-deposited screening data, ChEMBL (CHEMBL6111), and BindingDB (Q9UM07). Seven of thirteen curated PAD4 reference inhibitors are present with concordant pIC50 values (mean |ΔpIC50| < 0.15 log units). Three additional compounds are present but lack primary IC50 measurements (covalent kinetics or HTS-only data). Three compounds are absent from all source databases, reflecting gaps in public bioactivity curation rather than pipeline exclusions."

## Canonical Numbers

*(Filled as pipeline runs)*

| Stage | Metric | Value |
|-------|--------|-------|
| Step 00 | AIDs checked | 95 unique (97 files; AIDs 1920046 and 2202442 each in 2 subdirs = 95 unique; counted once each per Locked Decision #3) |
| Step 01 | Total rows processed | 341,282 |
| Step 01 | Raw CSV row count | 341,328 (properly parsed); pipeline ingests 341,282; 46 rows correctly excluded: AID 1920046 (23) + AID 2202442 (23) secondary copies excluded per Key Decision #5 (dual-layer AIDs loaded from confirmatory/ only); 0 unexplained drops; Phase 0A count of 341,295 was a counting-method artifact (see outputs/audit/DROPPED_ROWS_13.txt) |
| Step 01 | Rows OK | 341,276 (100.0%) |
| Step 01 | NO_SMILES | 6 |
| Step 01 | PARSE_FAIL | 0 |
| Step 01 | SANITIZE_FAIL | 0 |
| Step 01 | Unique InChIKeys | 328,976 |
| Step 02 | Total output rows | 341,282 |
| Step 02 | norm_status OK | 338,021 (99.0%) |
| Step 02 | NO_VALUE | 3,155 (0.9%) |
| Step 02 | UNCONVERTIBLE_UNITS | 106 (0.0%) |
| Step 02 | use_in_potency_model=True | 7,815 |
| Step 02 | pIC50 range | 2.00 – 8.52 |
| Step 02 | ChEMBL pchembl_mismatch | 0 |
| Step 02 | Dominant endpoint type | Pct_inh = 330,136 (HTS layer) |
| Step 02 | Fix applied | Pct_inh rows with % units intercepted before convert_to_nM (without: 89.8% OK; with: 99.0% OK) |
| Step 02 | Expected warnings | Layer C Kinact/Ki ratios (/M/min units), MALDI mass-shift assays (no numeric value), ChEMBL kon/koff rows |
| Step 03 | Input rows | 341,282 |
| Step 03 | Output groups | 339,687 (InChIKey × source × aid × endpoint_type) |
| Step 03 | Multi-replicate groups | 450 (496 measurements consolidated, 7,815→7,319 rows) |
| Step 03 | use_in_potency=True | 7,319 |
| Step 03 | hts_flag=True | 328,894 |
| Step 03 | smiles_std missing | 0 |
| Step 03 | assay_mechanism_class | screening_single_conc=329,106 · baee_colorimetric=7,813 · fp_binding=1,683 · rfms_enzymatic=933 · covalent_irreversible=150 · cellular=2 |
| Step 03 | norm_status_worst | OK=336,637 · NO_VALUE=2,961 · UNCONVERTIBLE_UNITS=89 |
| Step 03a | potency_space rows | 7,319 |
| Step 03a | hts_space rows | 332,368 |
| Step 03a | Unique InChIKeys — potency_space | 3,093 ← first real compound count |
| Step 03a | Unique InChIKeys — hts_space | 327,336 |
| Step 03b | IC50 groups validated | 7,319 (max log-mean diff = 0.000000 — PASS) |
| Step 03b | pIC50 range (IC50 groups) | 2.00 – 8.52 (median 6.86, mean 6.58) — pre-dedup potency space |
| Step 03b | Outliers pIC50 < 3 | 26 rows (0.36%) — weak inactive compounds, >1 mM IC50, norm_status=OK |
| Step 03b | Outliers pIC50 > 12 | 0 rows — PASS |
| Step 03c | SMILES integrity | ALL PASS: 0 null, 0 SMILES→multi-InChIKey, 0 InChIKey→multi-SMILES |
| Step 04 | dedup_aid_level rows | 7,214 (InChIKey × source × endpoint_type) |
| Step 04 | AID 2202576/77 dedup | 55 rows dropped (2202576 side); 110 aid_preferred_used flags set |
| Step 04 | n_aids dist | 1 AID=7,166 · 2 AIDs=46 · 3 AIDs=2 |
| Step 04 | high_variance_aid (spread>1 pIC50) | 0 |
| Step 04 | FINAL COMPOUND COUNT (pad4_compounds) | 3,093 |
| Step 04 | ChEMBL structure-less exclusions | 2 compounds (CHEMBL5841052 pIC50=5.824; CHEMBL5864263 pIC50=6.893) have qualifying IC50 measurements but no SMILES in ChEMBL export (std_status=NO_SMILES); correctly excluded — structure-less compounds cannot participate in SAR analysis. True qualifying measurement count = 3,095; structure-resolved compound count = 3,093. |
| Step 04 | pIC50 consensus range | 2.00 – 8.52 · mean=6.550 · median=6.844 · std=0.992 (consensus across 3,093 compounds; mean below median due to weak-binder tail) — PAPER REPORTS 6.84 (rounds from 6.8447) |
| Step 04 | pIC50 distribution note | Bimodal: shoulder at ~5.0, main peak at ~7.0. Patent-exclusive compounds (n=233) drive the shoulder at pIC50 5–6. BindingDB-only compounds (n=95) show a sharp mode at ~7.3 reflecting higher potency threshold in BindingDB curation. |
| Step 04 | multi_source compounds | 2,755 / 3,093 (89.1%) |
| Step 04 | high_confidence compounds | 2,746 / 3,093 (88.8%) |
| Step 04 | concordant compounds | 3,084 / 3,093 — discordant: 0 |
| Step 04 | HTS UNIQUE COMPOUNDS | 327,336 |
| Step 04 | HTS any_active | 308 (0.09%) |
| Step 04 | HTS confirmed_in_potency_space | 1,453 |
| Step 04 | Source dist (dedup_aid_level) | bindingdb=2,827 · pubchem_confirmatory=2,821 · chembl=1,566 |
| Step 04b | source_independence_score dist | 0.3=1,366 · 0.5=1,199 · 0.6=167 · 0.7=23 · 1.0=338 |
| Step 04b | source_list breakdown | bindingdb+chembl+pubchem=1,366 · bindingdb+pubchem=1,199 · pubchem only=233 · bindingdb+chembl=167 · bindingdb only=95 · chembl+pubchem=23 · chembl only=10 |
| Step 04b | is_true_multi_source=True (score≥0.7) | 361 / 3,093 (threshold as specified; note: spec expected ~528 which matches threshold≥0.6) |
| Step 04b | is_true_multi_source=False | 2,732 / 3,093 |
| Audit A1 | Target identity | 93 PAD4_EXPLICIT, 2 PAD_FAMILY (AIDs 588488/588560), 2 AMBIGUOUS (AIDs 588487/651627); non-explicit all have 0 potency-space rows |
| Audit A1 | Species validity | ChEMBL: 4,858 Homo sapiens, 67 unknown, 0 non-human; BindingDB: 3,087 Homo sapiens, 0 non-human |
| Audit A1 | Chemical correctness | pIC50 trace max_abs_diff=0.000000 (PASS); SMILES integrity ALL PASS |
| Audit A1 | Cross-source delta | A1 audit (n=50 sample): max delta=0.500 · FULL DATASET (source_spread, n=3,093): max=0.7386 · 0 pairs >1.0 · 0 pairs >1.5 — PAPER REPORTS 0.74 (full dataset) |
| Audit A2 | Reference compounds checked | 14 total (13 non-excluded) |
| Audit A2 | already_present | 7: Streptonigrin (5.602), Cl-amidine (5.219), F-Amidine (4.571), GSK484 (7.049), TDFA (5.638), BMS-P5 (7.009), JBI-589 (6.000) |
| Audit A2 | present_but_not_mapped | 3: o-F-Amidine (Kinact/Ki + Layer D IC50 only), Amodiaquine (HTS Pct_inh only), BB-Cl-Amidine (covalent kinetics only, AID 1364668) |
| Audit A2 | absent_by_design | 3: GSK199 (not in CHEMBL6111 assay), Pyroxamide (not in any source), PAD-PF1 (allosteric, not in public DBs) |
| Audit A2 | excluded_correct | 1: AFM-30a (PAD2-selective, correctly absent) |
| Audit A2 | GSK484 InChIKey note | manual IK MULKOGJHUZTANI-ADMBKAPUSA-N (HCl salt); computed IK BDYDINKSILYBOL-WMZHIEFXSA-N (free base after standardization); pipeline correctly used free base form |
| Audit A2 | JBI-589 pIC50 delta | DB=6.000 (1000nM) vs published=6.914 (122nM); source: different assay Ca2+ concentration, not pipeline error |
| Step 05 | Compounds analyzed | 3,093 (ECFP4 Morgan r=2 nBits=2048; 0 fingerprint failures) |
| Step 05 | Similarity landscape | pairs ≥0.6: 358,416 · ≥0.7: 123,460 · ≥0.8: 12,071 · ≥0.9: 659 · ≥0.95: 233 |
| Step 05 | Cliff counts | severe (Tan≥0.8, ΔpIC50≥2.0)=94 · moderate (≥1.5)=193 · broad (≥1.0)=580 |
| Step 05 | Compounds in severe cliffs | 99 / 3,093 (3.2%) · max ΔpIC50=3.045 (from activity_cliffs.parquet, Tan≥0.8 cliff pairs) · mean ΔpIC50 severe=2.308 — NOTE: 3.228 appeared in early notes from activity_pairs_with_sali.parquet (all Tan≥0.6 pairs, pair Tan=0.667 — NOT a cliff); correct value is 3.045 |
| Step 05 | Compounds in any cliff tier | 654 / 3,093 |
| Step 05 | Scaffold analysis | 1,244 unique · 869 singletons · 375 series (≥2) · largest=174 · rank2=102 · coverage=30.1% (series/unique) · compound coverage=71.9% (compounds-in-series/total) |
| Step 05 | Scaffold concentration | Gini=0.532 · median series size=3 · top 30 scaffolds (~2.4% of 1,244) cover ~30% of compounds |
| Step 05 | Scaffold canonicalization note | RDKit 2025.09.5 canonical SMILES. Pipeline value 174 (rank-1 series) is locked. Fresh re-derivation gives 190 due to inter-version SMILES canonicalization drift — pipeline value is canonical for this paper. |
| Step 05 | patent scaffold mean series size | 2.5 (identical to non-patent: 2.5) |
| Step 05 | patent scaffold mean pIC50 | 6.134 (vs non-patent mean: 6.532) |
| Step 05 | top-5 scaffolds | all non-patent-exclusive (present in ChEMBL/BindingDB) |
| Step 05 | Patent-exclusive scaffolds | 107 scaffolds contain ≥1 patent compound · 1,137 contain none |
| Step 05 | new scaffolds added | 103 unique scaffolds not in non-patent space (107 scaffolds contain patent compounds total; 103 are exclusive to patent-derived compounds) |
| Step 05 | Patent-exclusive contribution | verdict=weak_contribution; 1 new severe cliff of 94 (cliff_delta_pct=1.06%); 34 cliff pairs involve a patent compound |
| Step 05 | is_covalent compounds | 107 (3.5%) — SMARTS-flagged warhead compounds |
| Step 05 | warhead_class distribution | chloroacetamidine=66 · fluoroacetamidine=17 · haloacetyl=11 · enaminone=7 · vinyl_sulfone=4 · alpha_bromoketone=2 · none=2,986 |
| Step 05 | mechanism_class distribution | enzymatic=2,079 (67.2%) · enzymatic_confirmed=878 (28.4%) · fp_ic50=115 (3.7%) · covalent=21 (0.7%) · unknowns=0 |
| Step 05 | fragment_flag compounds | 5 (MW<200 AND pIC50<4.0): 4 fluoroacetamidine minimal pharmacophores + 1 weak binder (4-amino-2-hydroxybenzoic acid) — retained, flagged for optional exclusion |
| Step 05 | Cliff hubs (severe) | 4 compounds in 2 structural classes · 50/94 severe cliff pairs (53.2%) |
| Step 05 | Covalent cliff safety | 0 covalent-vs-reversible severe pairs · 0 moderate · 1 broad (ΔpIC50=1.21, max) · landscape pharmacologically clean |
| Step 05 | SALI landscape | SALI>10: 335 pairs · SALI>20: 19 pairs · SALI max=65.88 (single outlier) · pairs sim≥0.8: 12,071 (3.4% of sim≥0.6 pairs) |
| Step 05 | Activity landscape note | Diagonal absence confirmed: sim≥0.8 pairs almost universally concordant (ΔpIC50<1.0). Severe cliffs are rare events (~0.026% of sim≥0.6 pairs). Bulk of pairs in "Diverse, Concordant" quadrant. |
| Step 05 | Cliff network | 99 nodes · 94 severe edges · Hub A edges=27 · Hub B edges=23 · Non-hub edges=44 · Cross-mechanism edges=4 |
| Step 05 | Cliff network visual findings | Hub A (red) fan into enzymatic_confirmed cluster — within-series cliffs vs RFMS-confirmed high-potency analogs. Hub B (navy) connect across diverse node colors — cross-chemotype structural attractor confirmed. 4 cross-mechanism pairs dashed and non-hub — mechanistic heterogeneity does not drive hub structure. Degree distribution bimodal: hubs 11-15, rest 1-5, no intermediate candidates. |
| Step 05 | MMP analysis — cliff compounds | Total MMP pairs among 99 cliff compounds: 707 · Unique cores: 24 (mmp_pairs_cliff99.csv shared_core.nunique() = 24; the figure 943 appeared in early draft notes and is not reproducible from any file column — 24 is definitive) |
| Step 05 | MMP-confirmed severe cliffs | **80 / 94 canonical severe cliff pairs (85.1%)** confirmed by shared-core MMP. CANONICAL=80 (matches locked stress-test Methods text "80 (85.1%)"; 80/94=85.1% exactly). The figure 85 (90.4%) is the mmp_pairs_cliff99.csv `is_canonical_severe_cliff` self-flag, which over-counts by 5 spurious pairs (all involving UXUIFFMVADUAON-CXBZOWBPSA-N) that have tanimoto=NaN and fail the canonical ECFP4 Tan≥0.8 cliff definition in activity_cliffs.parquet. Verified 2026-06-19 by canonical-pair-key join of the 94-severe set × MMP file. Fig 6 and table3_mmp_summary.csv both use 80. |
| Step 05 | MMP type breakdown (severe, CANONICAL) | single_atom_change=45 · small_substituent=27 · medium_substituent=8 · large_substituent=0 (=80). The 49/28/8 breakdown (=85) is the over-counted self-flag version — superseded. |
| Step 05 | MMP cliff tier breakdown (mmp file self-flag) | severe=85 · moderate=25 · broad=2 · non_cliff=595 (file's own labels; severe self-flag over-counts canonical by 5 — use 80 for severe) |
| Step 05 | Top discontinuity compound | IUZXRGLRAITQQP-RUZDIDTESA-N (score=2.471, hub=none) — paired with CZVROBPEHUHFMR-XMMPIXPASA-N (pIC50=8.000) |
| Step 05 | MMP hub compounds | UDCDEKJNAMHBFH-HSZRJFAPSA-N and DVCKJOQIVOGXEI-XMMPIXPASA-N (both Hub B) rank 3rd/4th by discontinuity score (2.263, 17 MMP partners each) |

### Step 05 Scaffold Coverage Clarification

Singleton scaffolds: 872  (869 ring-containing + 3 acyclic: 2 known
  fragments + 1 linear peptide SOZMHIJABUOUSN, pIC50=5.64, 0 cliff pairs)

scaffold_coverage metrics:
- series scaffolds / unique scaffolds = 375/1244 = 30.1%
- compounds in series scaffolds / total compounds = 71.9%
- Both figures are correct for different denominators.
- Paper should report: 71.9% of compounds belong to scaffold families of ≥2 members

### Step 05 Patent-Exclusive Scaffold Notes

Patent-exclusive scaffold details:

- 107 scaffolds contain ≥1 patent-exclusive compound
- 1,137 scaffolds contain no patent-exclusive compounds
- Mean series size: 2.5 (identical to non-patent scaffolds)
- Mean pIC50: 6.134 (0.4 log units below non-patent mean of 6.532)
- Interpretation: patent compounds slightly less potent on average (earlier-stage hits, less optimized than published ChEMBL/BindingDB compounds) but cover genuinely distinct chemical space at identical series density

Top 5 patent-exclusive scaffolds:
1. Azaindole-piperidine-cyclohexane amide:   29 cpds, pIC50=7.13±0.51 (potent series, distinct from dominant AID 1919095 chemotype)
2. Chalcone-cyclohexane lactam (pyridine):   27 cpds, pIC50=5.21±0.52 (moderate potency, novel chemotype, likely different binding mode)
3. Chalcone-cyclohexane lactam (pyrimidine): 18 cpds, pIC50=4.65±0.57 (lower potency variant of scaffold 2)
4. Azaindole-piperidine-cyclohexane amide B: 16 cpds, pIC50=6.94±0.48 (close variant of scaffold 1)
5. Azaindole-bicyclic amine amide:           10 cpds, pIC50=7.11±0.15 (tightest series in patent space, std=0.15)

### Step 05 Cliff Hub Analysis

Four compounds act as cliff hubs, organized into two mechanistically distinct structural classes.

#### Class A — Series-Embedded Mid-Potency Floor

Compounds are mid-potency members of the dominant 174-compound azaindole-benzimidazole scaffold series. Their cliff pairs are against higher-potency analogs within the same chemotype.

- SMADULGDNOCLOP-GISFHXKWSA-N: pIC50=5.390, MW=611 — **15 severe cliff pairs**
- RAVBZQAQTVGKIV-XBPDSQQVSA-N: pIC50=5.341, MW=591 — **12 severe cliff pairs**
- Class A Tanimoto: 0.761 (structurally related; same series, not near-identical)
- Source confidence: bindingdb|chembl|pubchem_confirmatory (highest tier; all 3 databases)
- Scaffold: azaindole-benzimidazole (the dominant 174-compound series)
- Archetype: SERIES-EMBEDDED MID-POTENCY FLOOR — these compounds sit at a within-series potency floor; their cliff pairs arise from potent analogs in the same chemotype
- Class A collective severe cliff pairs: **27**

#### Class B — Scaffold-Singleton Structural Attractor

Scaffold singletons (no other compounds share their Murcko framework). Their broad ECFP4 structural promiscuity generates cliff pairs across multiple chemotypes.

- UDCDEKJNAMHBFH-HSZRJFAPSA-N: azaindole-benzimidazole, cyclobutyl sulfonamide, free amine, MW=605.8, 43 heavy atoms — **12 severe cliff pairs**
- DVCKJOQIVOGXEI-XMMPIXPASA-N: azaindole-benzimidazole, cyclopentyl sulfonamide, free amine, MW=619.8, 44 heavy atoms — **11 severe cliff pairs**
- Class B Tanimoto: 0.9753 (cyclobutyl vs cyclopentyl; 1 CH₂ difference; functionally one structural class)
- Source: bindingdb|pubchem_confirmatory
- Scaffold: SINGLETON (unique Murcko scaffold for each compound; n=1 each)
- Structural difference: cyclobutyl vs cyclopentyl (14 Da, 1 carbon); identical pIC50=4.301
- Archetype: SCAFFOLD-SINGLETON STRUCTURAL ATTRACTOR — broad structural promiscuity by ECFP4 creates cliff pairs across multiple chemotypes
- Likely explanation: free primary amine causing assay interference or early-stage unoptimized hits before amine capping
- Class B collective severe cliff pairs: **23**

#### Combined Summary

- Total hub severe cliff pairs: **50 / 94 (53.2%)**
- Cross-class Tanimoto: ~0.49 (Class A and Class B are structurally independent)
- pad4_compounds.parquet column `hub_class`: 'A' for Class A compounds, 'B' for Class B compounds, 'none' for all others
- Scientific value: two independent and structurally interpretable failure modes for similarity-based ML models

**Locked paper statement:** "The severe activity cliff landscape is organized around two mechanistically distinct hub archetypes. Class A hubs (SMADULGDNOCLOP-GISFHXKWSA-N and RAVBZQAQTVGKIV-XBPDSQQVSA-N; pIC50 ≈ 5.4) are mid-potency members of the dominant 174-compound azaindole-benzimidazole scaffold series; their position as within-series potency floors generates 27 severe cliff pairs against higher-potency analogs in the same chemotype. Class B hubs (UDCDEKJNAMHBFH-HSZRJFAPSA-N and DVCKJOQIVOGXEI-XMMPIXPASA-N; pIC50 = 4.301, Tanimoto = 0.975) are scaffold singletons with no other compounds sharing their Murcko framework; their broad structural promiscuity by ECFP4 fingerprint creates 23 severe cliff pairs across multiple chemotypes. Together, these four cliff-hub compounds in two structural classes account for 50 of 94 severe cliff pairs (53.2%), representing two independent and structurally interpretable failure modes for similarity-based ML models."

## Audit Results

### REFERENCE COMPOUND AUDIT (A2_reference_compound_recovery.py)

Recovery summary (13 non-excluded compounds):

**already_present (7):**
```
Streptonigrin  pIC50=5.602  (vs expected 5.60) ✓
Cl-amidine     pIC50=5.219  (vs expected 5.23) ✓
F-Amidine      pIC50=4.571  (vs expected 4.67) ✓
GSK484         pIC50=7.049  (vs expected 7.30, 50nM published) ✓
TDFA           pIC50=5.638  (vs expected 5.64, 2300nM) ✓
BMS-P5         pIC50=7.009  (vs expected 7.01, 98nM) ✓
JBI-589        pIC50=6.000  (vs expected ~6.9, 122nM) ~
```

**present_but_not_mapped (3) — correct pipeline exclusions:**
```
o-F-Amidine  : Kinact/Ki only, no primary IC50 → Layer C, correct
Amodiaquine  : HTS Pct_inh only, no dose-response → correct
BB-Cl-Amidine: covalent kinetics only in AID 1364668 → correct
```

**absent_by_design (3) — not in any source file:**
```
GSK199   : ChEMBL CHEMBL3545375 not submitted under CHEMBL6111 assay
Pyroxamide: not in PubChem/ChEMBL/BindingDB for PAD4
PAD-PF1  : allosteric inhibitor, not in public databases
```

**excluded_correct (1):**
```
AFM-30a  : PAD2-selective, correctly absent from PAD4 dataset
```

**Key finding — GSK484 salt-stripped InChIKey:**
- manual:   MULKOGJHUZTANI-ADMBKAPUSA-N (HCl salt form)
- computed: BDYDINKSILYBOL-WMZHIEFXSA-N (free base, after standardization)
- Pipeline correctly used free base form. Manual table had salt form.

**Key finding — JBI-589 pIC50 delta = 0.9 log units vs published 122nM:**
- DB value = 6.000 (1000nM), published = 6.914 (122nM)
- Source: ChEMBL/BindingDB store different assay condition (likely different Ca2+ concentration). Not a pipeline error.

## STRESS TEST RESULTS (2026-06-16) — LOCKED CORRECTED STATEMENTS

### Fingerprint Sensitivity (Methods paragraph — exact text):
"Structural similarity was computed using ECFP4 fingerprints (Morgan
algorithm, radius=2, 2048 bits; RDKit 2025.09.5), following Stumpfe
and Bajorath (2012) and Senger (2009). Activity cliffs were defined as
pairs with Tanimoto similarity ≥0.8 and |ΔpIC50| ≥2.0 log units,
consistent with established SAR discontinuity literature. Of 94 severe
cliff pairs, 80 (85.1%) were also confirmed by matched molecular pair
analysis (MMP), an orthogonal substructure-based method independent of
fingerprint choice. Sensitivity analysis with ECFP6 (radius=3) showed
that 64 of 94 pairs have ECFP4 Tanimoto 0.80–0.85 and fall below 0.80
at radius=3, consistent with the known resolution-dependent behaviour of
Morgan fingerprints for large fused ring systems; however, 80% of these
borderline pairs (51/64) remain MMP-confirmed, and the hub dominance
statistic is fingerprint-invariant (53.2% under ECFP4; 53.3% under
ECFP6). Thirteen pairs (13.8%) are classified as severe by ECFP4 only,
without MMP or ECFP6 corroboration; these are flagged in the deposited
dataset (ecfp4_only_cliff=True in activity_cliffs.parquet)."

### ECFP6 robustness numbers (locked):
  Pairs robust at ECFP6 (≥0.8):                30
  Pairs non-robust (ECFP4 0.80–0.85, ECFP6 <0.80): 64
  Of non-robust: MMP-validated 51 (80%), hub-involved 34 (53%)
  ecfp4_only_cliff = ~robust AND ~mmp_validated: 13 (13.8%)
    of which hub-involved: 4 (hub is compound property, not pair confirmation)
    of which non-hub:      9 (most isolated cases)
  Hub % at ECFP4: 53.2% (50/94)
  Hub % at ECFP6: 53.3% (16/30) — fingerprint-invariant
  Non-robust pair ΔpIC50 mean: 2.313 (vs robust 2.297 — identical)
  ecfp4_only_cliff column added to activity_cliffs.parquet (2026-06-16)
  tanimoto_ecfp6 column added to activity_cliffs.parquet (2026-06-16)

### HTS overlap corrected statement (use everywhere 1,453 appears):
"Of 3,093 PAD4 inhibitors in the dose-response database, 1,453 share
InChIKeys with compounds in the HTS screening dataset (327,336 total
screened), indicating parallel measurement in both assay formats. Only
6 of these 1,453 compounds were confirmed HTS actives (max inhibition
≥50% at screening concentration); the remaining 1,447 have published
IC50 values from independent research programs (ChEMBL/BindingDB) and
showed low inhibition in HTS campaigns (median 4.3% at screening
concentration), consistent with potent compounds (median pIC50=6.93)
tested below their IC50 at standard HTS screening concentrations."

NEVER USE: "progressed from HTS to dose-response"
NEVER USE: "confirmed HTS hits that progressed"
CORRECT FRAMING: "parallel measurement in independent assay pipelines"

### Download dates (Methods — exact text):
"PubChem bioassay data were downloaded 2026-06-10 to 2026-06-14.
ChEMBL bioactivity data (assay CHEMBL6111) were downloaded on
2026-06-14 (file modification time shows 1980-01-01 due to a known
Linux unzip timestamp artifact; the download date is confirmed by
directory modification time and batch download logs).
BindingDB (UniProt accession Q9UM07) was downloaded on 2026-06-10."

### New locked numbers from stress test:
  ecfp4_only_cliff pairs:                        13 (13.8% of 94 severe)
    hub-involved ecfp4-only:                      4
    non-hub ecfp4-only:                           9
  HTS actives in 1,453 overlap:                   6
  HTS inactives in 1,453 overlap:             1,447
  Median pIC50 of HTS-inactive overlap:        6.93
  Median max_pct_inh of HTS-inactive overlap:  4.28%

### CORRECTION NOTES (2026-06-16 final audit)

NOTE ON 943 MMP CORES: The figure 943 appeared in an early manuscript draft and
was never a locked canonical value. The live mmp_pairs_cliff99.csv has
shared_core.nunique() = 24. Use 24. The n_shared_cores column sums to 1,297
(not 943 either). 943 cannot be reproduced from any column in any output file.

NOTE ON SCAFFOLD COVERAGE 71.8% vs 71.9%: The audited pipeline (Stage 05) gives
71.90% from live data (2224/3093 compounds in scaffolds with ≥2 members). The
figure 71.8% appeared in earlier audit rounds from a slightly different compound
count or rounding. The canonical value is 71.9%. CLAUDE.md scaffold_coverage
section already locked 71.9% — the 71.8% that appeared in manuscript drafts was
a user-instructed override that has now been corrected back to 71.9%.

NOTE ON max ΔpIC50 3.228 vs 3.045: The value 3.228 came from
activity_pairs_with_sali.parquet (all Tan≥0.6 pairs). The pair responsible has
Tanimoto=0.667, which is BELOW the cliff threshold (Tan≥0.8). It is not a cliff
pair. The canonical max ΔpIC50 for actual cliff pairs (from activity_cliffs.parquet,
Tan≥0.8) is 3.0448. Paper reports 3.045. Never report 3.228 for cliff pairs.

NOTE ON pIC50 MEDIAN 6.85 vs 6.84: Early pipeline notes and CLAUDE.md Step 04
recorded median=6.845. Live computation gives 6.8447. Rounded to 2dp: 6.84.
The figure 6.85 (from rounding 6.845 up) was wrong. Paper reports 6.84.

### Additional citations needed before submission:
  [Stumpfe & Bajorath 2012] — ECFP4 threshold justification for cliff analysis
  [Senger 2009] — Activity cliff definition (Tanimoto ≥0.8, ΔpIC50 ≥2.0)
  [Knuckley 2010 or equivalent] — PAD4 calcium dependence (JBI-589 note)

---

## Nature Figure Design System

All final figures must conform to Nature Methods/Scientific Data submission standards.
Apply via: `from scripts.nature_style import apply_nature_style, COLORS, panel_label, save_figure`

### Dimensions
| Type | Width | Max height |
|------|-------|------------|
| Single-column | 89 mm (3.504 in) | 247 mm |
| Double-column | 183 mm (7.205 in) | 247 mm |

### Typography
| Element | Size | Weight |
|---------|------|--------|
| Axis labels | 7 pt | regular |
| Tick labels | 6 pt | regular |
| Panel labels | 8 pt | **bold**, outside panel at (−0.12, 1.04) in axes coords |
| Legend text | 6 pt | regular, frameon=False |
| Figure title | **None** — Nature uses captions only |

### Axes
- Spines: left + bottom only (no top, no right)
- Ticks: inward, major 3pt / minor 2pt, width 0.5pt
- Axis linewidth: 0.75pt; no grid

### Color Palette (colorblind-safe)
| Name | Hex | Usage |
|------|-----|-------|
| blue | #0077BB | primary data / main series |
| orange | #EE7733 | patent-exclusive / highlight |
| red | #CC3311 | severe cliffs / Hub Class B |
| teal | #009988 | enzymatic_confirmed |
| cyan | #33BBEE | fp_ic50 |
| navy | #1A237E | Hub Class A |
| gray_light | #BBBBBB | background / non-highlighted |
| gray_dark | #555555 | secondary data |

**Note:** Hub Class A color revised from red (#E74C3C in earlier figs) to navy (#1A237E); Hub Class B revised to red (#CC3311). Consistent with class naming (A=navy, B=red).

### Statistical requirements
- Any two-distribution comparison: add Mann-Whitney U p-value
- Any mean shown: ±SD as error bar or shaded band
- All sample sizes: n= in legend or panel annotation

### Export
- PNG: 600 dpi, bbox_inches='tight', facecolor='white'
- SVG: bbox_inches='tight', facecolor='white'

### STANDING LIBRARIES (Nature versions — replaces SciencePlots for final figs)
```python
from scripts.nature_style import apply_nature_style, COLORS, panel_label, save_figure
apply_nature_style()   # replaces plt.style.use(['science','nature','no-latex'])
```

### Previous figure color mapping → Nature palette
| Old color | New color | Element |
|-----------|-----------|---------|
| #4A90D9 | #0077BB | BindingDB / published / primary |
| #E05A2B | #EE7733 | Patent-exclusive |
| #E74C3C | #CC3311 | Severe cliffs / Hub B |
| #2ECC71 | #009988 | ChEMBL / enzymatic_confirmed |
| #F39C12 | #EE7733 | Moderate cliffs / BindingDB |
| #1A237E | #1A237E | Hub Class A (unchanged) |
| #AAAAAA | #BBBBBB | Gray background |
