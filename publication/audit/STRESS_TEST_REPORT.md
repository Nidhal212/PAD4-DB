# PAD4-DB v2 — Pre-Submission Stress Test Report
Generated: 2026-06-16

## Summary Table

| Code | Check | Risk | Status |
|------|-------|------|--------|
| 01A | Source independence claim | MEDIUM | ✅ Done |
| 01B | BindingDB provenance | LOW | ✅ Done |
| 01C | Download date staleness | HIGH | ✅ Done |
| 02A | Hub threshold sensitivity | MEDIUM | ✅ Done |
| 02B | SALI max artifact | LOW | ✅ Done |
| 02C | Gini interpretation | LOW | ✅ Done |
| 02D | pIC50 statistics precision | LOW | ✅ Done |
| 03A | Fingerprint robustness | MEDIUM | ✅ Done |
| 03B | Replicate aggregation | LOW | ✅ Done |
| 03C | MMP scope limitation | MEDIUM | ✅ Done |
| 03D | Scaffold variant defense | LOW | ✅ Done |
| 04A | Hub scaffold verification | HIGH | ✅ Done |
| 04B | Cross-mechanism defense | MEDIUM | ✅ Done |
| 04C | JBI-589 Ca²⁺ claim | MEDIUM | ✅ Done |
| 04D | Patent potency confound | MEDIUM | ✅ Done |
| 05A | Known compound coverage | MEDIUM | ✅ Done |
| 05B | HTS overlap framing | HIGH | ✅ Done |

## ⛔ STOP CONDITIONS TRIGGERED

The following conditions require immediate attention before submission:

- **FINGERPRINT SENSITIVITY: 64/94 severe pairs drop below Tanimoto=0.8 under ECFP6 (>20% threshold)**
- **HTS OVERLAP: 1447 of 1453 overlap compounds are HTS INACTIVE — paper framing '1,453 progressed from HTS' is wrong**

## HIGH RISK — Must address before submission

### 01C — Download Dates
ChEMBL file shows mtime 1980-01-01 (filesystem artifact). Record actual query date
in Methods. JCheminf requires explicit download dates.

### 04A — Hub Scaffold Verification
Hub A compounds claim to be in the 174-compound scaffold series.
RDKit canonicalization drift (pipeline: 174, fresh: 190) means the live check
may show different numbers. Use scaffold_family_map.csv as ground truth.

### 05B — HTS Overlap Framing
1,453 compounds appear in both HTS and dose-response datasets.
Must check whether any are HTS-inactive. Do not frame as 'HTS progression'
without verifying all are HTS actives.

## MEDIUM RISK — Should address

| Code | Action Required |
|------|----------------|
| 01A | Add explicit caveat: concordance ≠ independent replication |
| 02A | Cite Senger 2009 or Stumpfe & Bajorath 2012 for threshold justification |
| 03A | Add fingerprint sensitivity table as supplementary if >20% pairs drop |
| 03C | Clarify MMP scope: validation of cliff detection, not full SAR coverage |
| 04B | State cross-mechanism pairs are within enzymatic IC50 family |
| 04C | Cite PAD4 calcium-dependence paper for JBI-589 discrepancy explanation |
| 04D | Frame patent potency difference as optimization stage, not compound quality |
| 05A | Add limitations paragraph for post-2023 patent series |

## LOW RISK — Mention in limitations

01B, 02B, 02C, 02D, 03B, 03D

## Individual Report Files

- 01A_source_independence_defense.txt
- 01B_bindingdb_provenance.txt
- 01C_download_date_assessment.txt
- 02A_hub_threshold_sensitivity.txt
- 02B_sali_max_artifact.txt
- 02C_gini_interpretation.txt
- 02D_pic50_statistics_precision.txt
- 03A_fingerprint_robustness.txt
- 03B_replicate_aggregation_defense.txt
- 03C_mmp_scope_defense.txt
- 03D_scaffold_variant_defense.txt
- 04A_hub_scaffold_verification.txt
- 04B_cross_mechanism_defense.txt
- 04C_jbi589_discrepancy.txt
- 04D_patent_potency_confound.txt
- 05A_known_compound_coverage.txt
- 05B_hts_overlap_framing.txt