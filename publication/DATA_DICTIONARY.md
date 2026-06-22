# PAD4-DB v2 — Data Dictionary

Complete column descriptions for all deposited data files.

---

## pad4_compounds.parquet

Primary compound file. 3,093 rows × 25 columns.

| Column | Type | Description |
|--------|------|-------------|
| `inchi_key` | str | Standard InChIKey (IUPAC International Chemical Identifier key); primary compound identifier |
| `smiles_std` | str | Canonical SMILES after RDKit 2025.09.5 standardization (salt-stripped, charge-neutralized) |
| `pic50_consensus` | float | Consensus pIC50 = −log₁₀(IC50[M]); mean of per-assay median pIC50 values across all Layer A sources |
| `source_list` | str | Pipe-delimited list of source databases contributing IC50 data (e.g. `bindingdb\|chembl\|pubchem_confirmatory`) |
| `source_independence_score` | float | Independence score: 0.3 (BDB+ChEMBL+PC), 0.5 (BDB+PC), 0.6 (BDB+ChEMBL), 0.7 (ChEMBL+PC), 1.0 (single-source). Use ≥0.6 for genuinely independent multi-source data. |
| `multi_source` | bool | True if compound appears in ≥2 source databases (89.1% of dataset; reflects pipeline overlap, not independent replication) |
| `concordant` | bool | True if max cross-source pIC50 delta ≤ 1.0 log unit (99.7% of dataset) |
| `high_confidence` | bool | True if source_independence_score ≥ 0.3 AND concordant (88.8% of dataset) |
| `source_spread` | float | Maximum pIC50 difference across sources for this compound; 0 for single-source compounds |
| `log_value_std_global` | float | Standard deviation of all pIC50 values across all assays for this compound; measure of within-compound variability |
| `mechanism_class` | str | Assay mechanism class: `enzymatic` (BAEE colorimetric, n=2,079), `enzymatic_confirmed` (RFMS fluorescence, n=878), `fp_ic50` (fluorescence polarization, n=115), `covalent` (irreversible kinetics, n=21) |
| `is_covalent` | bool | True if compound contains at least one SMARTS-flagged reactive warhead (n=107, 3.5%) |
| `warhead_class` | str | Chemical class of covalent warhead: `chloroacetamidine`, `fluoroacetamidine`, `haloacetyl`, `enaminone`, `vinyl_sulfone`, `alpha_bromoketone`, or None |
| `fragment_flag` | bool | True if MW < 200 Da AND pIC50 < 4.0 (n=5; retained but flagged for optional exclusion) |
| `scaffold_smiles` | str | Bemis-Murcko heteroatom-preserving scaffold SMILES (RDKit 2025.09.5 canonical) |
| `scaffold_family_size` | int | Number of compounds sharing this scaffold in the dataset |
| `patent_exclusive` | bool | True if compound is present only in PubChem confirmatory data without ChEMBL or BindingDB coverage (n=233) |
| `hub_class` | str | Cliff-hub structural archetype: `A` (series-embedded potency floor, 2 compounds), `B` (scaffold-singleton attractor, 2 compounds), or `none` (3,089 compounds) |
| `cliff_tier_max` | str | Highest cliff tier this compound participates in: `severe`, `moderate`, `broad`, or `none` |
| `ecfp4_only_cliff` | bool | True if compound participates in at least one ECFP4-only severe cliff pair (not confirmed by ECFP6 ≥0.8 or MMP) |

---

## activity_cliffs.parquet

All cliff pairs at three tiers. 867 rows × 13 columns.

| Column | Type | Description |
|--------|------|-------------|
| `inchi_key_a` | str | InChIKey of compound A (lower pIC50 of the pair) |
| `inchi_key_b` | str | InChIKey of compound B (higher pIC50 of the pair) |
| `tanimoto` | float | ECFP4 Tanimoto similarity (Morgan radius=2, 2048 bits, RDKit 2025.09.5) |
| `tanimoto_ecfp6` | float | ECFP6 Tanimoto similarity (Morgan radius=3, 2048 bits) for fingerprint sensitivity analysis |
| `pic50_a` | float | Consensus pIC50 of compound A |
| `pic50_b` | float | Consensus pIC50 of compound B |
| `delta_pic50` | float | pIC50_b − pIC50_a (always positive; compound B is more potent) |
| `cliff_tier` | str | `severe` (Tan≥0.8, ΔpIC50≥2.0), `moderate` (Tan≥0.8, ΔpIC50≥1.5), or `broad` (Tan≥0.8, ΔpIC50≥1.0) |
| `source_combination` | str | Combined source annotations for the pair |
| `patent_exclusive_a` | bool | True if compound A is patent-exclusive |
| `patent_exclusive_b` | bool | True if compound B is patent-exclusive |
| `any_patent_exclusive` | bool | True if either compound is patent-exclusive |
| `ecfp4_only_cliff` | bool | True if this severe pair is confirmed by ECFP4 only (not by ECFP6 ≥0.8 or MMP); applicable to severe tier only |

**Cliff tier counts:** severe=94, moderate=193, broad=580

---

## activity_pairs_with_sali.parquet

All pairwise comparisons with Tanimoto ≥ 0.6. 358,416 rows × 10 columns.

| Column | Type | Description |
|--------|------|-------------|
| `inchi_key_a` | str | InChIKey of compound A |
| `inchi_key_b` | str | InChIKey of compound B |
| `tanimoto` | float | ECFP4 Tanimoto similarity |
| `pic50_a` | float | Consensus pIC50 of compound A |
| `pic50_b` | float | Consensus pIC50 of compound B |
| `delta_pic50` | float | Absolute pIC50 difference |
| `source_a` | str | Source database(s) for compound A |
| `source_b` | str | Source database(s) for compound B |
| `source_combination` | str | Combined source label for the pair |
| `same_source` | bool | True if both compounds come from the same source database |

**Note:** SALI (Structure-Activity Landscape Index) = ΔpIC50 / (1 − Tanimoto).
For the activity cliff analysis, only pairs with Tanimoto ≥ 0.8 are used as cliffs.

---

## hts_compound_index.parquet

HTS structural reference compounds. 327,336 rows × 8 columns.

| Column | Type | Description |
|--------|------|-------------|
| `inchi_key` | str | Standard InChIKey |
| `smiles_std` | str | Canonical SMILES |
| `any_active` | bool | True if compound showed ≥50% inhibition in any HTS campaign (n=308, 0.09%) |
| `hts_outcome` | str | `Active` or `Inactive` based on max_pct_inh threshold |
| `hts_consensus_confidence` | float | Confidence score based on replication across HTS campaigns (0–1) |
| `max_pct_inh` | float | Maximum percent inhibition observed across all HTS assays |
| `confirmed_in_potency_space` | bool | True if compound also appears in the SAR dose-response dataset (n=1,453) |
| `n_hts_aids` | int | Number of HTS AIDs in which this compound was tested |

**Note on 1,453 overlap compounds:** These share InChIKeys between the HTS structural
reference and the SAR dose-response dataset. Only 6 are confirmed HTS actives;
the remaining 1,447 have IC50 values from independent ChEMBL/BindingDB sources
and showed low HTS inhibition (median 4.3%), consistent with potent compounds
tested below their IC50 at screening concentration.

---

## mmp_pairs_cliff99.csv

Matched molecular pairs among the 99 severe cliff compounds. 707 rows × 12 columns.

| Column | Type | Description |
|--------|------|-------------|
| `inchi_key_a` | str | InChIKey of compound A |
| `inchi_key_b` | str | InChIKey of compound B |
| `shared_core` | str | SMILES of the shared MMP scaffold core (24 unique values) |
| `n_shared_cores` | int | Number of MMP pairs sharing this core across the dataset |
| `delta_pic50` | float | pIC50 difference for this pair |
| `tanimoto` | float | ECFP4 Tanimoto for this pair |
| `cliff_tier` | str | Cliff tier of this pair (severe/moderate/broad/non_cliff) |
| `mmp_large_delta` | bool | True if ΔpIC50 ≥ 2.0 for this MMP pair |
| `hub_a` | bool | True if compound A is a cliff hub (Class A or B) |
| `hub_b` | bool | True if compound B is a cliff hub (Class A or B) |
| `mmp_type` | str | Change type: `single_atom_change`, `small_substituent`, `medium_substituent`, `large_substituent` |
| `is_canonical_severe_cliff` | bool | True if this MMP pair is also in the Tanimoto-defined severe cliff set (85 pairs, 90.4%) |

---

## mmp_discontinuity_scores.csv

Per-compound MMP discontinuity scores for the 99 severe cliff compounds. 99 rows × 5 columns.

| Column | Type | Description |
|--------|------|-------------|
| `inchi_key` | str | InChIKey |
| `pic50` | float | Consensus pIC50 |
| `hub_class` | str | `A`, `B`, or `none` |
| `mmp_discontinuity_score` | float | Mean ΔpIC50 across all MMP partners; measures structural attractor strength |
| `n_mmp_partners` | int | Number of MMP partners among the cliff compound set |

---

## fingerprint_sensitivity_94pairs.csv / tables/table_s6_fingerprint_sensitivity.csv

ECFP4 vs ECFP6 fingerprint sensitivity analysis for all 94 severe cliff pairs.
94 rows. Full Supplementary Table S6.

| Column | Type | Description |
|--------|------|-------------|
| `inchi_key_a` | str | InChIKey of compound A |
| `inchi_key_b` | str | InChIKey of compound B |
| `tanimoto_ecfp4` | float | ECFP4 Tanimoto (radius=2, 2048 bits) |
| `tanimoto_ecfp6` | float | ECFP6 Tanimoto (radius=3, 2048 bits) |
| `delta_pic50` | float | pIC50 difference |
| `mmp_confirmed` | bool | True if confirmed as MMP pair |
| `hub_involved` | bool | True if either compound is a cliff hub |
| `ecfp4_only_cliff` | bool | True if not confirmed by ECFP6 ≥0.8 or MMP |
| `robust_at_ecfp6` | bool | True if Tanimoto ≥ 0.8 under both ECFP4 and ECFP6 |

**Summary:** 30 pairs robust at ECFP6; 64 pairs non-robust (ECFP4 Tan 0.80–0.85, ECFP6 <0.80);
13 pairs ECFP4-only (not confirmed by ECFP6 or MMP).
Hub dominance statistic is fingerprint-invariant: 53.2% under ECFP4, 53.3% under ECFP6.
