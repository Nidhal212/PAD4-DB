# ARCHIVE INDEX — PAD4-DB v1.0
**Date archived:** 2026-06-22  
**Policy:** Files are preserved, not deleted. This index documents what was archived and why.

---

## Superseded Manuscript Drafts

| Original path | Archive location | Reason |
|--------------|-----------------|--------|
| `publication/manuscript/PAD4_DB_manuscript_DRAFT_v7.md` | Source preserved in place | Superseded by FINAL.md — same content plus Discussion revision and M5/M18 hardening additions |
| `publication/manuscript/PAD4_DB_manuscript_DRAFT_v7.docx` | Source preserved in place | Derived from DRAFT_v7.md — superseded |
| `publication/manuscript/PAD4_DB_manuscript_DRAFT_v7.pdf` | Source preserved in place | Derived from DRAFT_v7.md — superseded |
| `publication/manuscript/PAD4_DB_v2_manuscript_integrated.docx` | Source preserved in place | Pre-v7 integrated version — older manuscript with different structure |
| `publication/manuscript/PAD4_DB_v2_manuscript_integrated.pdf` | Source preserved in place | Derived from integrated — superseded |

---

## Obsolete Manuscript Builders

| Original path | Archive location | Reason |
|--------------|-----------------|--------|
| `publication/scripts/build_manuscript_v7_docx.py` | `archive/obsolete_builders/build_manuscript_v7_docx.py` | Hand-duplicated text; drift risk; superseded by build_manuscript_from_md.py |
| `publication/scripts/build_manuscript_docx.py` | `archive/obsolete_builders/build_manuscript_docx.py` | Older builder for superseded manuscript version |
| `scripts/build_manuscript_docx.py` (root) | Root legacy scripts | Pre-nature-style builder; superseded |

Note: `archive/obsolete_builders/` already existed with these files. Archive action confirmed and indexed here.

---

## Superseded Supplementary Figures (Pre-v7)

Location: `publication/figures/supplementary/_archive_pre_v7/` (already archived)

| File | Reason |
|------|--------|
| `fig_s01_pipeline.png/.pdf` | Pipeline flowchart — removed in v7 restructure; not in FINAL manuscript |
| `fig_s02_sali.png/.pdf` | SALI visualization — merged into SAS map (S3) |
| `fig_s03_patent.png/.pdf` | Patent scaffold analysis — moved to main text Fig 4 |
| `fig_s04_independence.png/.pdf` | Independence scores — superseded |
| `fig_s04_reference_recovery.png/.pdf` | Reference recovery — from audit (A2), not in FINAL |
| `fig_s05_scaffold_structures.png/.pdf` | Scaffold structures — removed in v7 |
| `fig_s06_permutation.png/.pdf` | Old permutation figure — superseded by S5 (null models) |
| `fig_s07_physicochemical.png/.pdf` | Physicochemical properties — removed from FINAL |

---

## Legacy Root Scripts (`scripts/` directory)

These scripts pre-date the `publication/scripts/` reorganization and have been superseded. They remain in place but should not be confused with the canonical versions.

| Path | Reason |
|------|--------|
| `scripts/fig1_pipeline_workflow.py` through `fig11_*.py` | Early-development figure scripts, pre-nature-style palette |
| `scripts/nature/fig01_*.py` through `fig11_*.py` | Nature-style versions that were superseded by `publication/scripts/figures/` |
| `scripts/nature/generate_tables.py`, `tables_nature.py` | Superseded by `publication/scripts/figures/generate_tables.py` |
| `scripts/audit/E3_permutation_test.py`, `E3_render_figure_only.py` | Superseded by `publication/scripts/analysis/audit_constrained_permutation.py` |
| `scripts/diagnostic_mechanism.py` | Debugging script — development artifact |
| `scripts/repro/run_repro_pipeline.py` | Reproducibility runner — internal audit tool |
| `scripts/stress_test/run_stress_test.py` | Stress test runner — internal audit tool |
| `scripts/0000.py` | Empty scratch script |
| `scripts/final_audit.py`, `generate_pdf.py`, `generate_project_state_report.py` | Utility scripts from earlier development |
| `scripts/05_cliffs/05e_golden_set.py` | Superseded by publication version |

---

## Dryrun Data Files

All `*_dryrun.parquet` files throughout `data/` are development-phase test runs. They are not needed for the final release and should not be confused with canonical outputs.

| Path pattern | Count | Reason |
|-------------|-------|--------|
| `data/processed/*_dryrun.parquet` | 4 files | Pipeline development test runs |
| `data/interim/normalized/*_dryrun.parquet` | 3 files | Pipeline development test runs |
| `data/interim/standardized/*_dryrun.parquet` | 1 file | Pipeline development test runs |
| `outputs/tables/*_dryrun.*` | 3 files | Table QC development runs |

---

## Dated Backup Snapshots

| File | Reason |
|------|--------|
| `data/processed/pad4_compounds_pre_columns_2026-06-15.parquet` | Pre-column-addition backup; superseded by current pad4_compounds.parquet |
| `data/processed/pad4_compounds_pre_remediation_2026-06-16.parquet` | Pre-remediation snapshot; superseded |

---

## Internal Process Documents (PRIVATE — NOT FOR PUBLIC RELEASE)

These files are retained locally but must never appear in a public repository or Zenodo deposition.

| File | Reason |
|------|--------|
| `DEPLOYMENT_REPORT.md` | Internal deployment process notes |
| `FINAL_REVIEWER_RESPONSE_MATRIX.md` | Reviewer simulation — pre-submission strategy |
| `FIGURE_READINESS_REMEDIATION_REPORT.md` | Internal figure audit |
| `HANDOFF.md` | Project handoff with unpublished interpretation |
| `REPOSITORY_AND_DEPOSITION_AUDIT.md` | Internal deposition strategy |
| `FREEZE_REPORT.md` | Pre-freeze audit |
| `FREEZE_DECISION.md` | Pre-freeze decision notes |
| `builder_audit_report.md` | Internal builder comparison |
| `git_safety_report.md` | Internal git safety analysis |
| `smoothness_harmonization_report.md` | Internal editorial decision |
| `SUBMISSION_READINESS_REPORT.md` | Internal readiness notes |
| `deploy_commands.sh` | Internal deployment commands |
| `repro_audit/REPRODUCTION_AUDIT_REPORT.md` | Internal reproducibility audit |
| `outputs/audit/PEER_REVIEW_SIMULATION.md` | Reviewer simulation — never public |
| `outputs/audit/PAD4_DB_ANALYSIS_EXPANSION_REPORT.md` | Internal expansion planning |
| `outputs/audit/MANUSCRIPT_REVISION_CHANGELOG.md` | Internal revision history |

---

## Scratch / Utility Files

| File | Reason |
|------|--------|
| `covaltest.py` | Covalent chemistry debugging script |
| `AID_dowloed.py` | PubChem download utility (note: typo in filename) |
| `scripts/0000.py` | Empty scratch file |
| `E4_cliff_pairs_v2.json` | Development version superseded by E4_cliff_pairs.json |
| `outputs/audit/E5_golden_set.json` | Unused golden set — golden set concept abandoned (audit A2 used instead) |

---

## Repro Audit Copies

| Directory | Reason |
|-----------|--------|
| `data/interim_repro/` | Reproducibility audit pipeline re-run. Correct, but redundant with `data/interim/` |
| `data/processed_repro/` | Reproducibility audit pipeline re-run. Verified identical to `data/processed/` |
| `repro_audit/*.log` | Run logs from reproducibility audit. Internal reference |
