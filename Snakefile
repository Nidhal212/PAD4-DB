# PAD4-DB v2 — Snakemake Workflow
# ============================================================
# Run from project root with conda env management:
#
#   snakemake --cores 4 --use-conda
#
# Dry-run (show DAG without executing):
#   snakemake --dry-run --use-conda
#
# DAG visualisation:
#   snakemake --dag | dot -Tpdf > dag.pdf
# ============================================================

configfile: "config/config.yaml"

from pathlib import Path

# ── Directory setup ──────────────────────────────────────────────────────────

DIRS = [
    "data/interim/standardized",
    "data/interim/normalized",
    "data/processed",
    "outputs/tables",
    "outputs/audit",
    "logs",
]

onstart:
    for d in DIRS:
        Path(d).mkdir(parents=True, exist_ok=True)


# ── Target rule ──────────────────────────────────────────────────────────────

rule all:
    input:
        # Core outputs
        "data/processed/pad4_compounds.parquet",
        "data/processed/hts_compound_index.parquet",
        "logs/independence_scores.done",
        # SAR analysis
        "data/processed/activity_cliffs.parquet",
        "data/processed/activity_pairs_sim_ge06.parquet",
        "outputs/tables/05_scaffold_summary.csv",
        "outputs/tables/05_cliff_summary.json",
        "outputs/tables/05_patent_exclusive_cliff_contribution.json",
        # Audits
        "outputs/audit/A1_audit_report.txt",
        "outputs/audit/A2_reference_recovery.json",


# ── Rule 00: Raw inventory QC ────────────────────────────────────────────────
# Validates that all expected raw files are present and non-empty.
# Produces a sentinel so downstream rules don't re-run unnecessarily.

rule qc_inventory:
    input:
        ancient("data/raw")
    output:
        touch("logs/qc_inventory.done")
    conda:
        "workflow/envs/pad4bench.yaml"
    log:
        "logs/00_qc_inventory.log"
    shell:
        "python 00_check_raw_inventory.py > {log} 2>&1"


# ── Rule 01: SMILES standardization ─────────────────────────────────────────
# Reads all raw PubChem CSVs, ChEMBL CSV, and BindingDB TSV.
# Strips salts, chooses largest fragment, assigns InChIKey.

rule standardize:
    input:
        sentinel = "logs/qc_inventory.done"
    output:
        parquet = "data/interim/standardized/standardized_compounds.parquet",
        report  = "outputs/tables/01_standardization_report.csv",
    conda:
        "workflow/envs/pad4bench.yaml"
    log:
        "logs/01_standardize.log"
    script:
        "scripts/01_standardize/01_standardize_smiles.py"


# ── Rule 02: Activity normalization ─────────────────────────────────────────
# Converts all raw bioactivity values to pIC50 (nM scale).
# Percent-inhibition HTS rows are intercepted before unit conversion.

rule normalize:
    input:
        standardized = "data/interim/standardized/standardized_compounds.parquet",
    output:
        parquet = "data/interim/normalized/normalized_activities.parquet",
        report  = "outputs/tables/02_normalization_report.csv",
    conda:
        "workflow/envs/pad4bench.yaml"
    log:
        "logs/02_normalize.log"
    script:
        "scripts/02_normalize/02_normalize_activities.py"


# ── Rule 03: Replicate aggregation ──────────────────────────────────────────
# Groups by InChIKey × source × aid × endpoint_type; aggregates replicates
# by log-mean; assigns assay_mechanism_class.

rule aggregate:
    input:
        normalized   = "data/interim/normalized/normalized_activities.parquet",
        standardized = "data/interim/standardized/standardized_compounds.parquet",
    output:
        parquet = "data/interim/normalized/replicate_aggregated.parquet",
    conda:
        "workflow/envs/pad4bench.yaml"
    log:
        "logs/03_aggregate.log"
    script:
        "scripts/03_aggregate/03_replicate_aggregate.py"


# ── Rule 03a: Potency / HTS space split ─────────────────────────────────────
# Splits replicate_aggregated.parquet into two non-overlapping subsets:
#   potency_space  (use_in_potency_model == True)
#   hts_space      (use_in_potency_model == False)

rule split_spaces:
    input:
        aggregated = "data/interim/normalized/replicate_aggregated.parquet",
    output:
        potency = "data/interim/normalized/potency_space.parquet",
        hts     = "data/interim/normalized/hts_space.parquet",
    conda:
        "workflow/envs/pad4bench.yaml"
    log:
        "logs/03a_split.log"
    script:
        "scripts/03_aggregate/03a_split_spaces.py"


# ── Rule 04: Cross-AID dedup + compound assembly ────────────────────────────
# Deduplicates by InChIKey × source × endpoint_type.
# Applies AID 2202576/77 preference rule (RFMS overlap).
# Assembles pad4_compounds.parquet (one row per unique compound).

rule dedup_assemble:
    input:
        potency      = "data/interim/normalized/potency_space.parquet",
        hts          = "data/interim/normalized/hts_space.parquet",
        standardized = "data/interim/standardized/standardized_compounds.parquet",
        normalized   = "data/interim/normalized/normalized_activities.parquet",
    output:
        compounds = "data/processed/pad4_compounds.parquet",
        hts_index = "data/processed/hts_compound_index.parquet",
        aid_level = "data/interim/normalized/dedup_aid_level.parquet",
    conda:
        "workflow/envs/pad4bench.yaml"
    log:
        "logs/04_dedup.log"
    script:
        "scripts/04_dedup/04_dedup_and_assemble.py"


# ── Rule 04b: Source independence scores ────────────────────────────────────
# Adds source_independence_score and is_true_multi_source columns.
# Updates pad4_compounds.parquet in-place.
# Uses a sentinel output to avoid Snakemake input/output path collision.

rule independence_scores:
    input:
        compounds = "data/processed/pad4_compounds.parquet",
    output:
        touch("logs/independence_scores.done"),
    conda:
        "workflow/envs/pad4bench.yaml"
    log:
        "logs/04b_independence.log"
    script:
        "scripts/04_dedup/04b_add_independence_scores.py"


# ── Rule 05: Scaffold + activity-cliff analysis ──────────────────────────────
# Computes ECFP4 fingerprints, Murcko scaffolds, all-pairs Tanimoto,
# cliff classification, SALI, MMP, and hub detection.

rule sar_analysis:
    input:
        compounds = "data/processed/pad4_compounds.parquet",
        sentinel  = "logs/independence_scores.done",
    output:
        cliffs      = "data/processed/activity_cliffs.parquet",
        pairs       = "data/processed/activity_pairs_sim_ge06.parquet",
        scaffold    = "outputs/tables/05_scaffold_summary.csv",
        cliff_json  = "outputs/tables/05_cliff_summary.json",
        patent_json = "outputs/tables/05_patent_exclusive_cliff_contribution.json",
    conda:
        "workflow/envs/pad4bench.yaml"
    log:
        "logs/05_sar.log"
    script:
        "scripts/05_cliffs/05_scaffold_and_cliffs.py"


# ── Rule A1: Biological & chemical audit ─────────────────────────────────────
# Verifies target identity, species, assay mechanism, and chemical correctness.
# VERIFICATION ONLY — does not modify any data files.

rule audit_biological:
    input:
        compounds  = "data/processed/pad4_compounds.parquet",
        aid_level  = "data/interim/normalized/dedup_aid_level.parquet",
        normalized = "data/interim/normalized/normalized_activities.parquet",
        reagg      = "data/interim/normalized/replicate_aggregated.parquet",
        sentinel   = "logs/independence_scores.done",
    output:
        report     = "outputs/audit/A1_audit_report.txt",
        aid_csv    = "outputs/audit/A1_aid_audit.csv",
        bio_csv    = "outputs/audit/A1_compound_bio_audit.csv",
        chem_csv   = "outputs/audit/A1_chemical_correctness.csv",
        summary    = "outputs/audit/A1_audit_summary.json",
    conda:
        "workflow/envs/pad4bench.yaml"
    log:
        "logs/A1_audit.log"
    script:
        "scripts/audit/A1_biological_audit.py"


# ── Rule A2: Reference compound recovery ─────────────────────────────────────
# Traces 14 curated PAD4 inhibitors through the full pipeline.
# VERIFICATION ONLY — does not modify any data files.

rule audit_reference:
    input:
        compounds = "data/processed/pad4_compounds.parquet",
        sentinel  = "logs/independence_scores.done",
    output:
        report    = "outputs/audit/A2_reference_recovery.json",
    conda:
        "workflow/envs/pad4bench.yaml"
    log:
        "logs/A2_audit.log"
    script:
        "scripts/audit/A2_reference_compound_recovery.py"
