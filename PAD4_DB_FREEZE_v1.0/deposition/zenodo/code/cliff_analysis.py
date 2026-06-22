#!/usr/bin/env python3
"""
PAD4-DB v2 — Step 05: Scaffold Analysis & Activity Cliff Detection
scripts/05_cliffs/05_scaffold_and_cliffs.py

Inputs:  data/processed/pad4_compounds.parquet
Outputs: data/processed/activity_cliffs.parquet
         data/processed/activity_pairs_sim_ge06.parquet
         outputs/tables/05_scaffold_summary.csv
         outputs/tables/05_cliff_summary.json
         outputs/tables/05_patent_exclusive_cliff_contribution.json
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

# ── Paths ──────────────────────────────────────────────────────────────────
COMP_PARQ    = Path("data/processed/pad4_compounds.parquet")
OUT_CLIFFS   = Path("data/processed/activity_cliffs.parquet")
OUT_PAIRS    = Path("data/processed/activity_pairs_sim_ge06.parquet")
OUT_SCAFFOLD = Path("outputs/tables/05_scaffold_summary.csv")
OUT_CLIFF_JSON  = Path("outputs/tables/05_cliff_summary.json")
OUT_PATENT_JSON = Path("outputs/tables/05_patent_exclusive_cliff_contribution.json")

BATCH_SIZE    = 500
SIM_THRESHOLD = 0.6


# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(description="Step 05: scaffold analysis & cliff detection")
    p.add_argument("--dry-run", action="store_true",
                   help="Use first 300 compounds; suffix outputs with _dryrun")
    p.add_argument("--no-patent-analysis", action="store_true",
                   help="Skip Section 6 patent-exclusive contribution analysis")
    return p.parse_args()


# =============================================================================
def load_compounds(args) -> pd.DataFrame:
    if not COMP_PARQ.exists():
        sys.exit(f"ERROR: {COMP_PARQ} not found")

    cols = ["inchi_key", "smiles_std", "pic50_consensus", "source_list",
            "n_sources", "n_heavy_atoms", "mol_weight", "use_in_potency_model"]
    df = pd.read_parquet(COMP_PARQ, columns=cols)

    before = len(df)
    df = df[
        (df["use_in_potency_model"] == True) &
        df["smiles_std"].notna() &
        df["pic50_consensus"].notna()
    ].copy().reset_index(drop=True)
    print(f"Loaded {before:,} rows → {len(df):,} pass filter "
          f"(use_in_potency_model + smiles_std + pic50_consensus)")

    if args.dry_run:
        df = df.head(300).copy().reset_index(drop=True)
        print(f"DRY RUN: truncated to {len(df)} compounds")

    return df


# =============================================================================
def compute_fingerprints(df: pd.DataFrame) -> tuple[np.ndarray, list, list]:
    """ECFP4 (Morgan r=2, nBits=2048). Returns (fp_matrix, valid_row_indices, failure_iks)."""
    fp_rows:   list = []
    valid_idx: list = []
    failures:  list = []

    for i, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["smiles_std"])
        if mol is None:
            print(f"  WARNING: SMILES parse failed — {row['inchi_key']}")
            failures.append(row["inchi_key"])
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        fp_rows.append(np.array(fp, dtype=np.uint8))
        valid_idx.append(i)

    if fp_rows:
        fp_matrix = np.array(fp_rows, dtype=np.uint8)
    else:
        fp_matrix = np.zeros((0, 2048), dtype=np.uint8)

    print(f"Fingerprints: {len(fp_rows):,} OK, {len(failures)} failures")
    return fp_matrix, valid_idx, failures


# =============================================================================
def compute_pairwise_similarity(
        fp_matrix: np.ndarray, ik_list: list
) -> tuple[list, list, list]:
    """
    Vectorized batch Tanimoto for all (i < j) pairs.
    Returns (ika_list, ikb_list, tan_list) for pairs with Tanimoto >= 0.6.
    """
    n = len(fp_matrix)
    fp_f32     = fp_matrix.astype(np.float32)
    bit_counts = fp_f32.sum(axis=1)   # (n,) — popcount per compound

    ika_list: list = []
    ikb_list: list = []
    tan_list: list = []

    n_batches = (n + BATCH_SIZE - 1) // BATCH_SIZE

    for b, i0 in enumerate(range(0, n, BATCH_SIZE)):
        i1    = min(i0 + BATCH_SIZE, n)
        batch = fp_f32[i0:i1]              # (bs, 2048)

        dot   = batch @ fp_f32.T           # (bs, n)
        bc_b  = bit_counts[i0:i1].reshape(-1, 1)
        union = bc_b + bit_counts.reshape(1, -1) - dot  # (bs, n)

        with np.errstate(divide="ignore", invalid="ignore"):
            tan = np.where(union > 0, dot / union, 0.0)  # (bs, n)

        # Upper-triangle mask: for local row k, valid j > i0+k
        row_idx   = np.arange(i1 - i0).reshape(-1, 1)  # (bs, 1)
        col_idx   = np.arange(n).reshape(1, -1)         # (1, n)
        upper     = col_idx > (i0 + row_idx)            # (bs, n) bool

        hit_r, hit_c = np.where((tan >= SIM_THRESHOLD) & upper)

        for r, c in zip(hit_r.tolist(), hit_c.tolist()):
            ika_list.append(ik_list[i0 + r])
            ikb_list.append(ik_list[c])
            tan_list.append(float(tan[r, c]))

        print(f"  Batch {b+1}/{n_batches} (rows {i0}–{i1-1}): "
              f"{len(ika_list):,} pairs found so far")

    print(f"Total pairs with Tanimoto >= {SIM_THRESHOLD}: {len(ika_list):,}")
    return ika_list, ikb_list, tan_list


# =============================================================================
def build_pair_df(ika: list, ikb: list, tan: list, df: pd.DataFrame) -> pd.DataFrame:
    """Construct pair DataFrame from flat similarity lists."""
    pic50_map  = df.set_index("inchi_key")["pic50_consensus"].to_dict()
    source_map = df.set_index("inchi_key")["source_list"].to_dict()

    p_a = [pic50_map[k] for k in ika]
    p_b = [pic50_map[k] for k in ikb]
    s_a = [source_map[k] for k in ika]
    s_b = [source_map[k] for k in ikb]

    pair_df = pd.DataFrame({
        "inchi_key_a":        ika,
        "inchi_key_b":        ikb,
        "tanimoto":           np.array(tan, dtype=np.float32),
        "pic50_a":            p_a,
        "pic50_b":            p_b,
        "delta_pic50":        [abs(a - b) for a, b in zip(p_a, p_b)],
        "source_a":           s_a,
        "source_b":           s_b,
        "source_combination": ["|".join(sorted([a, b])) for a, b in zip(s_a, s_b)],
        "same_source":        [a == b for a, b in zip(s_a, s_b)],
    })
    print(f"Pair DataFrame: {len(pair_df):,} rows")
    return pair_df


# =============================================================================
def classify_cliffs(pair_df: pd.DataFrame) -> pd.DataFrame:
    """Vectorized cliff tier assignment + patent-exclusive flags."""
    t = pair_df["tanimoto"]
    d = pair_df["delta_pic50"]

    pair_df = pair_df.copy()
    pair_df["cliff_tier"] = np.select(
        [(t >= 0.8) & (d >= 2.0),
         (t >= 0.8) & (d >= 1.5),
         (t >= 0.8) & (d >= 1.0)],
        ["severe", "moderate", "broad"],
        default="non_cliff",
    )
    pair_df["patent_exclusive_a"] = pair_df["source_a"] == "pubchem_confirmatory"
    pair_df["patent_exclusive_b"] = pair_df["source_b"] == "pubchem_confirmatory"
    pair_df["any_patent_exclusive"] = (
        pair_df["patent_exclusive_a"] | pair_df["patent_exclusive_b"]
    )
    return pair_df


# =============================================================================
def compute_scaffold_analysis(
        df: pd.DataFrame,
) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    """
    Bemis-Murcko decomposition.
    Returns (scaffold_summary_df, stats_dict, df_with_scaffold_column).
    """
    scaffolds: list = []
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["smiles_std"])
        if mol is None:
            scaffolds.append("NO_SCAFFOLD")
            continue
        try:
            scaf = MurckoScaffold.GetScaffoldForMol(mol)
            smi  = Chem.MolToSmiles(scaf) if scaf is not None else ""
            scaffolds.append(smi if smi else "NO_SCAFFOLD")
        except Exception:
            scaffolds.append("NO_SCAFFOLD")

    df = df.copy()
    df["scaffold_smiles"] = scaffolds

    patent_iks = set(df[df["source_list"] == "pubchem_confirmatory"]["inchi_key"])

    grp = df.groupby("scaffold_smiles").agg(
        n_compounds=("inchi_key",        "count"),
        mean_pic50 =("pic50_consensus",   "mean"),
        std_pic50  =("pic50_consensus",   "std"),
        pic50_min  =("pic50_consensus",   "min"),
        pic50_max  =("pic50_consensus",   "max"),
    ).reset_index()
    grp["pic50_range"] = grp["pic50_max"] - grp["pic50_min"]
    grp["std_pic50"]   = grp["std_pic50"].fillna(0.0)

    patent_scafs = set(df[df["inchi_key"].isin(patent_iks)]["scaffold_smiles"])
    grp["contains_patent_exclusive"] = grp["scaffold_smiles"].isin(patent_scafs)
    grp = grp.sort_values("n_compounds", ascending=False).reset_index(drop=True)

    n_unique = len(grp)
    n_single = int((grp["n_compounds"] == 1).sum())
    n_series = int((grp["n_compounds"] >= 2).sum())
    largest  = int(grp["n_compounds"].max())
    top_n    = int(grp["n_compounds"].iloc[0])

    print(f"Scaffolds: {n_unique:,} unique | {n_single:,} singletons | "
          f"{n_series:,} series (≥2) | largest={largest}")

    stats = {
        "n_unique_scaffolds":       n_unique,
        "n_singleton_scaffolds":    n_single,
        "n_series_scaffolds":       n_series,
        "largest_series_size":      largest,
        "scaffold_coverage":        round(n_series / n_unique, 4) if n_unique else 0.0,
        "top_scaffold_n_compounds": top_n,
    }
    return grp, stats, df


# =============================================================================
def _cliff_counts(pair_df: pd.DataFrame, n_total: int) -> dict:
    """Aggregate cliff statistics from a classified pair DataFrame."""
    severe_df = pair_df[pair_df["cliff_tier"] == "severe"]
    cliff_df  = pair_df[pair_df["cliff_tier"] != "non_cliff"]

    cliff_cpds  = set(cliff_df["inchi_key_a"])  | set(cliff_df["inchi_key_b"])
    severe_cpds = set(severe_df["inchi_key_a"]) | set(severe_df["inchi_key_b"])

    max_d     = float(pair_df["delta_pic50"].max())    if len(pair_df)   else 0.0
    mean_d_sv = float(severe_df["delta_pic50"].mean()) if len(severe_df) else 0.0

    return {
        "n_severe_cliffs_sim08_delta20":   int(len(severe_df)),
        "n_moderate_cliffs_sim08_delta15": int((pair_df["cliff_tier"] == "moderate").sum()),
        "n_broad_cliffs_sim08_delta10":    int((pair_df["cliff_tier"] == "broad").sum()),
        "n_compounds_in_severe_cliffs":    len(severe_cpds),
        "n_compounds_in_any_cliff":        len(cliff_cpds),
        "pct_compounds_in_severe_cliffs":  round(len(severe_cpds) / n_total * 100, 2)
                                           if n_total else 0.0,
        "max_delta_pic50":                 round(max_d, 4),
        "mean_delta_pic50_severe":         round(mean_d_sv, 4),
    }


# =============================================================================
def compute_patent_contribution(
        full_pair_df: pd.DataFrame, df: pd.DataFrame,
) -> dict:
    """
    Re-run pairwise similarity on dataset B (patent-exclusive compounds removed).
    Returns patent_contribution_json dict.
    """
    patent_mask = df["source_list"] == "pubchem_confirmatory"
    n_patent    = int(patent_mask.sum())
    df_b        = df[~patent_mask].copy().reset_index(drop=True)
    n_b         = len(df_b)

    print(f"\nDataset B: {n_b:,} compounds (excl. {n_patent} patent-exclusive)")

    fp_b, valid_b, _ = compute_fingerprints(df_b)
    ik_b      = df_b.iloc[valid_b]["inchi_key"].tolist()
    df_b_val  = df_b.iloc[valid_b].copy().reset_index(drop=True)

    print("  Pairwise similarity for dataset B...")
    ika_b, ikb_b, tan_b = compute_pairwise_similarity(fp_b, ik_b)

    if ika_b:
        pair_b = build_pair_df(ika_b, ikb_b, tan_b, df_b_val)
        pair_b = classify_cliffs(pair_b)
    else:
        pair_b = pd.DataFrame(columns=full_pair_df.columns)

    counts_a = _cliff_counts(full_pair_df, len(df))
    counts_b = _cliff_counts(pair_b, n_b)

    severe_a = counts_a["n_severe_cliffs_sim08_delta20"]
    severe_b = counts_b["n_severe_cliffs_sim08_delta20"]
    delta    = severe_a - severe_b
    pct      = round(delta / severe_a * 100, 2) if severe_a > 0 else 0.0

    verdict = ("strong_contribution"   if pct >= 20 else
               "moderate_contribution" if pct >= 5  else
               "weak_contribution")

    return {
        "dataset_a_n_compounds": len(df),
        "dataset_b_n_compounds": n_b,
        "severe_cliffs_a":       severe_a,
        "severe_cliffs_b":       severe_b,
        "cliff_delta":           delta,
        "cliff_delta_pct":       pct,
        "verdict":               verdict,
    }


# =============================================================================
def write_outputs(
        pair_df:   pd.DataFrame,
        cliff_df:  pd.DataFrame,
        scaf_df:   pd.DataFrame,
        cliff_json:  dict,
        patent_json: dict,
        args,
) -> None:
    sfx = "_dryrun" if args.dry_run else ""

    def _p(base: Path) -> Path:
        return base.parent / (base.stem + sfx + base.suffix)

    pair_cols = [
        "inchi_key_a", "inchi_key_b", "tanimoto", "pic50_a", "pic50_b",
        "delta_pic50", "source_a", "source_b", "source_combination",
        "same_source", "cliff_tier", "patent_exclusive_a",
        "patent_exclusive_b", "any_patent_exclusive",
    ]
    cliff_cols = [
        "inchi_key_a", "inchi_key_b", "tanimoto", "pic50_a", "pic50_b",
        "delta_pic50", "cliff_tier", "source_combination",
        "patent_exclusive_a", "patent_exclusive_b", "any_patent_exclusive",
    ]
    scaf_cols = [
        "scaffold_smiles", "n_compounds", "mean_pic50", "std_pic50",
        "pic50_range", "contains_patent_exclusive",
    ]

    pair_df[pair_cols].to_parquet(_p(OUT_PAIRS), index=False)
    print(f"Written → {_p(OUT_PAIRS)}")

    cliff_df[cliff_cols].to_parquet(_p(OUT_CLIFFS), index=False)
    print(f"Written → {_p(OUT_CLIFFS)}")

    scaf_df[scaf_cols].to_csv(_p(OUT_SCAFFOLD), index=False)
    print(f"Written → {_p(OUT_SCAFFOLD)}")

    _p(OUT_CLIFF_JSON).write_text(json.dumps(cliff_json, indent=2))
    print(f"Written → {_p(OUT_CLIFF_JSON)}")

    _p(OUT_PATENT_JSON).write_text(json.dumps(patent_json, indent=2))
    print(f"Written → {_p(OUT_PATENT_JSON)}")


# =============================================================================
def print_summary(cliff_json: dict, patent_json: dict) -> None:
    print("\n" + "=" * 60)
    print("STEP 05 COMPLETE")
    print("=" * 60)
    print(f"Compounds analyzed: {cliff_json['n_compounds_analyzed']:,}")
    print(f"FP failures:        {cliff_json['n_fingerprint_failures']}")
    print("\nSimilarity landscape:")
    for k, v in cliff_json["similarity_landscape"].items():
        print(f"  {k}: {v:,}")
    print("\nCliff counts:")
    for k, v in cliff_json["cliff_counts"].items():
        print(f"  {k}: {v}")
    print("\nScaffold analysis:")
    for k, v in cliff_json["scaffold_analysis"].items():
        print(f"  {k}: {v}")
    print("\nPatent exclusive contribution:")
    for k, v in cliff_json["patent_exclusive_contribution"].items():
        print(f"  {k}: {v}")
    if patent_json:
        print("\nPatent contribution JSON:")
        for k, v in patent_json.items():
            print(f"  {k}: {v}")


# =============================================================================
def main():
    args = parse_args()

    Path("data/processed").mkdir(parents=True, exist_ok=True)
    Path("outputs/tables").mkdir(parents=True, exist_ok=True)

    # ── Load ─────────────────────────────────────────────────────────────────
    df       = load_compounds(args)
    n_analy  = len(df)

    # ── Fingerprints ──────────────────────────────────────────────────────────
    print(f"\nComputing ECFP4 fingerprints for {n_analy:,} compounds...")
    fp_matrix, valid_idx, failures = compute_fingerprints(df)
    df_valid  = df.iloc[valid_idx].copy().reset_index(drop=True)
    ik_list   = df_valid["inchi_key"].tolist()
    n_fps     = len(fp_matrix)

    # ── Pairwise similarity ───────────────────────────────────────────────────
    n_pairs_total = n_fps * (n_fps - 1) // 2
    print(f"\nComputing pairwise Tanimoto: {n_fps:,} compounds → ~{n_pairs_total:,} pairs...")
    ika, ikb, tan = compute_pairwise_similarity(fp_matrix, ik_list)

    # ── Pair DataFrame + cliff classification ─────────────────────────────────
    if ika:
        pair_df = build_pair_df(ika, ikb, tan, df_valid)
    else:
        pair_df = pd.DataFrame(columns=[
            "inchi_key_a","inchi_key_b","tanimoto","pic50_a","pic50_b",
            "delta_pic50","source_a","source_b","source_combination","same_source",
        ])
    pair_df  = classify_cliffs(pair_df)
    cliff_df = pair_df[pair_df["cliff_tier"] != "non_cliff"].copy()

    # ── Similarity landscape ──────────────────────────────────────────────────
    sim_counts = {
        "n_pairs_sim_ge_06":  int((pair_df["tanimoto"] >= 0.60).sum()),
        "n_pairs_sim_ge_07":  int((pair_df["tanimoto"] >= 0.70).sum()),
        "n_pairs_sim_ge_08":  int((pair_df["tanimoto"] >= 0.80).sum()),
        "n_pairs_sim_ge_09":  int((pair_df["tanimoto"] >= 0.90).sum()),
        "n_pairs_sim_ge_095": int((pair_df["tanimoto"] >= 0.95).sum()),
    }

    cliff_counts = _cliff_counts(pair_df, n_fps)

    # ── Scaffold analysis ──────────────────────────────────────────────────────
    print(f"\nComputing Bemis-Murcko scaffolds for {n_fps:,} compounds...")
    scaf_df, scaffold_stats, df_with_scaffold = compute_scaffold_analysis(df_valid)

    # ── Patent contribution (Section 6) ───────────────────────────────────────
    patent_json: dict  = {}
    n_patent_cpds = int((df_valid["source_list"] == "pubchem_confirmatory").sum())

    patent_stats = {
        "n_patent_exclusive_compounds":                  n_patent_cpds,
        "severe_cliffs_all_compounds":                   cliff_counts["n_severe_cliffs_sim08_delta20"],
        "severe_cliffs_without_patent_exclusive":        None,
        "severe_cliffs_added_by_patent":                 None,
        "pct_severe_cliffs_involving_patent_exclusive":  None,
        "new_scaffolds_from_patent_exclusive":           None,
    }

    if not args.no_patent_analysis:
        patent_json = compute_patent_contribution(pair_df, df_with_scaffold)

        # % severe cliffs where either compound is patent-exclusive
        n_severe = cliff_counts["n_severe_cliffs_sim08_delta20"]
        n_severe_involving = int(
            cliff_df[cliff_df["cliff_tier"] == "severe"]["any_patent_exclusive"].sum()
        )
        pct_involving = round(n_severe_involving / n_severe * 100, 2) if n_severe > 0 else 0.0

        # Scaffolds that exist only among patent-exclusive compounds
        patent_iks = set(
            df_with_scaffold[df_with_scaffold["source_list"] == "pubchem_confirmatory"]["inchi_key"]
        )
        non_patent_scafs = set(
            df_with_scaffold[~df_with_scaffold["inchi_key"].isin(patent_iks)]["scaffold_smiles"]
        )
        patent_only_scafs = (
            set(df_with_scaffold[df_with_scaffold["inchi_key"].isin(patent_iks)]["scaffold_smiles"])
            - non_patent_scafs
        )
        patent_only_scafs.discard("NO_SCAFFOLD")
        n_new_scafs = len(patent_only_scafs)

        patent_stats.update({
            "severe_cliffs_without_patent_exclusive":       patent_json["severe_cliffs_b"],
            "severe_cliffs_added_by_patent":                patent_json["cliff_delta"],
            "pct_severe_cliffs_involving_patent_exclusive": pct_involving,
            "new_scaffolds_from_patent_exclusive":          n_new_scafs,
        })

    # ── Assemble cliff_summary JSON ───────────────────────────────────────────
    cliff_json = {
        "n_compounds_analyzed":          n_analy,
        "n_fingerprints_computed":       n_fps,
        "n_fingerprint_failures":        len(failures),
        "similarity_landscape":          sim_counts,
        "cliff_counts":                  cliff_counts,
        "scaffold_analysis":             scaffold_stats,
        "patent_exclusive_contribution": patent_stats,
    }

    # ── Write ─────────────────────────────────────────────────────────────────
    print("\nWriting outputs...")
    write_outputs(pair_df, cliff_df, scaf_df, cliff_json, patent_json, args)

    print_summary(cliff_json, patent_json)
    print("\nStep 05 done.")


if __name__ == "__main__":
    main()
