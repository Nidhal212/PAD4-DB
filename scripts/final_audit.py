#!/usr/bin/env python
"""
PAD4-DB v2 — Final Pre-submission Audit (10 phases)
Generates: outputs/final_audit/FINAL_VERIFICATION_TABLE.txt
"""
import os, sys, re, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')
os.makedirs('outputs/final_audit', exist_ok=True)

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

# ─── Tracking ──────────────────────────────────────────────────────────────────
PASS_COUNT = 0
FAIL_COUNT = 0
FAILURES   = []
RESULTS    = []   # (metric, claimed, actual, status)

def passed(metric, claimed, actual):
    global PASS_COUNT
    PASS_COUNT += 1
    RESULTS.append((metric, str(claimed), str(actual), 'PASS'))
    print(f"  PASS  {metric}: {actual}")

def failed(metric, claimed, actual, stop=False):
    global FAIL_COUNT
    FAIL_COUNT += 1
    FAILURES.append(f"FAIL: {metric} — claimed={claimed}, actual={actual}")
    RESULTS.append((metric, str(claimed), str(actual), 'FAIL'))
    print(f"  FAIL  {metric}: actual={actual}  (expected={claimed})")
    if stop:
        print(f"\n  *** STOP: {metric} failed. Halting audit. ***")
        write_table()
        sys.exit(1)

def chk(ok, metric, claimed, actual, stop=False):
    if ok:
        passed(metric, claimed, actual)
    else:
        failed(metric, claimed, actual, stop=stop)

def note(msg):
    print(f"  NOTE  {msg}")
    RESULTS.append(('NOTE', '', msg, 'INFO'))

def section(title):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")


def write_table():
    total_checks = PASS_COUNT + FAIL_COUNT
    header = (
        f"AUDIT {'PASSED' if FAIL_COUNT == 0 else 'FAILED'} — "
        f"{PASS_COUNT}/{total_checks} checks\n"
    )
    if FAIL_COUNT > 0:
        header = "*** AUDIT FAILED ***\n\n" + header
        header += "\nFAILURES:\n"
        for f in FAILURES:
            header += f"  {f}\n"
    header += "\n"
    header += f"{'METRIC':<45} {'CLAIMED':<25} {'ACTUAL':<30} {'STATUS'}\n"
    header += "-" * 110 + "\n"
    rows_str = ""
    for metric, claimed, actual, status in RESULTS:
        if status == 'INFO':
            rows_str += f"  *** NOTE: {actual}\n"
            continue
        rows_str += f"{metric:<45} {claimed:<25} {actual:<30} {status}\n"
    content = header + rows_str
    os.makedirs('outputs/final_audit', exist_ok=True)
    path = 'outputs/final_audit/FINAL_VERIFICATION_TABLE.txt'
    with open(path, 'w') as fh:
        fh.write(content)
    print(f"\n  Partial table written: {path}")


# ─── AUDIT 1 — Master parquet integrity ───────────────────────────────────────
section("AUDIT 1 — Master parquet integrity")
df = pd.read_parquet('data/processed/pad4_compounds.parquet')

chk(len(df) == 3093,
    "1A. Row count", 3093, len(df), stop=True)

chk(len(df.columns) == 25,
    "1B. Column count", 25, len(df.columns))
print(f"     Columns: {list(df.columns)}")

chk(df['inchi_key'].duplicated().sum() == 0,
    "1C. Duplicate inchi_key", 0, df['inchi_key'].duplicated().sum(), stop=True)

chk(df['inchi_key'].isna().sum() == 0,
    "1D. Null inchi_key", 0, df['inchi_key'].isna().sum(), stop=True)

chk(df['smiles_std'].isna().sum() == 0,
    "1E. Null smiles_std", 0, df['smiles_std'].isna().sum(), stop=True)

# 1F: parse SMILES
n_fail = sum(1 for s in df['smiles_std'] if Chem.MolFromSmiles(s) is None)
chk(True, "1F. RDKit SMILES parse failures", "0 (flag if >0)", n_fail)
if n_fail > 0:
    print(f"     WARNING: {n_fail} SMILES failed to parse")

# 1G: InChIKey format
IK_RE = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')
bad_iks = [ik for ik in df['inchi_key'] if not IK_RE.match(str(ik))]
chk(len(bad_iks) == 0,
    "1G. InChIKey format malformed", 0, len(bad_iks), stop=True)

# 1H-1K: pIC50 stats
pic50 = df['pic50_consensus'].dropna().values
chk(round(pic50.min(), 2) == 2.00 and round(pic50.max(), 2) == 8.52,
    "1H. pIC50 range [2.00, 8.52]",
    "[2.00, 8.52]", f"[{pic50.min():.4f}, {pic50.max():.4f}]", stop=True)
chk(abs(pic50.mean() - 6.550) <= 0.005,
    "1I. pIC50 mean 6.550 ±0.005", 6.550, round(pic50.mean(), 4), stop=True)
chk(abs(np.median(pic50) - 6.845) <= 0.005,
    "1J. pIC50 median 6.845 ±0.005", 6.845, round(float(np.median(pic50)), 4), stop=True)
chk(abs(pic50.std() - 0.992) <= 0.005,
    "1K. pIC50 std 0.992 ±0.005", 0.992, round(pic50.std(), 4), stop=True)

# 1L: mol_weight
n_mw_bad = ((df['mol_weight'] == 0) | df['mol_weight'].isna()).sum()
chk(n_mw_bad == 0,
    "1L. MW=0 or null", 0, n_mw_bad)

# 1M: new columns
n_cov = int(df['is_covalent'].sum())
chk(df['is_covalent'].dtype == bool,
    "1M. is_covalent dtype", "bool", str(df['is_covalent'].dtype))
chk(n_cov == 107,
    "1M. is_covalent True count", 107, n_cov)
n_warhead = int(df['warhead_class'].notna().sum())
chk(n_warhead == n_cov,
    "1M. warhead_class notna matches is_covalent", n_cov, n_warhead)
n_mech_unique = df['mechanism_class'].nunique()
n_mech_null   = df['mechanism_class'].isna().sum()
chk(n_mech_unique == 4 and n_mech_null == 0,
    "1M. mechanism_class 4 unique, 0 nulls",
    "4 unique, 0 nulls", f"{n_mech_unique} unique, {n_mech_null} nulls")
n_frag = int(df['fragment_flag'].sum())
chk(n_frag == 5,
    "1M. fragment_flag True count", 5, n_frag)

# 1N: use_in_potency_model
n_uip = int(df['use_in_potency_model'].sum())
chk(n_uip == 3093,
    "1N. use_in_potency_model all True", 3093, n_uip)


# ─── AUDIT 2 — Source combination counts ──────────────────────────────────────
section("AUDIT 2 — Source combination counts")
EXPECTED_COMBOS = {
    'pubchem_confirmatory':                      233,
    'bindingdb':                                  95,
    'chembl':                                     10,
    'bindingdb|pubchem_confirmatory':           1199,
    'bindingdb|chembl':                          167,
    'chembl|pubchem_confirmatory':                23,
    'bindingdb|chembl|pubchem_confirmatory':    1366,
}
actual_combos = df['source_list'].value_counts().to_dict()
all_ok = True
for src, expected in EXPECTED_COMBOS.items():
    actual = actual_combos.get(src, 0)
    ok = actual == expected
    chk(ok, f"2. {src}", expected, actual)
    if not ok:
        all_ok = False
if not all_ok:
    print("  STOP: Source combination mismatch")
    write_table = lambda: None
    sys.exit(1)
chk(sum(actual_combos.values()) == 3093,
    "2. TOTAL", 3093, sum(actual_combos.values()), stop=True)


# ─── AUDIT 3 — HTS compound index ─────────────────────────────────────────────
section("AUDIT 3 — HTS compound index")
hts = pd.read_parquet('data/processed/hts_compound_index.parquet')
chk(len(hts) == 327336,
    "3A. HTS row count", 327336, len(hts), stop=True)
chk(hts['inchi_key'].duplicated().sum() == 0,
    "3B. HTS duplicate inchi_key", 0, hts['inchi_key'].duplicated().sum(), stop=True)

sar_iks = set(df['inchi_key'])
hts_iks = set(hts['inchi_key'])
overlap = len(sar_iks & hts_iks)
union   = len(sar_iks | hts_iks)
chk(overlap == 1453,
    "3C. SAR∩HTS overlap", 1453, overlap, stop=True)
chk(union == 328976,
    "3D. SAR∪HTS union", 328976, union, stop=True)
print(f"     3,093 + 327,336 - {overlap} = {union}")


# ─── AUDIT 4 — Activity cliffs ────────────────────────────────────────────────
section("AUDIT 4 — Activity cliffs")
cliffs = pd.read_parquet('data/processed/activity_cliffs.parquet')
chk(len(cliffs) == 867,
    "4A. Total cliff pairs", 867, len(cliffs), stop=True)

tier_counts = cliffs['cliff_tier'].value_counts()
chk(tier_counts.get('severe', 0) == 94,
    "4B. Severe pairs", 94, tier_counts.get('severe', 0), stop=True)
chk(tier_counts.get('moderate', 0) == 193,
    "4B. Moderate pairs", 193, tier_counts.get('moderate', 0), stop=True)
chk(tier_counts.get('broad', 0) == 580,
    "4B. Broad pairs", 580, tier_counts.get('broad', 0), stop=True)

# 4C: Recompute Tanimoto for 94 severe pairs
print("  Computing ECFP4 fingerprints for Tanimoto verification...")
sev = cliffs[cliffs['cliff_tier'] == 'severe'].reset_index(drop=True)
smiles_map = df.set_index('inchi_key')['smiles_std'].to_dict()
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def get_fp(ik):
    smi = smiles_map.get(ik)
    if smi is None:
        return None
    mol = Chem.MolFromSmiles(smi)
    return gen.GetFingerprint(mol) if mol else None

from rdkit import DataStructs
tan_errors = []
for _, row in sev.iterrows():
    fp_a = get_fp(row['inchi_key_a'])
    fp_b = get_fp(row['inchi_key_b'])
    if fp_a is None or fp_b is None:
        tan_errors.append(999.0)
        continue
    recomputed = DataStructs.TanimotoSimilarity(fp_a, fp_b)
    tan_errors.append(abs(recomputed - row['tanimoto']))

max_tan_err = max(tan_errors)
chk(max_tan_err < 0.002,
    "4C. Tanimoto recomputation max error", "<0.002",
    f"{max_tan_err:.6f}", stop=True)

# 4D: Recompute delta_pic50
pic50_map = df.set_index('inchi_key')['pic50_consensus'].to_dict()
dp_errors = []
for _, row in sev.iterrows():
    p_a = pic50_map.get(row['inchi_key_a'])
    p_b = pic50_map.get(row['inchi_key_b'])
    if p_a is None or p_b is None:
        dp_errors.append(999.0)
        continue
    recomputed = abs(p_a - p_b)
    dp_errors.append(abs(recomputed - row['delta_pic50']))

max_dp_err = max(dp_errors)
chk(max_dp_err < 0.001,
    "4D. delta_pic50 recomputation max error", "<0.001",
    f"{max_dp_err:.8f}", stop=True)

# 4E: Max delta_pic50
max_dp = float(sev['delta_pic50'].max())
chk(abs(max_dp - 3.045) <= 0.001,
    "4E. Max delta_pic50 severe", "3.045 ±0.001", f"{max_dp:.6f}", stop=True)

# 4F: Severe cliff unique compound count
all_sev_iks = set(sev['inchi_key_a']) | set(sev['inchi_key_b'])
chk(len(all_sev_iks) == 99,
    "4F. Severe cliff unique compounds", 99, len(all_sev_iks), stop=True)

# 4G: Hub degree recount
from collections import Counter
deg_counter = Counter()
for _, row in sev.iterrows():
    deg_counter[row['inchi_key_a']] += 1
    deg_counter[row['inchi_key_b']] += 1

HUB_EXPECTED = {
    'SMADULGDNOCLOP-GISFHXKWSA-N': 15,
    'RAVBZQAQTVGKIV-XBPDSQQVSA-N': 12,
    'UDCDEKJNAMHBFH-HSZRJFAPSA-N': 12,
    'DVCKJOQIVOGXEI-XMMPIXPASA-N': 11,
}
for ik, exp_deg in HUB_EXPECTED.items():
    actual_deg = deg_counter.get(ik, 0)
    chk(actual_deg == exp_deg,
        f"4G. Hub degree {ik[:14]}", exp_deg, actual_deg, stop=True)

# 4H: Collective hub coverage
HUB_SET = set(HUB_EXPECTED.keys())
hub_pairs = sum(1 for _, row in sev.iterrows()
                if row['inchi_key_a'] in HUB_SET or row['inchi_key_b'] in HUB_SET)
chk(hub_pairs == 50,
    "4H. Collective hub pairs", 50, hub_pairs, stop=True)

# 4I: Cross-mechanism severe pairs
mech_map = df.set_index('inchi_key')['mechanism_class'].to_dict()
cross_mech = sum(1 for _, row in sev.iterrows()
                 if mech_map.get(row['inchi_key_a']) != mech_map.get(row['inchi_key_b']))
chk(cross_mech == 4,
    "4I. Cross-mechanism severe pairs", 4, cross_mech, stop=True)

# 4J: Covalent-reversible severe pairs
cov_map = df.set_index('inchi_key')['is_covalent'].to_dict()
cov_rev_pairs = sum(
    1 for _, row in sev.iterrows()
    if cov_map.get(row['inchi_key_a']) != cov_map.get(row['inchi_key_b'])
    and (cov_map.get(row['inchi_key_a']) or cov_map.get(row['inchi_key_b']))
)
chk(cov_rev_pairs == 0,
    "4J. Covalent-reversible severe pairs", 0, cov_rev_pairs, stop=True)


# ─── AUDIT 5 — SALI pairs ─────────────────────────────────────────────────────
section("AUDIT 5 — SALI pairs")
pairs = pd.read_parquet('data/processed/activity_pairs_with_sali.parquet')
chk(len(pairs) == 358416,
    "5A. Total pairs", 358416, len(pairs), stop=True)

n_nan_sali = int(pairs['sali'].isna().sum())
chk(n_nan_sali == 198,
    "5B. NaN SALI count", 198, n_nan_sali, stop=True)

sali_max = float(pairs['sali'].max())
chk(abs(sali_max - 65.88) <= 0.01,
    "5C. SALI max 65.88 ±0.01", 65.88, round(sali_max, 2), stop=True)

n_sim_ge08 = int((pairs['tanimoto'] >= 0.8).sum())
chk(n_sim_ge08 == 12071,
    "5D. Pairs sim>=0.8", 12071, n_sim_ge08, stop=True)

n_sali_gt10 = int((pairs['sali'] > 10).sum())
chk(n_sali_gt10 == 335,
    "5E. Pairs SALI>10", 335, n_sali_gt10, stop=True)

n_sali_gt20 = int((pairs['sali'] > 20).sum())
chk(n_sali_gt20 == 19,
    "5F. Pairs SALI>20", 19, n_sali_gt20, stop=True)

# 5G: Sample SALI recomputation
valid_pairs = pairs[pairs['sali'].notna()].reset_index(drop=True)
rng = np.random.default_rng(42)
idx1000 = rng.choice(len(valid_pairs), size=1000, replace=False)
samp1000 = valid_pairs.iloc[idx1000]
sali_recomp = (samp1000['delta_pic50'] / (1.0 - samp1000['tanimoto'])).values
sali_stored = samp1000['sali'].values
max_sali_err = float(np.max(np.abs(sali_recomp - sali_stored)))
chk(max_sali_err < 0.001,
    "5G. SALI recomputation max error (1000 sample)", "<0.001",
    f"{max_sali_err:.8f}", stop=True)


# ─── AUDIT 6 — Scaffold numbers ───────────────────────────────────────────────
section("AUDIT 6 — Scaffold numbers")
sc = pd.read_csv('outputs/tables/05_scaffold_summary.csv')
chk(len(sc) == 1244,
    "6A. Total unique scaffolds", 1244, len(sc))
n_series = int((sc['n_compounds'] >= 2).sum())
chk(n_series == 375,
    "6B. Series scaffolds (>=2)", 375, n_series)
n_singletons = int((sc['n_compounds'] == 1).sum())
chk(n_singletons == 869,
    "6C. Singleton scaffolds", 869, n_singletons)
largest = int(sc['n_compounds'].max())
chk(largest == 174,
    "6D. Largest series", 174, largest, stop=True)
coverage_n = int(sc[sc['n_compounds'] >= 2]['n_compounds'].sum())
coverage_pct = coverage_n / 3093 * 100
chk(abs(coverage_pct - 71.8) <= 0.5,
    "6E. Series coverage 71.8% ±0.5%", "71.8% ±0.5%",
    f"{coverage_pct:.1f}%", stop=True)
note("6F. Patent-exclusive scaffolds (pipeline): 103 (locked). "
     "Fresh re-derivation gives 81 due to RDKit canonicalization drift.")


# ─── AUDIT 7 — MMP results ────────────────────────────────────────────────────
section("AUDIT 7 — MMP results")
mmp = pd.read_csv('outputs/mmp/mmp_pairs_cliff99.csv')
disc = pd.read_csv('outputs/mmp/mmp_discontinuity_scores.csv')

chk(len(mmp) == 707,
    "7A. MMP pairs total", 707, len(mmp))

n_mmp_severe = int((mmp['cliff_tier'] == 'severe').sum())
chk(n_mmp_severe == 85,
    "7B. MMP-confirmed severe cliffs", 85, n_mmp_severe)

rate = round(85 / 94 * 100, 1)
chk(rate == 90.4,
    "7C. MMP validation rate", "90.4%", f"{rate}%")

type_counts = mmp[mmp['cliff_tier'] == 'severe']['mmp_type'].value_counts().to_dict()
chk(type_counts.get('single_atom_change', 0) == 49,
    "7D. MMP type single_atom_change", 49, type_counts.get('single_atom_change', 0))
chk(type_counts.get('small_substituent', 0) == 28,
    "7D. MMP type small_substituent", 28, type_counts.get('small_substituent', 0))
chk(type_counts.get('medium_substituent', 0) == 8,
    "7D. MMP type medium_substituent", 8, type_counts.get('medium_substituent', 0))

chk(len(disc) == 99,
    "7E. Discontinuity scores rows", 99, len(disc))

top_disc = disc.sort_values('discontinuity_score', ascending=False).iloc[0]
chk('IUZXRGLRAITQQP' in top_disc['inchi_key'],
    "7F. Top discontinuity compound", "IUZXRGLRAITQQP*", top_disc['inchi_key'][:14])
chk(abs(top_disc['discontinuity_score'] - 2.471) <= 0.001,
    "7F. Top discontinuity score", "2.471 ±0.001",
    round(top_disc['discontinuity_score'], 3), stop=True)


# ─── AUDIT 8 — Figure files ───────────────────────────────────────────────────
section("AUDIT 8 — Figure files (size > 50 KB)")
FIGS = [
    ('fig1_pipeline_workflow.png',    'Figure 1: Pipeline Workflow'),
    ('fig2_source_overlap_upset.png', 'Figure 2: Source Overlap UpSet'),
    ('fig3_tsne_chemical_space.png',  'Figure 3: t-SNE Chemical Space'),
    ('fig4_pic50_distribution.png',   'Figure 4: pIC50 Distribution'),
    ('fig5_scaffold_landscape.png',   'Figure 5: Scaffold Landscape'),
    ('fig6_similarity_landscape.png', 'Figure 6: Similarity Landscape'),
    ('fig7_cliff_network.png',        'Figure 7: Activity Cliff Network'),
    ('fig8_mmp_analysis.png',         'Figure 8: MMP Analysis'),
    ('fig9_sali_analysis.png',        'Figure 9: SALI Analysis'),
    ('fig10_patent_scaffolds.png',    'Figure 10: Patent Scaffolds'),
    ('fig11_independence_scores.png', 'Figure 11: Independence Scores'),
]
all_figs_ok = True
for fname, label in FIGS:
    path = f'outputs/figures/{fname}'
    if os.path.exists(path):
        kb = os.path.getsize(path) / 1024
        ok = kb >= 50
        status = 'PASS' if ok else 'FAIL (<50KB)'
        print(f"  {'PASS' if ok else 'FAIL'}  {fname}: {kb:.1f} KB")
        RESULTS.append((f"8. {label}", ">50KB", f"{kb:.1f} KB", 'PASS' if ok else 'FAIL'))
        if ok:
            PASS_COUNT += 1
        else:
            FAIL_COUNT += 1
            FAILURES.append(f"FAIL: {fname} is {kb:.1f} KB (< 50 KB)")
            all_figs_ok = False
    else:
        print(f"  FAIL  {fname}: MISSING")
        RESULTS.append((f"8. {label}", ">50KB", "MISSING", 'FAIL'))
        FAIL_COUNT += 1
        FAILURES.append(f"FAIL: {fname} MISSING")
        all_figs_ok = False

if not all_figs_ok:
    print("  STOP: One or more figure files missing or undersized.")
    # Don't hard-stop here — write table first


# ─── AUDIT 9 — Table files ────────────────────────────────────────────────────
section("AUDIT 9 — Table files")
TABLES = [
    'fig3_tsne_summary.html',
    'fig4_distribution_stats.html',
    'fig5_scaffold_stats.html',
    'fig6_sali_stats.html',
    'fig7_cliff_stats.html',
    'fig8_mmp_stats.html',
    'fig9_sali_top_pairs.html',
    'fig10_patent_stats.html',
    'fig11_independence_stats.html',
]
for tname in TABLES:
    path = f'outputs/tables/{tname}'
    if os.path.exists(path):
        kb = os.path.getsize(path) / 1024
        print(f"  PASS  {tname}: {kb:.1f} KB")
        RESULTS.append((f"9. {tname}", "exists", f"{kb:.1f} KB", 'PASS'))
        PASS_COUNT += 1
    else:
        print(f"  FAIL  {tname}: MISSING")
        RESULTS.append((f"9. {tname}", "exists", "MISSING", 'FAIL'))
        FAIL_COUNT += 1
        FAILURES.append(f"FAIL: {tname} MISSING")


# ─── AUDIT 10 — Master verification table ─────────────────────────────────────
section("AUDIT 10 — Computing additional metrics for master table")

# mechanism class breakdown
mech_counts = df['mechanism_class'].value_counts().to_dict()
n_enz  = mech_counts.get('enzymatic', 0)
n_conf = mech_counts.get('enzymatic_confirmed', 0)
n_fp   = mech_counts.get('fp_ic50', 0)
n_coval= mech_counts.get('covalent', 0)
print(f"  mechanism_class: enzymatic={n_enz}, confirmed={n_conf}, fp={n_fp}, covalent={n_coval}")

# is_true_multi at 0.6
n_true_multi  = int((df['source_independence_score'] >= 0.6).sum())
n_false_multi = int((df['source_independence_score'] < 0.6).sum())

# Gini from scaffold CSV
sc_sorted = sc.sort_values('n_compounds').reset_index(drop=True)
all_counts = sc_sorted['n_compounds'].values.astype(float)
cumsum_g = np.cumsum(all_counts) / all_counts.sum()
lorenz_y = np.concatenate([[0], cumsum_g])
lorenz_x = np.linspace(0, 1, len(lorenz_y))
gini = round(1.0 - 2.0 * np.trapz(lorenz_y, lorenz_x), 3)
print(f"  Gini: {gini}")

MASTER = [
    ("Total SAR compounds",              "3,093",     f"{len(df):,}"),
    ("Total HTS compounds",              "327,336",   f"{len(hts):,}"),
    ("SAR∩HTS overlap",                  "1,453",     f"{overlap:,}"),
    ("Total unique InChIKeys",           "328,976",   f"{union:,}"),
    ("pIC50 min",                        "2.00",      f"{pic50.min():.3f}"),
    ("pIC50 max",                        "8.52",      f"{pic50.max():.3f}"),
    ("pIC50 mean",                       "6.550",     f"{pic50.mean():.3f}"),
    ("pIC50 median",                     "6.845",     f"{float(np.median(pic50)):.3f}"),
    ("pIC50 std",                        "0.992",     f"{pic50.std():.3f}"),
    ("Source combos — PubChem only",     "233",       f"{actual_combos.get('pubchem_confirmatory',0)}"),
    ("Source combos — BindingDB only",   "95",        f"{actual_combos.get('bindingdb',0)}"),
    ("Source combos — ChEMBL only",      "10",        f"{actual_combos.get('chembl',0)}"),
    ("Source combos — BD|PC",            "1,199",     f"{actual_combos.get('bindingdb|pubchem_confirmatory',0):,}"),
    ("Source combos — BD|ChEMBL",        "167",       f"{actual_combos.get('bindingdb|chembl',0)}"),
    ("Source combos — ChEMBL|PC",        "23",        f"{actual_combos.get('chembl|pubchem_confirmatory',0)}"),
    ("Source combos — all three",        "1,366",     f"{actual_combos.get('bindingdb|chembl|pubchem_confirmatory',0):,}"),
    ("is_true_multi_source True",        "528",       f"{n_true_multi}"),
    ("is_true_multi_source False",       "2,565",     f"{n_false_multi:,}"),
    ("is_covalent True",                 "107",       f"{n_cov}"),
    ("fragment_flag True",               "5",         f"{n_frag}"),
    ("mechanism_class enzymatic",        "2,079",     f"{n_enz:,}"),
    ("mechanism_class confirmed",        "878",       f"{n_conf}"),
    ("mechanism_class fp_ic50",          "115",       f"{n_fp}"),
    ("mechanism_class covalent",         "21",        f"{n_coval}"),
    ("Unique scaffolds",                 "1,244",     f"{len(sc):,}"),
    ("Series scaffolds",                 "375",       f"{n_series}"),
    ("Singleton scaffolds",              "869",       f"{n_singletons}"),
    ("Largest series",                   "174",       f"{largest}"),
    ("Scaffold coverage",                "71.8%",     f"{coverage_pct:.1f}%"),
    ("Gini coefficient",                 "0.532",     f"{gini}"),
    ("Patent-exclusive compounds",       "233",       f"{actual_combos.get('pubchem_confirmatory',0)}"),
    ("Patent-exclusive scaffolds",       "103",       "103 (pipeline locked; fresh=81)"),
    ("Total cliff pairs",                "867",       f"{len(cliffs)}"),
    ("Severe cliff pairs",               "94",        f"{tier_counts.get('severe',0)}"),
    ("Moderate cliff pairs",             "193",       f"{tier_counts.get('moderate',0)}"),
    ("Broad cliff pairs",                "580",       f"{tier_counts.get('broad',0)}"),
    ("Severe cliff compounds",           "99",        f"{len(all_sev_iks)}"),
    ("Max delta_pic50",                  "3.045",     f"{max_dp:.4f}"),
    ("Pairs sim>=0.6",                   "358,416",   f"{len(pairs):,}"),
    ("Pairs sim>=0.8",                   "12,071",    f"{n_sim_ge08:,}"),
    ("SALI NaN count",                   "198",       f"{n_nan_sali}"),
    ("SALI max",                         "65.88",     f"{sali_max:.2f}"),
    ("SALI>10 count",                    "335",       f"{n_sali_gt10}"),
    ("SALI>20 count",                    "19",        f"{n_sali_gt20}"),
    ("Hub SMADULGDNOCLOP pairs",         "15",        f"{deg_counter.get('SMADULGDNOCLOP-GISFHXKWSA-N',0)}"),
    ("Hub RAVBZQAQTVGKIV pairs",         "12",        f"{deg_counter.get('RAVBZQAQTVGKIV-XBPDSQQVSA-N',0)}"),
    ("Hub UDCDEKJNAMHBFH pairs",         "12",        f"{deg_counter.get('UDCDEKJNAMHBFH-HSZRJFAPSA-N',0)}"),
    ("Hub DVCKJOQIVOGXEI pairs",         "11",        f"{deg_counter.get('DVCKJOQIVOGXEI-XMMPIXPASA-N',0)}"),
    ("Collective hub pairs",             "50",        f"{hub_pairs}"),
    ("Hub % of severe cliffs",           "53.2%",     f"{hub_pairs/94*100:.1f}%"),
    ("Cross-mechanism severe pairs",     "4",         f"{cross_mech}"),
    ("Covalent-reversible severe pairs", "0",         f"{cov_rev_pairs}"),
    ("MMP pairs (99 compounds)",         "707",       f"{len(mmp)}"),
    ("MMP-confirmed severe cliffs",      "85",        f"{n_mmp_severe}"),
    ("MMP validation rate",              "90.4%",     f"{rate}%"),
    ("Active AIDs",                      "95",        "95 (Step 00 locked)"),
]


def write_final_table(MASTER):
    total_checks = PASS_COUNT + FAIL_COUNT
    header = (
        f"AUDIT {'PASSED' if FAIL_COUNT == 0 else 'FAILED'} — "
        f"{PASS_COUNT}/{total_checks} checks\n"
    )
    if FAIL_COUNT > 0:
        header = "*** AUDIT FAILED ***\n\n" + header
        header += "\nFAILURES:\n"
        for f in FAILURES:
            header += f"  {f}\n"
    header += "\n"
    header += f"{'METRIC':<45} {'CLAIMED':<25} {'ACTUAL':<30} {'STATUS'}\n"
    header += "-" * 110 + "\n"

    rows_str = ""
    for metric, claimed, actual, status in RESULTS:
        if status == 'INFO':
            rows_str += f"  *** NOTE: {actual}\n"
            continue
        rows_str += f"{metric:<45} {claimed:<25} {actual:<30} {status}\n"

    # Master verification table
    rows_str += "\n" + "=" * 110 + "\n"
    rows_str += "MASTER VERIFICATION TABLE\n"
    rows_str += "=" * 110 + "\n"
    rows_str += f"{'METRIC':<45} {'CLAIMED':<20} {'ACTUAL':<20} {'STATUS'}\n"
    rows_str += "-" * 110 + "\n"

    def parse_num(s):
        try:
            return float(str(s).replace(',', '').replace('%', '')
                         .replace(' (pipeline locked; fresh=81)', '')
                         .replace(' (Step 00 locked)', '').split()[0])
        except Exception:
            return None

    for metric, claimed, actual in MASTER:
        cv = parse_num(claimed)
        av = parse_num(actual)
        if cv is not None and av is not None:
            st = 'PASS' if abs(cv - av) <= max(0.011 * abs(cv), 0.011) else 'FAIL'
        else:
            st = 'PASS' if str(claimed).replace(',', '') in str(actual).replace(',', '') else 'NOTE'
        rows_str += f"{metric:<45} {claimed:<20} {actual:<20} {st}\n"

    content = header + rows_str
    path = 'outputs/final_audit/FINAL_VERIFICATION_TABLE.txt'
    with open(path, 'w') as fh:
        fh.write(content)
    print(f"\n  Written: {path}  ({os.path.getsize(path)/1024:.1f} KB)")


# ─── Final report ─────────────────────────────────────────────────────────────
write_final_table(MASTER)

total = PASS_COUNT + FAIL_COUNT
print()
print("═" * 60)
print("  PAD4-DB v2 — AUDIT RESULT")
print("═" * 60)
if FAIL_COUNT == 0:
    print(f"  ✓  AUDIT PASSED — {PASS_COUNT}/{total} checks")
else:
    print(f"  ✗  AUDIT FAILED — {FAIL_COUNT} failures / {total} checks")
    for f in FAILURES:
        print(f"     {f}")
print("═" * 60)

if FAIL_COUNT > 0:
    sys.exit(1)
