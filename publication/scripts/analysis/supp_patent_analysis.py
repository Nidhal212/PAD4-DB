"""
supp_patent_analysis.py — Supplementary Analysis: Patent vs Literature Cliff Contribution

Computes odds ratio and effect size of patent compounds appearing in severe cliff pairs.
Framed as exploratory; statistical power is reported where limited.

Outputs:
  outputs/tables/supp_patent_cliff_odds.csv
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

ROOT = Path('/home/nidhal/PAD4-db_V2')

print("=" * 60)
print("SUPP: Patent vs Literature Cliff Enrichment (Exploratory)")
print("=" * 60)

# ── Load data ────────────────────────────────────────────────────────────────
pad  = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
ac   = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
sev  = ac[ac['cliff_tier'] == 'severe'].copy()

# Patent flag per InChIKey
patent_map = dict(zip(pad['inchi_key'], pad['patent_flag']))

sev_iks = set(sev['inchi_key_a'].tolist() + sev['inchi_key_b'].tolist())
n_total       = len(pad)
n_cliff_total = len(sev_iks)

# Patent status
n_patent       = int(pad['patent_flag'].sum())
n_nonpatent    = n_total - n_patent
n_pat_cliff    = int(sum(1 for ik in sev_iks if patent_map.get(ik, False)))
n_nonpat_cliff = n_cliff_total - n_pat_cliff

print(f"  Total compounds:          {n_total}")
print(f"  Patent-exclusive:         {n_patent} ({n_patent/n_total*100:.1f}%)")
print(f"  Severe cliff compounds:   {n_cliff_total}")
print(f"  Patent in cliff:          {n_pat_cliff} ({n_pat_cliff/n_cliff_total*100:.1f}%)")
print(f"  Non-patent in cliff:      {n_nonpat_cliff}")

# 2×2 contingency table for Fisher's exact test
#           | In cliff | Not in cliff |
# Patent    |   A      |      B       |
# Non-patent|   C      |      D       |
A = n_pat_cliff
B = n_patent - n_pat_cliff
C = n_nonpat_cliff
D = n_nonpatent - n_nonpat_cliff

contingency = [[A, B], [C, D]]
odds_ratio, p_val = stats.fisher_exact(contingency)
p_str = '< 0.001' if p_val < 0.001 else f'{p_val:.4f}'

print(f"\n  Contingency table:")
print(f"    Patent/cliff={A} | Patent/non-cliff={B}")
print(f"    NonPat/cliff={C} | NonPat/non-cliff={D}")
print(f"  Odds Ratio: {odds_ratio:.3f}")
print(f"  Fisher's p: {p_str}")

# Cliff pair-level analysis: any_patent_exclusive column
if 'any_patent_exclusive' in sev.columns:
    n_pairs_patent = int(sev['any_patent_exclusive'].sum())
    pct_pairs = n_pairs_patent / len(sev) * 100
    print(f"\n  Cliff pairs involving ≥1 patent compound: {n_pairs_patent} / {len(sev)} ({pct_pairs:.1f}%)")
    print(f"  Expected by chance (patent={n_patent/n_total*100:.1f}% of compounds, "
          f"both compounds independent): {n_patent/n_total*100:.1f}%")
else:
    n_pairs_patent = 0
    pct_pairs = 0.0

# Effect size (phi coefficient for 2x2)
n_total_4 = A + B + C + D
chi2 = (A*D - B*C)**2 * n_total_4 / ((A+B) * (C+D) * (A+C) * (B+D))
phi  = np.sqrt(chi2 / n_total_4)
print(f"\n  Effect size (phi): {phi:.3f}")
print(f"  Note: phi < 0.1 = negligible, 0.1–0.3 = small, 0.3–0.5 = medium, > 0.5 = large")

# Power note
print(f"\n  Statistical power note: n_patent_cliff = {n_pat_cliff} compounds.")
if n_pat_cliff < 20:
    print("  WARNING: Low cell count (< 20). Fisher's test is valid but power is limited.")
    print("  Interpret odds ratio with caution; report as exploratory.")

# ── Save output ──────────────────────────────────────────────────────────────
result = {
    'n_total_compounds':        n_total,
    'n_patent_compounds':       n_patent,
    'n_severe_cliff_compounds': n_cliff_total,
    'n_patent_in_cliff':        n_pat_cliff,
    'n_nonpatent_in_cliff':     n_nonpat_cliff,
    'n_cliff_pairs_total':      len(sev),
    'n_cliff_pairs_any_patent': int(n_pairs_patent),
    'pct_cliff_pairs_patent':   round(pct_pairs, 1),
    'odds_ratio':               round(odds_ratio, 4),
    'fishers_exact_p':          p_str,
    'phi_effect_size':          round(phi, 4),
    'note':                     'Exploratory; statistical power limited by small patent cliff cell count.',
}
df_out = pd.DataFrame([result])
out_dir = ROOT / 'outputs/tables'
out_dir.mkdir(parents=True, exist_ok=True)
df_out.to_csv(out_dir / 'supp_patent_cliff_odds.csv', index=False)

print(f"\nTable saved: {out_dir / 'supp_patent_cliff_odds.csv'}")
print("DONE")
