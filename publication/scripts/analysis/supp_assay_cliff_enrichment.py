"""
supp_assay_cliff_enrichment.py — Supplementary Analysis: Assay Class Cliff Enrichment

For each mechanism class: compute % of all compounds vs % of severe cliff compounds.
Performs Fisher's Exact Test to check for significant enrichment/depletion.

Outputs:
  publication/figures/supplementary/fig_s02_assay_enrichment.{png,pdf}
  outputs/tables/supp_assay_enrichment.csv
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, save_fig, SEM, C, DOUBLE

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

set_style()

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT  = ROOT / 'publication/figures/supplementary'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("SUPP: Assay Class Cliff Enrichment Analysis")
print("=" * 60)

# ── Load data ────────────────────────────────────────────────────────────────
pad  = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
ac   = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
sev  = ac[ac['cliff_tier'] == 'severe'].copy()

sev_iks = set(sev['inchi_key_a'].tolist() + sev['inchi_key_b'].tolist())
n_total       = len(pad)
n_cliff_total = len(sev_iks)
print(f"  Total compounds:        {n_total}")
print(f"  Severe cliff compounds: {n_cliff_total}")

mech_order  = ['enzymatic', 'enzymatic_confirmed', 'fp_ic50', 'covalent']
mech_labels = ['Enzymatic\n(BAEE)', 'Enz. confirmed\n(RFMS)', 'FP binding', 'Covalent']
mech_colors = [SEM['enzymatic'], SEM['enzymatic_confirmed'], SEM['fp_ic50'], SEM['covalent']]

pad['in_cliff'] = pad['inchi_key'].isin(sev_iks)

results = []
for mech in mech_order:
    in_mech       = (pad['mechanism_class'] == mech)
    n_mech        = int(in_mech.sum())
    n_mech_cliff  = int((in_mech & pad['in_cliff']).sum())
    n_not_mech    = n_total - n_mech
    n_not_mech_cliff = n_cliff_total - n_mech_cliff

    # Fisher's exact test: 2x2 contingency table
    # [cliff in mech,    cliff NOT in mech]
    # [noncliff in mech, noncliff NOT in mech]
    contingency = [
        [n_mech_cliff,          n_not_mech_cliff],
        [n_mech - n_mech_cliff, n_not_mech - n_not_mech_cliff],
    ]
    odds_ratio, p_val = stats.fisher_exact(contingency)
    p_str = '< 0.001' if p_val < 0.001 else f'{p_val:.3f}'

    pct_all   = n_mech / n_total * 100
    pct_cliff = n_mech_cliff / n_cliff_total * 100 if n_cliff_total else 0

    print(f"  {mech}: n={n_mech} ({pct_all:.1f}% of all) | "
          f"in cliff={n_mech_cliff} ({pct_cliff:.1f}% of cliff) | OR={odds_ratio:.2f} p={p_str}")

    results.append({
        'Mechanism':              mech,
        'n all compounds':        n_mech,
        '% all compounds':        round(pct_all, 1),
        'n cliff compounds':      n_mech_cliff,
        '% cliff compounds':      round(pct_cliff, 1),
        'Odds Ratio':             round(odds_ratio, 3),
        "Fisher's Exact p-value": p_str,
    })

result_df = pd.DataFrame(results)
(ROOT / 'outputs/tables').mkdir(parents=True, exist_ok=True)
result_df.to_csv(ROOT / 'outputs/tables/supp_assay_enrichment.csv', index=False)

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(DOUBLE * 0.55, 2.8), constrained_layout=True)

x = np.arange(len(mech_order))
width = 0.35

pct_all_vals   = [r['% all compounds'] for r in results]
pct_cliff_vals = [r['% cliff compounds'] for r in results]

bars_all   = ax.bar(x - width/2, pct_all_vals,   width, color=[SEM['background']] * 4,
                     edgecolor='white', lw=0.4, label='% of all compounds (n=3,093)', alpha=0.85)
bars_cliff = ax.bar(x + width/2, pct_cliff_vals, width, color=mech_colors,
                     edgecolor='white', lw=0.4, label='% of severe cliff compounds (n=99)', alpha=0.85)

# Annotate p-values above each pair
for i, r in enumerate(results):
    p_str = r["Fisher's Exact p-value"]
    ax.text(i, max(pct_all_vals[i], pct_cliff_vals[i]) + 1.5,
            f'p={p_str}', ha='center', va='bottom', fontsize=5.5)

ax.set_xticks(x)
ax.set_xticklabels(mech_labels, fontsize=6)
ax.set_ylabel('Percentage (%)')
ax.set_title("Assay class representation: all compounds vs severe cliff compounds", fontsize=6, pad=3)
ax.legend(fontsize=5.5, loc='upper right', framealpha=0.85, edgecolor='none')
ax.set_ylim(0, max(max(pct_all_vals), max(pct_cliff_vals)) * 1.25)

save_fig(fig, str(OUT / 'fig_s02_assay_enrichment'))
plt.close(fig)

print(f"\nFigure saved: {OUT / 'fig_s02_assay_enrichment'}")
print(f"Table saved:  {ROOT / 'outputs/tables/supp_assay_enrichment.csv'}")
print("DONE")
