"""
fig_s03_sas_map.py — Supplementary Figure S10: Structure-Activity Similarity (SAS) map

The canonical activity-cliff landscape (Maggiora / SAS map):
  Panel a: pairwise Tanimoto similarity (x) vs |ΔpIC50| (y) for all 358,416
           structurally related pairs (Tanimoto ≥ 0.6). Hexbin density (log scale).
           Cliff thresholds (Tanimoto = 0.8, |ΔpIC50| = 2.0) drawn; the upper-right
           "activity cliff" quadrant is shaded. Severe cliff pairs highlighted.
  Panel b: |ΔpIC50| distribution stratified by similarity bin — demonstrates
           "diagonal absence": highly similar pairs are overwhelmingly concordant,
           with only a thin activity-cliff tail.

Outputs: publication/figures/supplementary/fig_s03_sas_map.{png,pdf}
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, panel_label, save_fig, SEM, C, DOUBLE

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm
from pathlib import Path

set_style()

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT  = ROOT / 'publication/figures/supplementary'
OUT.mkdir(parents=True, exist_ok=True)

TAN_CUT, DELTA_CUT = 0.8, 2.0

print("=" * 60)
print("FIGURE S10 — STRUCTURE-ACTIVITY SIMILARITY (SAS) MAP")
print("=" * 60)

pairs = pd.read_parquet(ROOT / 'data/processed/activity_pairs_with_sali.parquet')
pairs['abs_delta'] = pairs['delta_pic50'].abs()
print(f"  Total related pairs (Tan≥0.6): {len(pairs):,}")

sev = pairs[(pairs['tanimoto'] >= TAN_CUT) & (pairs['abs_delta'] >= DELTA_CUT)]
print(f"  Severe cliff pairs (Tan≥0.8, Δ≥2.0): {len(sev)}")

# ── SAS quadrant summary table ────────────────────────────────────────────────
hi_sim = pairs['tanimoto'] >= TAN_CUT
hi_dlt = pairs['abs_delta'] >= DELTA_CUT
quad = {
    'Activity cliffs (high sim, high Δ)':       int((hi_sim & hi_dlt).sum()),
    'Smooth/continuous SAR (high sim, low Δ)':  int((hi_sim & ~hi_dlt).sum()),
    'Discontinuous (low sim, high Δ)':          int((~hi_sim & hi_dlt).sum()),
    'Non-descript (low sim, low Δ)':            int((~hi_sim & ~hi_dlt).sum()),
}
tot = len(pairs)
qdf = pd.DataFrame([
    {'SAS quadrant': k, 'N pairs': v, '% of related pairs': round(v / tot * 100, 3)}
    for k, v in quad.items()
])
(ROOT / 'outputs/tables').mkdir(parents=True, exist_ok=True)
qdf.to_csv(ROOT / 'outputs/tables/supp_sas_quadrants.csv', index=False)
print("  SAS quadrant counts:")
for k, v in quad.items():
    print(f"    {k}: {v:,} ({v/tot*100:.3f}%)")

# ── Layout ────────────────────────────────────────────────────────────────────
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(DOUBLE, 3.2),
                                  gridspec_kw={'width_ratios': [1.35, 1.0]},
                                  constrained_layout=True)

# ── Panel a: SAS map (hexbin) ─────────────────────────────────────────────────
print("[Panel a] SAS hexbin ...")
xmax = pairs['abs_delta'].max() * 1.02
hb = ax_a.hexbin(pairs['tanimoto'], pairs['abs_delta'],
                 gridsize=55, bins='log', cmap='Greys',
                 mincnt=1, linewidths=0.0, extent=(0.6, 1.0, 0, xmax))

# Shade the activity-cliff quadrant
ax_a.axhspan(DELTA_CUT, xmax, xmin=(TAN_CUT - 0.6) / 0.4, xmax=1.0,
             color=SEM['cliff'], alpha=0.07, zorder=1)
ax_a.axvline(TAN_CUT, color=SEM['cliff'], lw=0.8, ls='--', zorder=3)
ax_a.axhline(DELTA_CUT, color=SEM['cliff'], lw=0.8, ls='--', zorder=3)

# Overlay severe cliff pairs as distinct points
ax_a.scatter(sev['tanimoto'], sev['abs_delta'], s=9,
             c=SEM['cliff'], edgecolors='white', linewidths=0.25,
             alpha=0.9, zorder=4, rasterized=True,
             label=f'Severe cliffs (n={len(sev)})')

# Quadrant annotations
ax_a.text(0.985, 0.30, 'smooth SAR\n(continuous)', fontsize=5.5,
          color=C['gray_dark'], ha='right', va='bottom', style='italic',
          bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.75))
ax_a.text(0.985, xmax * 0.94, 'ACTIVITY\nCLIFFS', fontsize=6.5,
          color=SEM['cliff'], ha='right', va='top', fontweight='bold')

ax_a.set_xlabel('Pairwise Tanimoto similarity (ECFP4)')
ax_a.set_ylabel('|ΔpIC50|')
ax_a.set_xlim(0.6, 1.005)
ax_a.set_ylim(0, xmax)
ax_a.legend(fontsize=5.5, loc='upper left', framealpha=0.0)

cb = fig.colorbar(hb, ax=ax_a, shrink=0.72, aspect=16, pad=0.02)
cb.set_label('Pairs per bin (log)', fontsize=6)
cb.ax.tick_params(labelsize=5.5)
panel_label(ax_a, 'a', x=-0.14, y=1.03)

# ── Panel b: |ΔpIC50| distribution by similarity bin ──────────────────────────
print("[Panel b] Delta distribution by similarity bin ...")
bins   = [(0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 1.01)]
labels = ['0.6–0.7', '0.7–0.8', '0.8–0.9', '0.9–1.0']
colors = [C['grey'], C['cyan'], C['blue'], C['navy']]

box_data = []
for lo, hi in bins:
    sub = pairs[(pairs['tanimoto'] >= lo) & (pairs['tanimoto'] < hi)]['abs_delta']
    box_data.append(sub.values)

parts = ax_b.violinplot(box_data, positions=range(len(bins)),
                        showmedians=False, showextrema=False, widths=0.85)
for pc, col in zip(parts['bodies'], colors):
    pc.set_facecolor(col); pc.set_alpha(0.55)
    pc.set_edgecolor(C['black']); pc.set_linewidth(0.4)

# median markers + cliff-rate annotation per bin
for i, (data, (lo, hi)) in enumerate(zip(box_data, bins)):
    med = np.median(data)
    ax_b.hlines(med, i - 0.30, i + 0.30, color=C['black'], lw=1.4, zorder=4)
    cliff_rate = np.mean(data >= DELTA_CUT) * 100
    ax_b.text(i, DELTA_CUT + 0.12, f'{cliff_rate:.2f}%', ha='center', va='bottom',
              fontsize=5.5, color=SEM['cliff'], fontweight='bold', zorder=5,
              bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.8))

ax_b.axhline(DELTA_CUT, color=SEM['cliff'], lw=0.8, ls='--', zorder=2)
ax_b.text(len(bins) - 0.5, DELTA_CUT - 0.18, 'cliff threshold', fontsize=5.5,
          color=SEM['cliff'], ha='right', va='top')

ax_b.set_xticks(range(len(bins)))
ax_b.set_xticklabels(labels, fontsize=6)
ax_b.set_xlabel('Tanimoto similarity bin')
ax_b.set_ylabel('|ΔpIC50|')
ax_b.set_title('% above cliff threshold annotated', fontsize=6, pad=3)
panel_label(ax_b, 'b', x=-0.18, y=1.03)

save_fig(fig, str(OUT / 'fig_s03_sas_map'))
plt.close(fig)
print("Figure S10 complete.")
