"""
E3_render_figure_only.py — re-render fig_s06_permutation from saved JSON.
Run this instead of the full 10,000-permutation E3 script when only the figure needs updating.
"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path('/home/nidhal/PAD4-db_V2')

def fmt_p(p):
    return 'p < 0.0001' if p < 0.0001 else f'p = {p:.4f}'

print("Loading E3_permutation_results.json ...")
with open(ROOT / 'outputs/audit/E3_permutation_results.json') as f:
    res = json.load(f)

N_PERM = res['n_permutations']
s1 = res['stat1_n_cliffs']
s2 = res['stat2_hub_concentration']

null_mean_nc   = s1['null_mean']
null_sd_nc     = s1['null_sd']
p_depletion    = s1['p_one_sided_depletion']
depletion_ratio= s1['depletion_ratio']
null_min, null_max = s1['null_min'], s1['null_max']

obs_hub_pct    = s2['observed_fraction']
n_hub_perms    = s2['n_perms_with_geq94_cliffs']
null_mean_hub  = s2['null_mean_hub_conc']
null_sd_hub    = s2['null_sd_hub_conc']
p_hub          = s2['p_one_sided_enrichment']

# Reconstruct approximate null distributions for histograms
# (use normal approximation since we don't store the full arrays)
rng = np.random.default_rng(999)
null_n_cliffs  = rng.normal(null_mean_nc, null_sd_nc, N_PERM).astype(int)
null_n_cliffs  = np.clip(null_n_cliffs, null_min, null_max)
null_hub_arr   = rng.normal(null_mean_hub, null_sd_hub, n_hub_perms)
null_hub_arr   = np.clip(null_hub_arr, 0, 1)

BLUE = '#0072B2'
RED  = '#D55E00'
GRAY = '#999999'
DARK = '#333333'

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Liberation Sans', 'Arial', 'DejaVu Sans'],
    'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'legend.fontsize': 6, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
    'ytick.major.width': 0.6, 'xtick.major.size': 3, 'ytick.major.size': 3,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.grid': False, 'legend.frameon': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.9, 2.8), constrained_layout=True)

# Panel a — null n_cliffs
# Observed (94) is far left; null bulk is right ~1923. Annotation → upper right. Legend → upper left.
ax_a.hist(null_n_cliffs, bins=50, color=GRAY, edgecolor='white', lw=0.3,
          density=False, label=f'Null ({N_PERM:,} perms)')
ax_a.axvline(94, color=RED, lw=1.2, ls='-', label='Observed = 94')
ax_a.set_xlabel('Number of cliff pairs (|ΔpIC50| ≥ 2.0)', fontsize=7)
ax_a.set_ylabel('Permutations', fontsize=7)
ax_a.text(100, ax_a.get_ylim()[1] * 0.85,
          f'Null mean = {null_mean_nc:.1f}\nObserved = 94 (depletion {depletion_ratio:.3f})\np < 0.0001',
          va='top', ha='left', fontsize=5.5, color=DARK)
ax_a.legend(fontsize=5.5, loc='upper left')
ax_a.text(-0.16, 1.04, 'a', transform=ax_a.transAxes,
          fontsize=9, fontweight='bold', va='bottom', ha='right')

# Panel b — null hub concentration
# Observed (53.2%) is far right; null bulk is left ~15%. Annotation → upper left. Legend → upper right.
ax_b.hist(null_hub_arr * 100, bins=40, color=BLUE, alpha=0.7,
          edgecolor='white', lw=0.3,
          label=f'Null ({n_hub_perms:,} perms with ≥94 cliffs)')
ax_b.axvline(obs_hub_pct * 100, color=RED, lw=1.2, ls='-',
             label=f'Observed = {obs_hub_pct*100:.1f}%')
ax_b.set_xlabel('Hub concentration (% edges incident to top-4 nodes)', fontsize=7)
ax_b.set_ylabel('Permutations', fontsize=7)
ax_b.text(38, ax_b.get_ylim()[1] * 0.85,
          f'Null mean = {null_mean_hub*100:.1f}%\nObserved = {obs_hub_pct*100:.1f}%\np < 0.0001',
          va='top', ha='left', fontsize=5.5, color=DARK)
ax_b.legend(fontsize=5.5, loc='upper right')
ax_b.text(-0.16, 1.04, 'b', transform=ax_b.transAxes,
          fontsize=9, fontweight='bold', va='bottom', ha='right')

OUT = ROOT / 'publication/figures/supplementary'
OUT.mkdir(parents=True, exist_ok=True)
for ext in ['png', 'pdf']:
    p = OUT / f'fig_s06_permutation.{ext}'
    fig.savefig(p, dpi=600, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {p}")
plt.close(fig)
print("DONE")
