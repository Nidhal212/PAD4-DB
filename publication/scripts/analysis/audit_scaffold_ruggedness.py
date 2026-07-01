"""
audit_scaffold_ruggedness.py — Phase 2 + Phase 7.1

Per-scaffold-series SAR ruggedness metrics, rugged/smooth ranking, scaffold
diversity metrics, and three publication-quality figures.

Outputs:
  outputs/tables/scaffold_ruggedness_table.csv
  outputs/tables/scaffold_diversity_metrics.csv
  outputs/audit/top20_rugged_scaffolds.csv
  outputs/audit/top20_smooth_scaffolds.csv
  publication/figures/supplementary/fig_s04_ruggedness_panels.{png,pdf}
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, panel_label, save_fig, SEM, C, DOUBLE
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from pathlib import Path

set_style()
ROOT = Path('/home/nidhal/PAD4-db_V2')
(ROOT / 'outputs/audit').mkdir(parents=True, exist_ok=True)

print("=" * 64)
print("PHASE 2 + 7.1 — SCAFFOLD RUGGEDNESS & DIVERSITY")
print("=" * 64)

assets = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
pairs  = pd.read_parquet(ROOT / 'data/processed/activity_pairs_with_sali.parquet')
ac     = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
sev    = ac[ac.cliff_tier == 'severe'].copy()

ik2s   = dict(zip(assets.inchi_key, assets.murcko_smiles))
ik2pic = dict(zip(assets.inchi_key, assets.pIC50))

# within-scaffold severe cliff counts
sev['sa'] = sev.inchi_key_a.map(ik2s); sev['sb'] = sev.inchi_key_b.map(ik2s)
within_cliffs = sev[sev.sa == sev.sb].groupby('sa').size()

# within-scaffold SALI (from full related-pairs set, both endpoints same scaffold)
pairs['sa'] = pairs.inchi_key_a.map(ik2s); pairs['sb'] = pairs.inchi_key_b.map(ik2s)
within_pairs = pairs[pairs.sa == pairs.sb]
sali_by_scaf = within_pairs.groupby('sa')['sali'].agg(['mean', 'max'])

# ── Per-scaffold-series metrics ───────────────────────────────────────────────
rows = []
for scaf, grp in assets.groupby('murcko_smiles'):
    n = len(grp)
    if n < 2:
        continue
    pic = grp.pIC50.dropna().values
    poss = n * (n - 1) / 2
    ncliff = int(within_cliffs.get(scaf, 0))
    dens = ncliff / poss if poss else 0.0
    rows.append({
        'scaffold_smiles': scaf,
        'series_size': n,
        'n_severe_cliffs': ncliff,
        'n_possible_pairs': int(poss),
        'cliff_density': round(dens, 5),
        'smoothness_score': round(1 - dens, 5),
        'mean_pic50': round(float(np.mean(pic)), 3),
        'median_pic50': round(float(np.median(pic)), 3),
        'potency_range': round(float(np.ptp(pic)), 3),
        'sd_pic50': round(float(np.std(pic, ddof=1)), 3) if n > 1 else 0.0,
        'mean_sali': round(float(sali_by_scaf['mean'].get(scaf, 0.0)), 3),
        'max_sali': round(float(sali_by_scaf['max'].get(scaf, 0.0)), 3),
        'contains_patent_exclusive': bool(grp.patent_flag.any()),
    })
rt = pd.DataFrame(rows).sort_values('series_size', ascending=False).reset_index(drop=True)
rt.insert(0, 'scaffold_rank', rt.index + 1)
rt.to_csv(ROOT / 'outputs/tables/scaffold_ruggedness_table.csv', index=False)
print(f"  scaffold series (n>=2): {len(rt)}")
print(f"  series with >=1 cliff: {(rt.n_severe_cliffs > 0).sum()}")

# ── Rugged / smooth ranking ───────────────────────────────────────────────────
# Rugged: rank by cliff_density then mean_sali then max_sali (require some signal)
rugged = rt.sort_values(['cliff_density', 'mean_sali', 'max_sali'], ascending=False).head(20)
# Smooth: large enough to be meaningful (n>=5), lowest cliff density & SALI
smooth = rt[rt.series_size >= 5].sort_values(['cliff_density', 'mean_sali']).head(20)
keep = ['scaffold_rank', 'series_size', 'n_severe_cliffs', 'cliff_density', 'smoothness_score',
        'mean_pic50', 'potency_range', 'mean_sali', 'max_sali', 'contains_patent_exclusive']
rugged[keep].to_csv(ROOT / 'outputs/audit/top20_rugged_scaffolds.csv', index=False)
smooth[keep].to_csv(ROOT / 'outputs/audit/top20_smooth_scaffolds.csv', index=False)
print(f"  top rugged max cliff_density={rugged.cliff_density.max():.3f}; "
      f"smooth (n>=5) all cliff_density=0: {(smooth.cliff_density == 0).all()}")

# ── Scaffold diversity metrics (Phase 7.1) ────────────────────────────────────
sizes = assets.groupby('murcko_smiles').size().values
N = int(sizes.sum()); S = len(sizes)
def gini(x):
    x = np.sort(np.asarray(x, float)); n = len(x); c = np.cumsum(x)
    return (n + 1 - 2 * np.sum(c) / c[-1]) / n
p = sizes / N
shannon = float(-np.sum(p * np.log(p)))
top10 = float(np.sort(sizes)[::-1][:10].sum() / N * 100)
div = pd.DataFrame([
    {'metric': 'n_compounds', 'value': N},
    {'metric': 'n_unique_scaffolds', 'value': S},
    {'metric': 'scaffold_to_compound_ratio', 'value': round(S / N, 4)},
    {'metric': 'singleton_scaffold_fraction', 'value': round(float((sizes == 1).mean()), 4)},
    {'metric': 'n_series_ge2', 'value': int((sizes >= 2).sum())},
    {'metric': 'top10_scaffold_coverage_pct', 'value': round(top10, 2)},
    {'metric': 'gini_coefficient', 'value': round(gini(sizes), 4)},
    {'metric': 'shannon_entropy', 'value': round(shannon, 4)},
    {'metric': 'normalized_shannon_evenness', 'value': round(shannon / np.log(S), 4)},
    {'metric': 'compounds_in_series_pct', 'value': round(float(sizes[sizes >= 2].sum() / N * 100), 2)},
])
div.to_csv(ROOT / 'outputs/tables/scaffold_diversity_metrics.csv', index=False)
print("  diversity:", {r.metric: r.value for r in div.itertuples()})

# ── Figures A/B/C ─────────────────────────────────────────────────────────────
fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(DOUBLE, 2.7), constrained_layout=True)
pnorm = Normalize(vmin=4.0, vmax=7.5)

# A: series size vs cliff density, colour=mean pIC50
sc = axA.scatter(rt.series_size, rt.cliff_density, c=rt.mean_pic50, cmap='viridis',
                 norm=pnorm, s=16, alpha=0.75, edgecolors='none', rasterized=True)
axA.set_xscale('log'); axA.set_xlabel('Series size (log)'); axA.set_ylabel('Cliff density')
axA.set_title('Size vs ruggedness', fontsize=6.5, pad=3)
cb = fig.colorbar(sc, ax=axA, shrink=0.7, aspect=14, pad=0.02); cb.set_label('mean pIC50', fontsize=6)
cb.ax.tick_params(labelsize=5.5)
panel_label(axA, 'a', x=-0.18, y=1.04)

# B: series size vs potency range, colour=cliff density
dnorm = Normalize(vmin=0, vmax=max(rt.cliff_density.max(), 1e-6))
sc2 = axB.scatter(rt.series_size, rt.potency_range, c=rt.cliff_density, cmap='magma_r',
                  norm=dnorm, s=16, alpha=0.8, edgecolors='#4d4d4d', linewidths=0.3,
                  rasterized=True)   # thin dark edge: cliff_density=0 maps to near-white magma_r fill, invisible on white bg without it
axB.set_xscale('log'); axB.set_xlabel('Series size (log)'); axB.set_ylabel('Intra-scaffold pIC50 range')
axB.set_title('Size vs potency range', fontsize=6.5, pad=3)
cb2 = fig.colorbar(sc2, ax=axB, shrink=0.7, aspect=14, pad=0.02); cb2.set_label('cliff density', fontsize=6)
cb2.ax.tick_params(labelsize=5.5)
panel_label(axB, 'b', x=-0.18, y=1.04)

# C: distribution of cliff density across series
axC.hist(rt.cliff_density, bins=np.linspace(0, rt.cliff_density.max() * 1.02, 25),
         color=SEM['cliff'], alpha=0.8, edgecolor='white', lw=0.3)
axC.set_yscale('log')
axC.set_xlabel('Cliff density'); axC.set_ylabel('Scaffold series (log count)')
axC.set_title(f'{(rt.cliff_density==0).mean()*100:.0f}% of series are perfectly smooth',
              fontsize=6.5, pad=3)
panel_label(axC, 'c', x=-0.18, y=1.04)

save_fig(fig, str(ROOT / 'publication/figures/supplementary/fig_s04_ruggedness_panels'))
plt.close(fig)
print("DONE")
