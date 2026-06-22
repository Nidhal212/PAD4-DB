"""
fig04_scaffold.py — Figure 4: Scaffold Landscape (2×2, DOUBLE width)

Panel a: Top-15 scaffold series (horizontal bars, coloured by mean pIC50)
Panel b: Scaffold concentration — Lorenz curve (Gini=0.532)
Panel c: SAR ruggedness — series size vs intra-scaffold pIC50 spread (std);
         directly visualises "scaffold-dependent SAR ruggedness"
Panel d: Cliff density — severe cliff pairs / possible pairs per scaffold (n>=4)

Outputs: publication/figures/main/fig04_scaffold.{png,pdf}
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
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib.lines import Line2D
from pathlib import Path

set_style()

CANON = {'n_compounds': 3093, 'n_scaffolds': 1244, 'n_series': 375, 'largest': 174}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT  = ROOT / 'publication/figures/main'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FIGURE 4 — SCAFFOLD LANDSCAPE (2x2)")
print("=" * 60)

assets = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
scaff  = pd.read_csv(ROOT / 'outputs/tables/05_scaffold_summary.csv')
scaff  = scaff.sort_values('n_compounds', ascending=False).reset_index(drop=True)
scaff['scaffold_rank'] = scaff.index + 1
assert len(scaff) == CANON['n_scaffolds'], f"scaffolds: {len(scaff)}"

# ── Cliff density per scaffold ────────────────────────────────────────────────
ac  = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
sev = ac[ac['cliff_tier'] == 'severe'].copy()
ik_to_scaffold = dict(zip(assets['inchi_key'], assets['murcko_smiles']))
sev['scaffold_a'] = sev['inchi_key_a'].map(ik_to_scaffold)
sev['scaffold_b'] = sev['inchi_key_b'].map(ik_to_scaffold)
same = sev[sev['scaffold_a'] == sev['scaffold_b']]
cliff_counts = same.groupby('scaffold_a').size().rename('n_cliff_pairs')

scaff_series = scaff[scaff['n_compounds'] >= 2].copy()
scaff_series = scaff_series.merge(
    cliff_counts.reset_index().rename(columns={'scaffold_a': 'scaffold_smiles'}),
    on='scaffold_smiles', how='left')
scaff_series['n_cliff_pairs'] = scaff_series['n_cliff_pairs'].fillna(0).astype(int)
scaff_series['n_possible_pairs'] = scaff_series['n_compounds'] * (scaff_series['n_compounds'] - 1) / 2
scaff_series['cliff_density'] = scaff_series['n_cliff_pairs'] / scaff_series['n_possible_pairs']

# ── Figure layout: 2×2 ────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(DOUBLE, 5.4), constrained_layout=True)
ax_a, ax_b = axes[0]
ax_c, ax_d = axes[1]

pic_norm = Normalize(vmin=4.0, vmax=7.5)
pic_cmap = cm.viridis

# ── Panel a: Top-15 scaffold series, coloured by mean pIC50 ────────────────────
print("[Panel a] Top-15 scaffold series ...")
top = scaff.head(15).copy()
y_pos = np.arange(len(top))[::-1]   # rank 1 at top

bar_colors = [pic_cmap(pic_norm(v)) for v in top['mean_pic50']]
ax_a.barh(y_pos, top['n_compounds'].values, color=bar_colors,
          height=0.74, edgecolor='white', lw=0.4)

# Mark patent-exclusive scaffolds with an orange edge
for i, row in enumerate(top.itertuples()):
    if getattr(row, 'patent_exclusive_scaffold', False):
        ax_a.barh(y_pos[i], row.n_compounds, height=0.74,
                  fill=False, edgecolor=SEM['patent'], lw=1.3, zorder=5)

for yp, row in zip(y_pos, top.itertuples()):
    ax_a.text(row.n_compounds + 2, yp, f'{row.n_compounds}',
              va='center', ha='left', fontsize=5.5, color=C['black'])

ax_a.set_yticks(y_pos)
ax_a.set_yticklabels([f'S{r}' for r in top['scaffold_rank']], fontsize=5.5)
ax_a.set_xlabel('Compounds in series')
ax_a.set_xlim(0, top['n_compounds'].max() * 1.16)
ax_a.set_ylabel('Scaffold rank')

sm = cm.ScalarMappable(norm=pic_norm, cmap=pic_cmap); sm.set_array([])
cb = fig.colorbar(sm, ax=ax_a, shrink=0.72, aspect=16, pad=0.02)
cb.set_label('Mean pIC50', fontsize=6); cb.ax.tick_params(labelsize=5.5)
ax_a.scatter([], [], marker='s', facecolor='none', edgecolor=SEM['patent'],
             lw=1.3, s=30, label='Patent-exclusive')
ax_a.legend(fontsize=5.5, loc='lower right', framealpha=0.0,
            handletextpad=0.4, borderpad=0.2)
panel_label(ax_a, 'a', x=-0.16, y=1.03)

# ── Panel b: Lorenz curve (scaffold concentration) ────────────────────────────
print("[Panel b] Lorenz curve ...")
vals = np.sort(scaff['n_compounds'].values)
cum  = np.cumsum(vals)
lx = np.concatenate([[0], np.arange(1, len(vals) + 1) / len(vals)])
ly = np.concatenate([[0], cum / cum[-1]])
gini = 1 - 2 * np.trapz(ly, lx)
print(f"  Gini: {gini:.3f}")

ax_b.plot(lx, ly, color=SEM['sar'], lw=1.6, label='Lorenz curve', zorder=3)
ax_b.plot([0, 1], [0, 1], color=C['gray_dark'], lw=0.8, ls='--', label='Perfect equality')
ax_b.fill_between(lx, ly, lx, alpha=0.16, color=SEM['sar'], zorder=2)
ax_b.text(0.97, 0.30, f'Gini = {gini:.3f}', transform=ax_b.transAxes,
          fontsize=7.5, fontweight='bold', color=C['black'], ha='right')
ax_b.text(0.97, 0.21, '71.9% of compounds in\nseries of ≥2 members',
          transform=ax_b.transAxes, fontsize=5.5, color=C['gray_dark'], ha='right')
ax_b.set_xlabel('Cumulative fraction of scaffolds')
ax_b.set_ylabel('Cumulative fraction of compounds')
ax_b.set_xlim(0, 1); ax_b.set_ylim(0, 1)
ax_b.set_aspect('equal')
ax_b.legend(fontsize=5.5, loc='upper left', framealpha=0.0)
panel_label(ax_b, 'b', x=-0.16, y=1.03)

# ── Panel c: SAR ruggedness — series size vs intra-scaffold pIC50 spread ───────
print("[Panel c] SAR ruggedness ...")
rug = scaff_series.dropna(subset=['std_pic50']).copy()
patent_mask = rug['contains_patent_exclusive'] == True

# size encodes mean potency; colour encodes patent status
ax_c.scatter(rug['n_compounds'][~patent_mask], rug['std_pic50'][~patent_mask],
             s=14, alpha=0.55, c=C['blue'], edgecolors='none',
             rasterized=True, label='Non-patent')
ax_c.scatter(rug['n_compounds'][patent_mask], rug['std_pic50'][patent_mask],
             s=16, alpha=0.75, c=SEM['patent'], edgecolors='none',
             rasterized=True, label='Patent-exclusive')

# Median ruggedness reference line
med_std = rug['std_pic50'].median()
ax_c.axhline(med_std, color=C['gray_dark'], lw=0.7, ls=':', zorder=1)
ax_c.text(rug['n_compounds'].max() * 0.62, med_std + 0.03,
          f'median σ = {med_std:.2f}', fontsize=5.5, color=C['gray_dark'])

# annotate the largest series
r1 = rug[rug['scaffold_rank'] == 1].iloc[0]
ax_c.scatter([r1['n_compounds']], [r1['std_pic50']], s=46, c=C['navy'],
             edgecolors='white', lw=0.6, zorder=6)
ax_c.annotate(f'S1 (n=174, σ={r1["std_pic50"]:.2f})',
              xy=(r1['n_compounds'], r1['std_pic50']),
              xytext=(r1['n_compounds'] - 96, r1['std_pic50'] + 0.34),
              fontsize=5.5, color=C['navy'],
              arrowprops=dict(arrowstyle='->', lw=0.5, color=C['navy']))

ax_c.set_xscale('log')
ax_c.set_xlabel('Scaffold series size (log scale)')
ax_c.set_ylabel('Intra-scaffold pIC50 spread (σ)')
ax_c.set_title('SAR ruggedness', fontsize=6.5, pad=3)
ax_c.legend(fontsize=5.5, loc='upper right', framealpha=0.0)
panel_label(ax_c, 'c', x=-0.16, y=1.03)

# ── Panel d: Cliff density per scaffold (series n>=4) ──────────────────────────
print("[Panel d] Cliff density (n>=4) ...")
cd = scaff_series[scaff_series['n_compounds'] >= 4].copy()
pm = cd['contains_patent_exclusive'] == True
print(f"  n>=4 scaffolds: {len(cd)} | with cliff pairs: {(cd.n_cliff_pairs>0).sum()} | max density {cd.cliff_density.max():.3f}")

ax_d.scatter(cd['n_compounds'][~pm], cd['cliff_density'][~pm],
             s=16, alpha=0.6, c=C['blue'], edgecolors='none',
             rasterized=True, label='Non-patent')
ax_d.scatter(cd['n_compounds'][pm], cd['cliff_density'][pm],
             s=18, alpha=0.8, c=SEM['patent'], edgecolors='none',
             rasterized=True, label='Patent-exclusive')

# label the scaffold with most cliff pairs
top_cliff = cd.sort_values('n_cliff_pairs', ascending=False).iloc[0]
ax_d.annotate(f'{int(top_cliff.n_cliff_pairs)} cliff pairs\n(n={int(top_cliff.n_compounds)})',
              xy=(top_cliff['n_compounds'], top_cliff['cliff_density']),
              xytext=(top_cliff['n_compounds'] + 12, top_cliff['cliff_density'] - 0.02),
              fontsize=5.5, color=C['black'],
              arrowprops=dict(arrowstyle='->', lw=0.5, color=C['gray_dark']))

ax_d.set_xlabel('Scaffold series size')
ax_d.set_ylabel('Severe cliff density')
ax_d.set_title('Cliff pairs / possible pairs (series ≥4)', fontsize=6.5, pad=3)
ax_d.set_ylim(-0.012, cd['cliff_density'].max() * 1.18)
ax_d.legend(fontsize=5.5, loc='upper right', framealpha=0.0)
panel_label(ax_d, 'd', x=-0.16, y=1.03)

save_fig(fig, str(OUT / 'fig04_scaffold'))
plt.close(fig)
print("Figure 4 complete.")
