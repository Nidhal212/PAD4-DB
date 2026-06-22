"""
fig01_headline.py — Figure 1: Layered UMAP / t-SNE landscape
ONEHALF width, single panel.

Layers (bottom to top):
  1. Grey background scatter (all 3,093 compounds, alpha=0.4)
  2. Contour density underlay (grey, alpha=0.3)
  3. Severe cliff compounds coloured by pIC50 (viridis)
  4. Hub A (navy star ★) + Hub B (red diamond ◆) with white borders

Outputs: publication/figures/main/fig01_headline.{png,pdf}
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, panel_label, save_fig, SEM, C, ONEHALF

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.colors import Normalize
from matplotlib import cm
from pathlib import Path
from scipy.stats import gaussian_kde

set_style()

CANON = {'n_compounds': 3093, 'n_in_severe': 99, 'n_severe': 94}

HUB_IKS = {
    'A1': 'SMADULGDNOCLOP-GISFHXKWSA-N',
    'A2': 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
    'B1': 'UDCDEKJNAMHBFH-HSZRJFAPSA-N',
    'B2': 'DVCKJOQIVOGXEI-XMMPIXPASA-N',
}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT  = ROOT / 'publication/figures/main'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FIGURE 1 — LAYERED UMAP/t-SNE LANDSCAPE")
print("=" * 60)

# ── Load data ──────────────────────────────────────────────────────────────────
df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
assert len(df) == CANON['n_compounds'], f"n={len(df)}"

tsne_coords = np.load(ROOT / 'data/interim/tsne_coords_3093.npy')
tsne_iks    = np.load(ROOT / 'data/interim/tsne_inchikeys_3093.npy', allow_pickle=True)
tsne_df     = pd.DataFrame({'inchi_key': tsne_iks,
                             'tx': tsne_coords[:, 0],
                             'ty': tsne_coords[:, 1]})
df_t = df.merge(tsne_df, on='inchi_key', how='inner')
print(f"  Aligned {len(df_t)} / {len(df)} compounds")

cliffs   = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
sev      = cliffs[cliffs['cliff_tier'] == 'severe']
sev_iks  = set(sev['inchi_key_a'].tolist() + sev['inchi_key_b'].tolist())
assert len(sev_iks) == CANON['n_in_severe'], f"sev n={len(sev_iks)}"

hub_all  = set(HUB_IKS.values())
hub_a    = {HUB_IKS['A1'], HUB_IKS['A2']}
hub_b    = {HUB_IKS['B1'], HUB_IKS['B2']}
cliff_nonhub = sev_iks - hub_all

df_bg    = df_t[~df_t['inchi_key'].isin(sev_iks)]
df_cliff = df_t[df_t['inchi_key'].isin(cliff_nonhub)]
df_hA    = df_t[df_t['inchi_key'].isin(hub_a)]
df_hB    = df_t[df_t['inchi_key'].isin(hub_b)]

# ── Build figure ──────────────────────────────────────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(ONEHALF, ONEHALF * 0.95), constrained_layout=True)

# Layer 1 — smooth filled density underlay (soft, single-tone; no busy line contours)
xy  = np.vstack([df_t['tx'].values, df_t['ty'].values])
kde = gaussian_kde(xy, bw_method=0.20)
xmin, xmax = df_t['tx'].min() - 3, df_t['tx'].max() + 3
ymin, ymax = df_t['ty'].min() - 3, df_t['ty'].max() + 3
xg, yg     = np.mgrid[xmin:xmax:200j, ymin:ymax:200j]
zg         = kde(np.vstack([xg.ravel(), yg.ravel()])).reshape(xg.shape)
# Filled contour, very light grey, few levels — recedes to pure context so cliffs/hubs pop
ax.contourf(xg, yg, zg, levels=8, cmap='Greys', alpha=0.10, zorder=1)

# Layer 2 — faint grey background points (all non-cliff compounds)
ax.scatter(df_bg['tx'], df_bg['ty'],
           s=1.6, alpha=0.16, c=SEM['background'],
           marker='.', rasterized=True, zorder=2)

# Layer 3 — severe cliff compounds coloured by pIC50 (viridis)
cmap     = cm.viridis
vmin, vmax = 2.0, 8.52
norm     = Normalize(vmin=vmin, vmax=vmax)
node_norm = norm   # alias used in legend

if len(df_cliff) > 0:
    sc = ax.scatter(df_cliff['tx'], df_cliff['ty'],
                    s=28, alpha=0.95,
                    c=df_cliff['pIC50'].values, cmap=cmap, norm=norm,
                    marker='o', edgecolors='white', linewidths=0.5,
                    zorder=3, rasterized=True)

# Small colorbar (pIC50)
cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                    shrink=0.38, aspect=14, pad=0.02,
                    orientation='vertical')
cbar.set_label('pIC50', fontsize=6)
cbar.ax.tick_params(labelsize=5.5)

# Layer 4 — Hub A (navy star) + Hub B (red diamond) with white borders
if len(df_hA) > 0:
    ax.scatter(df_hA['tx'], df_hA['ty'],
               s=180, c=SEM['classA'], marker='*',
               edgecolors='white', linewidths=1.1, zorder=5)
if len(df_hB) > 0:
    ax.scatter(df_hB['tx'], df_hB['ty'],
               s=70, c=SEM['classB'], marker='D',
               edgecolors='white', linewidths=1.0, zorder=5)

# ── Legend ───────────────────────────────────────────────────────────────────
# Use a viridis midpoint color (~pIC50 5.25) as representative cliff color
cliff_proxy_color = cm.viridis(node_norm(5.25))
legend_handles = [
    mpatches.Patch(color=SEM['background'], alpha=0.55, label=f'All compounds (n={len(df_t):,})'),
    mlines.Line2D([], [], color=SEM['classA'], marker='*', linestyle='None',
                  markersize=7, label='Hub A — series floor (n=2)'),
    mlines.Line2D([], [], color=SEM['classB'], marker='D', linestyle='None',
                  markersize=5, label='Hub B — singleton attractor (n=2)'),
    mlines.Line2D([], [], color='none', marker='o',
                  markerfacecolor=cliff_proxy_color,
                  markeredgecolor='white', markersize=5, markeredgewidth=0.4,
                  label=f'Severe cliff compound (n={len(sev_iks)}, coloured by pIC50)'),
]
ax.legend(handles=legend_handles, fontsize=5.5, loc='lower left',
          framealpha=0.88, edgecolor='none',
          handletextpad=0.4, labelspacing=0.35)

# Tighten limits to the data envelope (small margin) — removes dead whitespace
mx = (df_t['tx'].max() - df_t['tx'].min()) * 0.04
my = (df_t['ty'].max() - df_t['ty'].min()) * 0.04
ax.set_xlim(df_t['tx'].min() - mx, df_t['tx'].max() + mx)
ax.set_ylim(df_t['ty'].min() - my, df_t['ty'].max() + my)

ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_xlabel('t-SNE 1', labelpad=2)
ax.set_ylabel('t-SNE 2', labelpad=2)

# ── Save ──────────────────────────────────────────────────────────────────────
save_fig(fig, str(OUT / 'fig01_headline'))
plt.close(fig)
print("Figure 1 complete.")
