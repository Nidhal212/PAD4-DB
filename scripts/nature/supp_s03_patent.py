"""
supp_s03_patent.py — S3: Patent Scaffold Analysis (single panel)
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Liberation Sans', 'DejaVu Sans', 'Arial'],
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 3,    'ytick.major.size': 3,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,     'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

COLORS = {
    'blue': '#0077BB', 'orange': '#EE7733', 'navy': '#004488',
    'grey': '#BBBBBB', 'dark_grey': '#555555', 'light_grey': '#E8E8E8',
}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/figures/nature_v2'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("SUPPLEMENTARY S3 — PATENT SCAFFOLD ANALYSIS")
print("=" * 60)

df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
tsne_coords = np.load(ROOT / 'data/interim/tsne_coords_3093.npy')
tsne_iks = np.load(ROOT / 'data/interim/tsne_inchikeys_3093.npy', allow_pickle=True)
tsne_df = pd.DataFrame({'inchi_key': tsne_iks, 'tx': tsne_coords[:, 0], 'ty': tsne_coords[:, 1]})
df_tsne = df.merge(tsne_df, on='inchi_key', how='inner')

n_patent = df['patent_flag'].sum()
print(f"  Patent-only: {n_patent}")

fig, ax = plt.subplots(figsize=(6, 4))

published = df_tsne[~df_tsne['patent_flag']]
patent = df_tsne[df_tsne['patent_flag']]

ax.scatter(published['tx'], published['ty'], s=4, alpha=0.4, c=COLORS['blue'],
            marker='.', rasterized=True, zorder=1, label=f'Published (n={len(published):,})')
ax.scatter(patent['tx'], patent['ty'], s=15, alpha=0.8, c=COLORS['orange'],
            marker='o', rasterized=True, zorder=3, edgecolors='white', lw=0.3,
            label=f'Patent-exclusive (n={len(patent)})')

ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_title('Patent-exclusive compounds occupy peripheral chemical space',
              fontsize=9, fontweight='bold')
ax.legend(fontsize=8, loc='lower left', framealpha=0.8, edgecolor='none')

footnote = ('233 patent-exclusive compounds  ·  103 exclusive scaffolds  ·  '
            '1 severe cliff contribution')
ax.text(0.5, -0.08, footnote, transform=ax.transAxes, fontsize=7,
         ha='center', va='top', color=COLORS['dark_grey'],
         style='italic')

for ext in ['png', 'pdf']:
    outpath = OUT / f'supp_s03_patent.{ext}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"  Saved: {outpath}")

plt.close(fig)
print("S3 complete.")
