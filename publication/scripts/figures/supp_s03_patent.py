"""
supp_s03_patent.py — S3: Patent Compound t-SNE (single panel, ONEHALF width)
Outputs: publication/figures/supplementary/fig_s03_patent.{png,pdf}
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, save_fig, SEM, C, SINGLE, ONEHALF

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

set_style()

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'publication/figures/supplementary'
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

fig, ax = plt.subplots(figsize=(ONEHALF, 3.8))

published = df_tsne[~df_tsne['patent_flag']]
patent = df_tsne[df_tsne['patent_flag']]

ax.scatter(published['tx'], published['ty'], s=2, alpha=0.35, c=SEM['published'],
            marker='.', rasterized=True, zorder=1, label=f'Published (n={len(published):,})')
ax.scatter(patent['tx'], patent['ty'], s=10, alpha=0.75, c=SEM['patent'],
            marker='o', rasterized=True, zorder=3, edgecolors='white', lw=0.2,
            label=f'Patent-exclusive (n={len(patent)})')

ax.set_xticks([])
ax.set_yticks([])
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_xlabel('t-SNE 1')
ax.set_ylabel('t-SNE 2')
ax.legend(fontsize=6, loc='lower left', framealpha=0.8, edgecolor='none')

ax.set_title('n=233 patent-exclusive · 103 unique scaffolds · 1 severe cliff',
              fontsize=6, pad=4)

fig.tight_layout(pad=0.5)
save_fig(fig, str(OUT / 'fig_s03_patent'))
plt.close(fig)
print("S3 complete.")
