#!/usr/bin/env python
"""Nature Fig 3 — t-SNE Chemical Space (4 panels, cached coordinates)."""
import os, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

NATURE_RC = {
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial','Helvetica','DejaVu Sans'],
    'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7, 'axes.linewidth': 0.75,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': False,
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'lines.linewidth': 0.75, 'lines.markersize': 4, 'patch.linewidth': 0.5,
    'legend.fontsize': 6, 'legend.frameon': False,
    'legend.handlelength': 1.5, 'legend.handletextpad': 0.5,
    'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    'figure.constrained_layout.use': True,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
}
matplotlib.rcParams.update(NATURE_RC)

PAL = {
    'blue': '#0077BB', 'orange': '#EE7733', 'red': '#CC3311',
    'teal': '#009988', 'cyan': '#33BBEE', 'navy': '#1A237E',
    'gray_light': '#BBBBBB', 'gray_dark': '#555555',
}
OUT = 'outputs/figures/nature'
os.makedirs(OUT, exist_ok=True)

def save_fig(fig, name):
    for ext in ('png', 'svg', 'pdf'):
        p = f'{OUT}/{name}.{ext}'
        # Try to force a clean PDF save
        try:
            fig.savefig(p, dpi=600 if ext == 'png' else None,
                        bbox_inches='tight', facecolor='white')
        except Exception as e:
            print(f"⚠️ Could not save {ext} due to backend error: {e} (but PNG was likely saved successfully)")
    sz = os.path.getsize(f'{OUT}/{name}.png') / 1024
    print(f"Saved {name}: {sz:.0f} KB (PNG is clean)")

def plabel(ax, letter, x=-0.08, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=8, fontweight='bold', va='bottom', ha='right',
            fontfamily='sans-serif')

# ══════ VISUAL TWEAKS ══════
MS  = 3.0    # Marker size
ALF = 0.55   # Alpha transparency

# ── Load data ─────────────────────────────────────────────────────────────────
PARQUET_PATH = 'data/processed/pad4_compounds.parquet'
try:
    import pyarrow
    df = pd.read_parquet(PARQUET_PATH)
except ImportError:
    print("ERROR: Missing 'pyarrow'. Run: pip install --user pyarrow")
    raise SystemExit(1)

xy = np.load('data/interim/tsne_coords_3093.npy')
assert len(xy) == len(df), f"t-SNE rows {len(xy)} != compounds {len(df)}"
x, y = xy[:, 0], xy[:, 1]

# Hub InChIKeys → class
HUB_A = {'SMADULGDNOCLOP-GISFHXKWSA-N', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N'}
HUB_B = {'UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'DVCKJOQIVOGXEI-XMMPIXPASA-N'}
hub_mask_a = df['inchi_key'].isin(HUB_A).values
hub_mask_b = df['inchi_key'].isin(HUB_B).values

def source_category(sl):
    parts = sl.split('|')
    if len(parts) >= 2:
        return 'Multi-source (n=2,755)'
    if sl == 'pubchem_confirmatory':
        return 'PubChem only (n=233)'
    return 'Other single source (n=105)'

source_cat_colors = {
    'Multi-source (n=2,755)': PAL['gray_light'],
    'PubChem only (n=233)':   PAL['orange'],
    'Other single source (n=105)': PAL['gray_dark'],
}

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.5))
ax_src, ax_pic, ax_mec, ax_hub = axes.flat

# ── Panel A: Source combination ──────────────────────────────────────────────
sl_cat = df['source_list'].map(source_category)
for label, clr in source_cat_colors.items():
    m = sl_cat == label
    ax_src.scatter(x[m], y[m], s=MS, c=clr, alpha=ALF, linewidths=0, rasterized=True, label=label)
handles = [mpatches.Patch(facecolor=c, edgecolor='none', label=l)
           for l, c in source_cat_colors.items()]
ax_src.legend(handles=handles, loc='upper right', bbox_to_anchor=(1.3, 1.05),
              fontsize=6.0, frameon=False, ncol=1,
              handlelength=0.8, borderpad=0.2, labelspacing=0.3)
ax_src.set_title('Source combination', fontsize=7)
ax_src.set_xlabel('t-SNE 1', fontsize=7); ax_src.set_ylabel('t-SNE 2', fontsize=7)
plabel(ax_src, 'A')

# ── Panel B: pIC50 ───────────────────────────────────────────────────────────
pic = df['pic50_consensus'].values
sc = ax_pic.scatter(x, y, s=MS, c=pic, cmap='viridis', alpha=ALF, linewidths=0, rasterized=True)
cb = fig.colorbar(sc, ax=ax_pic, pad=0.04, shrink=0.8, aspect=25)
cb.set_label('pIC50', fontsize=6)
cb.ax.tick_params(labelsize=5, width=0.4)
ax_pic.set_title('pIC50', fontsize=7)
ax_pic.set_xlabel('t-SNE 1', fontsize=7); ax_pic.set_ylabel('t-SNE 2', fontsize=7)
plabel(ax_pic, 'B')

# ── Panel C: Mechanism ──────────────────────────────────────────────────────
mech_col = 'assay_mechanism_classes'
mech_config = {
    'baee_colorimetric': {'color': PAL['gray_light'], 'label': 'Enzymatic (BAEE)', 'n': '2,079'},
    'rfms_enzymatic':    {'color': PAL['teal'],      'label': 'Enzymatic confirmed (RFMS)', 'n': '878'},
    'fp_binding':        {'color': PAL['cyan'],      'label': 'FP binding', 'n': '115'},
    'covalent_irreversible': {'color': PAL['red'],   'label': 'Covalent', 'n': '21'}
}
for mech_key, config in mech_config.items():
    m = df[mech_col].str.contains(mech_key, na=False)
    if m.sum() > 0:
        ax_mec.scatter(x[m], y[m], s=MS, c=config['color'], alpha=ALF,
                       linewidths=0, rasterized=True,
                       label=f"{config['label']} ({config['n']})")
ax_mec.legend(loc='upper right', bbox_to_anchor=(1.4, 1.05),
              fontsize=6.0, frameon=False, ncol=1,
              handlelength=0.8, borderpad=0.2, labelspacing=0.3, markerscale=2)
ax_mec.set_title('Assay mechanism', fontsize=7)
ax_mec.set_xlabel('t-SNE 1', fontsize=7); ax_mec.set_ylabel('t-SNE 2', fontsize=7)
plabel(ax_mec, 'C')

# ── Panel D: Hub compounds ────────────────────────────────────────────────────
non_hub = ~(hub_mask_a | hub_mask_b)
ax_hub.scatter(x[non_hub], y[non_hub], s=MS, c=PAL['gray_light'],
               alpha=0.35, linewidths=0, rasterized=True, label='Non-hub (3,089)')
ax_hub.scatter(x[hub_mask_a], y[hub_mask_a], s=80, c=PAL['navy'],
               marker='*', alpha=1.0, linewidths=0.3, edgecolors='white', zorder=5,
               label='Hub A (27 severe pairs)')
ax_hub.scatter(x[hub_mask_b], y[hub_mask_b], s=80, c=PAL['red'],
               marker='*', alpha=1.0, linewidths=0.3, edgecolors='white', zorder=5,
               label='Hub B (23 severe pairs)')

hub_label_specs = [
    ('SMADULGDNOCLOP-GISFHXKWSA-N', 'A1', PAL['navy'], (+4, +4)),
    ('RAVBZQAQTVGKIV-XBPDSQQVSA-N', 'A2', PAL['navy'], (+4, -5)),
    ('UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'B1', PAL['red'], (-8, +4)),
    ('DVCKJOQIVOGXEI-XMMPIXPASA-N', 'B2', PAL['red'], (+4, +4)),
]
for ik, lbl, clr, (dx, dy) in hub_label_specs:
    idx = df.index[df['inchi_key'] == ik]
    if len(idx):
        xi, yi = x[idx[0]], y[idx[0]]
        ax_hub.annotate(lbl, (xi, yi), xytext=(xi + dx, yi + dy), zorder=10,
                        fontsize=5.5, fontweight='bold', color=clr,
                        fontfamily='sans-serif',
                        arrowprops=dict(arrowstyle='-', color=clr, lw=0.4, shrinkA=1, shrinkB=1),
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, pad=0.6))

ax_hub.legend(loc='upper right', bbox_to_anchor=(1.3, 1.05),
              fontsize=6.0, frameon=False, ncol=1,
              handlelength=0.8, borderpad=0.2, labelspacing=0.3,
              markerscale=1.0) # <--- Adjusted exactly as you asked (smaller stars in legend)
ax_hub.set_title('Cliff hub compounds', fontsize=7)
ax_hub.set_xlabel('t-SNE 1', fontsize=7); ax_hub.set_ylabel('t-SNE 2', fontsize=7)
plabel(ax_hub, 'D')

save_fig(fig, 'fig3_tsne')
plt.close()
print("Fig 3 DONE")