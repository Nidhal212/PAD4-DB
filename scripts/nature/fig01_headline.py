"""
fig01_headline.py — Figure 1: Headline 4-panel figure
Outputs: outputs/figures/nature/fig01_headline.{png,pdf}
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as FancyBbox
from pathlib import Path
from scipy import stats

# ── rcParams ──────────────────────────────────────────────────────────────────
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
    'blue':       '#0077BB',
    'orange':     '#EE7733',
    'teal':       '#009988',
    'cyan':       '#33BBEE',
    'magenta':    '#EE3377',
    'red':        '#CC3311',
    'navy':       '#004488',
    'grey':       '#BBBBBB',
    'dark_grey':  '#555555',
    'light_grey': '#E8E8E8',
    'olive':      '#999933',
    'purple':     '#AA4499',
}

fam_to_color = {
    'azaindole-benzimidazole biaryl amide derivatives': COLORS['navy'],
    'indazole-N-alkylindole biaryl amide derivatives': COLORS['blue'],
    'indazole-azaindole biaryl amide derivatives': COLORS['blue'],
    'indole-benzimidazole biaryl amide derivatives': COLORS['blue'],
    'chalcone-oxindole derivatives': COLORS['orange'],
    'chalcone-bicyclic lactam derivatives': COLORS['orange'],
    'benzimidazolyl-dihydroisoquinolinone derivatives': COLORS['teal'],
    'bis-benzimidazolyl biaryl diamide derivatives': COLORS['cyan'],
    'Other': COLORS['grey'],
}

CANON = {
    'n_compounds': 3093,
    'n_in_severe': 99,
    'n_patent': 233,
    'n_severe': 94,
    'n_scaffolds': 1244,
    'n_in_series': 2224,
    'n_mmp_validated': 85,
    'n_multi_06': 528,
}

HUB_IKS = {
    'A1': 'SMADULGDNOCLOP-GISFHXKWSA-N',
    'A2': 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
    'B1': 'UDCDEKJNAMHBFH-HSZRJFAPSA-N',
    'B2': 'DVCKJOQIVOGXEI-XMMPIXPASA-N',
}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/figures/nature_v2'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FIGURE 1 — HEADLINE (4-panel)")
print("=" * 60)

# ── Load data ──────────────────────────────────────────────────────────────────
print("\n[Load] shared_assets.parquet ...")
df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
assert len(df) == CANON['n_compounds'], f"n_compounds: {len(df)}"
print(f"  {len(df)} compounds ✓")

print("[Load] t-SNE coordinates ...")
tsne_coords = np.load(ROOT / 'data/interim/tsne_coords_3093.npy')
tsne_iks = np.load(ROOT / 'data/interim/tsne_inchikeys_3093.npy', allow_pickle=True)
print(f"  t-SNE shape: {tsne_coords.shape}")

# Align t-SNE to shared_assets
ik_to_idx = {ik: i for i, ik in enumerate(tsne_iks)}
tsne_df = pd.DataFrame({'inchi_key': tsne_iks, 'tx': tsne_coords[:, 0], 'ty': tsne_coords[:, 1]})
df_tsne = df.merge(tsne_df, on='inchi_key', how='inner')
print(f"  Aligned {len(df_tsne)} compounds to t-SNE ✓")

# Load cliffs for severe nodes
cliffs = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
severe = cliffs[cliffs['cliff_tier'] == 'severe']
severe_iks = set(severe['inchi_key_a'].tolist() + severe['inchi_key_b'].tolist())
assert len(severe_iks) == CANON['n_in_severe'], f"n_in_severe: {len(severe_iks)}"
print(f"  Severe cliff nodes: {len(severe_iks)} ✓")

# ── Build figure ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(10, 12), constrained_layout=True)
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

ax_a = fig.add_subplot(gs[0, 0])  # t-SNE
ax_b = fig.add_subplot(gs[0, 1])  # Scaffold cards
ax_c = fig.add_subplot(gs[1, 0])  # pIC50 KDE
ax_d = fig.add_subplot(gs[1, 1])  # Stats boxes

# ── Panel a: t-SNE ─────────────────────────────────────────────────────────────
print("\n[Panel a] t-SNE ...")

# Plot 'Other' first (background)
other = df_tsne[df_tsne['scaffold_family_group'] == 'Other']
ax_a.scatter(other['tx'], other['ty'], s=3, alpha=0.25, c=COLORS['grey'],
             marker='.', rasterized=True, zorder=1)

# Named families (sorted by size, largest last = on top)
named_fams = [f for f in fam_to_color if f != 'Other']
fam_counts = df_tsne[df_tsne['scaffold_family_group'].isin(named_fams)].groupby('scaffold_family_group').size()
fam_order = fam_counts.sort_values().index.tolist()

legend_handles = [mpatches.Patch(color=COLORS['grey'], label='Other', alpha=0.5)]
for fam in fam_order:
    sub = df_tsne[df_tsne['scaffold_family_group'] == fam]
    color = fam_to_color[fam]
    sz = 8 if 'azaindole-benzimidazole' in fam else 6
    alpha = 0.8 if 'azaindole-benzimidazole' in fam else 0.6
    ax_a.scatter(sub['tx'], sub['ty'], s=sz, alpha=alpha, c=color,
                 marker='o', rasterized=True, zorder=2,
                 label=fam[:20] + ('...' if len(fam) > 20 else ''))
    short = fam.replace(' derivatives', '').replace(' biaryl amide', '')[:25]
    legend_handles.append(mpatches.Patch(color=color, label=short, alpha=0.8))

# Patent compounds
patent_df = df_tsne[df_tsne['patent_flag']]
ax_a.scatter(patent_df['tx'], patent_df['ty'], s=20, alpha=0.7, c=COLORS['orange'],
             marker='x', zorder=4, lw=1.0, label=f'Patent-only (n={len(patent_df)})')
legend_handles.append(mpatches.Patch(color=COLORS['orange'], label=f'Patent-only (n={len(patent_df)})'))

# Severe cliff compounds
severe_sub = df_tsne[df_tsne['inchi_key'].isin(severe_iks)]
ax_a.scatter(severe_sub['tx'], severe_sub['ty'], s=40, alpha=0.9, c=COLORS['red'],
             marker='*', zorder=5, edgecolors='white', lw=0.5,
             label=f'Severe cliff (n={len(severe_sub)})')
legend_handles.append(mpatches.Patch(color=COLORS['red'], label=f'Severe cliff (n={len(severe_sub)})'))

ax_a.set_xlabel('t-SNE 1', fontsize=9)
ax_a.set_ylabel('t-SNE 2', fontsize=9)
ax_a.set_xticks([])
ax_a.set_yticks([])
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.legend(handles=legend_handles[:9], ncol=1, fontsize=6.5, loc='upper left',
            framealpha=0.7, edgecolor='none')
ax_a.text(-0.12, 1.02, 'a', transform=ax_a.transAxes, fontsize=11, fontweight='bold',
          va='bottom', ha='left')

# ── Panel b: Scaffold cards ────────────────────────────────────────────────────
print("\n[Panel b] Scaffold structure cards ...")
ax_b.set_xlim(0, 1)
ax_b.set_ylim(0, 1)
ax_b.axis('off')
ax_b.text(-0.08, 1.02, 'b', transform=ax_b.transAxes, fontsize=11, fontweight='bold',
          va='bottom', ha='left')

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import rdMolDescriptors
    from io import BytesIO
    from PIL import Image as PILImage

    # Top 5 families by compound count in dataset
    top_fams = (df_tsne[df_tsne['scaffold_family_group'] != 'Other']
                .groupby('scaffold_family_group').size()
                .sort_values(ascending=False).index[:5].tolist())

    # Get representative scaffold SMILES from scaffold summary
    scaffold_sum = pd.read_csv(ROOT / 'outputs/tables/05_scaffold_summary.csv')
    top50 = pd.read_csv(ROOT / 'outputs/tables/scaffold_top50_review.csv')

    fam_info = {}
    for fam in top_fams:
        rows = top50[top50['scaffold_family'] == fam]
        if len(rows) > 0:
            r = rows.iloc[0]
            n = df_tsne[df_tsne['scaffold_family_group'] == fam].shape[0]
            mean_pic50 = df_tsne[df_tsne['scaffold_family_group'] == fam]['pIC50'].mean()
            fam_info[fam] = {
                'smiles': r['scaffold_smiles'],
                'n': n,
                'mean_pic50': mean_pic50,
                'color': fam_to_color.get(fam, COLORS['grey'])
            }

    n_cards = len(fam_info)
    cols = 1
    rows_n = n_cards
    card_h = 1.0 / (rows_n + 0.5)
    card_w = 0.85

    for i, (fam, info) in enumerate(fam_info.items()):
        y_center = 1.0 - (i + 0.6) * card_h
        x_center = 0.45

        try:
            mol = Chem.MolFromSmiles(info['smiles'])
            if mol:
                img = Draw.MolToImage(mol, size=(160, 120))
                im_ax = OffsetImage(np.array(img), zoom=0.5)
                ab = AnnotationBbox(im_ax, (x_center - 0.15, y_center),
                                     xycoords='axes fraction',
                                     frameon=True,
                                     bboxprops=dict(edgecolor=info['color'], lw=2))
                ax_b.add_artist(ab)
        except Exception as e:
            # Fallback: colored rectangle
            rect = mpatches.FancyBboxPatch((x_center - 0.38, y_center - card_h*0.4),
                                            0.3, card_h*0.75,
                                            boxstyle='round,pad=0.01',
                                            facecolor=COLORS['light_grey'],
                                            edgecolor=info['color'], lw=2,
                                            transform=ax_b.transAxes)
            ax_b.add_patch(rect)

        short_name = fam.replace(' derivatives', '').replace(' biaryl amide', '')
        if len(short_name) > 28:
            short_name = short_name[:28] + '...'
        ax_b.text(x_center + 0.18, y_center + 0.02, short_name,
                  transform=ax_b.transAxes, fontsize=6.5, fontweight='bold',
                  ha='left', va='center', color=info['color'])
        ax_b.text(x_center + 0.18, y_center - 0.03,
                  f"n={info['n']} | pIC50={info['mean_pic50']:.2f}",
                  transform=ax_b.transAxes, fontsize=6.5,
                  ha='left', va='center', color=COLORS['dark_grey'])
        if 'azaindole-benzimidazole' in fam:
            ax_b.text(x_center + 0.18, y_center - 0.07, '* Hub scaffold',
                      transform=ax_b.transAxes, fontsize=6.5, fontweight='bold',
                      ha='left', va='center', color=COLORS['navy'])

    ax_b.text(0.5, 0.01, 'n = compounds per scaffold family\n(grouped Murcko series)',
              transform=ax_b.transAxes, fontsize=5.5, ha='center', va='bottom',
              color=COLORS['dark_grey'], style='italic')

except Exception as e:
    print(f"    Scaffold card error (using placeholders): {e}")
    # Placeholder boxes
    families_short = ['azaindole-benzimidazole\nbiaryl amide', 'indazole-N-alkylindole\nbiaryl amide',
                       'benzimidazolyl-\ndihydroisoquinolinone', 'chalcone-oxindole', 'indole-benzimidazole\nbiaryl amide']
    colors_list = [COLORS['navy'], COLORS['blue'], COLORS['teal'], COLORS['orange'], COLORS['blue']]
    ns = [471, 193, 99, 67, 79]
    for i, (fam, col, n) in enumerate(zip(families_short, colors_list, ns)):
        y = 0.88 - i * 0.19
        rect = mpatches.FancyBboxPatch((0.05, y - 0.07), 0.9, 0.14,
                                        boxstyle='round,pad=0.01',
                                        facecolor=COLORS['light_grey'],
                                        edgecolor=col, lw=2,
                                        transform=ax_b.transAxes)
        ax_b.add_patch(rect)
        ax_b.text(0.5, y, fam, transform=ax_b.transAxes,
                  fontsize=7, ha='center', va='center', color=col, fontweight='bold')
        ax_b.text(0.5, y - 0.04, f"n={n}", transform=ax_b.transAxes,
                  fontsize=6.5, ha='center', va='center', color=COLORS['dark_grey'])
    ax_b.text(0.5, 0.01, 'n = compounds per scaffold family\n(grouped Murcko series)',
              transform=ax_b.transAxes, fontsize=5.5, ha='center', va='bottom',
              color=COLORS['dark_grey'], style='italic')

# ── Panel c: pIC50 KDE by family ──────────────────────────────────────────────
print("\n[Panel c] pIC50 KDE by family ...")

# Top 5 families by compound count
top_fams_c = (df_tsne[df_tsne['scaffold_family_group'] != 'Other']
              .groupby('scaffold_family_group').size()
              .sort_values(ascending=False).index[:5].tolist())

x_range = np.linspace(2.0, 9.0, 300)

# Plot Other as grey background
other_vals = df_tsne[df_tsne['scaffold_family_group'] == 'Other']['pIC50'].dropna()
kde_other = stats.gaussian_kde(other_vals)
y_other = kde_other(x_range)
ax_c.fill_between(x_range, y_other, alpha=0.15, color=COLORS['grey'])
ax_c.plot(x_range, y_other, lw=1.0, color=COLORS['grey'], label='Other')

for fam in top_fams_c:
    vals = df_tsne[df_tsne['scaffold_family_group'] == fam]['pIC50'].dropna()
    if len(vals) < 5:
        continue
    color = fam_to_color.get(fam, COLORS['grey'])
    kde = stats.gaussian_kde(vals)
    y = kde(x_range)
    short = fam.replace(' derivatives', '').replace(' biaryl amide', '')[:22]
    ax_c.fill_between(x_range, y, alpha=0.2, color=color)
    ax_c.plot(x_range, y, lw=1.5, color=color, label=short)

ax_c.set_xlim(2.0, 9.0)
ax_c.set_xlabel('pIC50', fontsize=9)
ax_c.set_ylabel('Density', fontsize=9)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
ax_c.legend(fontsize=6.5, loc='upper left', ncol=1, framealpha=0.7, edgecolor='none')
ax_c.text(-0.12, 1.02, 'c', transform=ax_c.transAxes, fontsize=11, fontweight='bold',
          va='bottom', ha='left')

# ── Panel d: Stats boxes ──────────────────────────────────────────────────────
print("\n[Panel d] Stats boxes ...")
ax_d.set_xlim(0, 1)
ax_d.set_ylim(0, 1)
ax_d.axis('off')
ax_d.text(-0.08, 1.02, 'd', transform=ax_d.transAxes, fontsize=11, fontweight='bold',
          va='bottom', ha='left')

stats_data = [
    ("3,093", "curated\ninhibitors", COLORS['navy']),
    ("1,244", "unique\nscaffolds", COLORS['blue']),
    ("71.8%", "in scaffold\nseries", COLORS['teal']),
    ("94",    "severe\ncliffs", COLORS['red']),
    ("85/94", "MMP-\nvalidated", COLORS['orange']),
    ("528",   "truly\nindependent", COLORS['magenta']),
]

n_cols, n_rows = 3, 2
box_w = 0.30
box_h = 0.38
pad_x = 0.04
pad_y = 0.10

for idx, (num, label, color) in enumerate(stats_data):
    col = idx % n_cols
    row = idx // n_cols
    x0 = col * (box_w + pad_x) + 0.02
    y0 = 1.0 - (row + 1) * (box_h + pad_y) + pad_y * 0.5

    rect = mpatches.FancyBboxPatch(
        (x0, y0), box_w, box_h,
        boxstyle='round,pad=0.02',
        facecolor=color + '1A',  # ~10% alpha
        edgecolor=color, lw=2.0,
        transform=ax_d.transAxes
    )
    ax_d.add_patch(rect)

    cx = x0 + box_w / 2
    cy = y0 + box_h / 2

    ax_d.text(cx, cy + 0.045, num,
              transform=ax_d.transAxes, fontsize=18, fontweight='bold',
              ha='center', va='center', color=color)
    ax_d.text(cx, cy - 0.055, label,
              transform=ax_d.transAxes, fontsize=7.5,
              ha='center', va='center', color=COLORS['dark_grey'],
              multialignment='center')

# ── Save ──────────────────────────────────────────────────────────────────────
for ext in ['png', 'pdf']:
    outpath = OUT / f'fig01_headline.{ext}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"  Saved: {outpath}")

plt.close(fig)
print("\nFigure 1 complete.")
