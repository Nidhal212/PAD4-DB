#!/usr/bin/env python
"""Nature Fig 1 — Pipeline Workflow Diagram."""
import os, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

NATURE_RC = {
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial','Helvetica','DejaVu Sans'],
    'font.size': 7, 'axes.labelsize': 7, 'axes.linewidth': 0.75,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': False,
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'lines.linewidth': 0.75, 'patch.linewidth': 0.5,
    'legend.fontsize': 6, 'legend.frameon': False,
    'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    'figure.constrained_layout.use': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
}
matplotlib.rcParams.update(NATURE_RC)

PAL = {
    'blue': '#0077BB', 'orange': '#EE7733', 'red': '#CC3311',
    'teal': '#009988', 'navy': '#1A237E',
    'gray_light': '#BBBBBB', 'gray_dark': '#555555',
}
OUT = 'outputs/figures/nature'
os.makedirs(OUT, exist_ok=True)

def save_fig(fig, name):
    for ext in ('png', 'svg', 'pdf'):
        p = f'{OUT}/{name}.{ext}'
        fig.savefig(p, dpi=600 if ext == 'png' else None,
                    bbox_inches='tight', facecolor='white')
    sz = os.path.getsize(f'{OUT}/{name}.png') / 1024
    print(f"Saved {name}: {sz:.0f} KB")

# ── Build figure ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(7.2, 9))
ax  = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_axis_off()

def box(x, y, w, h, text, fc, tc='white', fs=7.0, bold=False, lw=0.6, ec='#DDDDDD'):
    patch = FancyBboxPatch(
        (x - w/2, y - h/2), w, h,
        boxstyle='round,pad=0.035',
        facecolor=fc, edgecolor=ec, linewidth=lw,
        transform=ax.transAxes, clip_on=False, zorder=2)
    ax.add_patch(patch)
    ax.text(x, y, text, transform=ax.transAxes,
            ha='center', va='center', fontsize=fs,
            color=tc, fontweight='bold' if bold else 'normal',
            multialignment='center', zorder=3,
            fontfamily='sans-serif')

def arrow(x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color=PAL['gray_dark'],
                                lw=0.75, mutation_scale=8))

# ── Source boxes (top row) ────────────────────────────────────────────────────
y_src = 0.92;  h_src = 0.065
box(0.18, y_src, 0.27, h_src,
    "PubChem\n95 AIDs  ·  341,328 records", PAL['blue'], bold=True, fs=6.5)
box(0.50, y_src, 0.22, h_src,
    "ChEMBL\nCHEMBL6111  ·  4,925 rows",   PAL['teal'], bold=True, fs=6.5)
box(0.82, y_src, 0.27, h_src,
    "BindingDB\nQ9UM07  ·  3,087 rows",     PAL['orange'], bold=True, fs=6.5)

# ── Stage 01 ─────────────────────────────────────────────────────────────────
y01 = 0.805; h01 = 0.065
box(0.50, y01, 0.60, h01,
    "Step 01 — SMILES Standardization & InChIKey\n341,328 input rows → 341,282 standardized",
    '#F5F5F5', tc='black', ec='#CCCCCC')
# three converging arrows meeting box top edge at equal x spacing
for sx in (0.18, 0.50, 0.82):
    tx = sx + (0.50 - sx) * 0.5   # midpoint x, then straight down
    ax.annotate('', xy=(sx + (0.50-sx)*0.65, y01 + h01/2),
                xytext=(sx, y_src - h_src/2),
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color=PAL['gray_dark'],
                                lw=0.65, mutation_scale=7,
                                connectionstyle='arc3,rad=0.0'))
# straighten center arrow
ax.annotate('', xy=(0.50, y01 + h01/2), xytext=(0.50, y_src - h_src/2),
            xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', color=PAL['gray_dark'],
                            lw=0.65, mutation_scale=7))

# ── Stage 02 ─────────────────────────────────────────────────────────────────
y02 = 0.685; h02 = 0.068
box(0.50, y02, 0.60, h02,
    "Step 02 — Activity Normalization\n6-layer architecture (A–F)  ·  IC50 → pIC50",
    '#F5F5F5', tc='black', ec='#CCCCCC')
arrow(0.50, y01 - h01/2, 0.50, y02 + h02/2)

# ── Stage 03 ─────────────────────────────────────────────────────────────────
y03 = 0.565; h03 = 0.068
box(0.50, y03, 0.60, h03,
    "Step 03 — HTS Extraction & Potency Filter\n7,815 pIC50-qualified measurements retained",
    '#F5F5F5', tc='black', ec='#CCCCCC')
arrow(0.50, y02 - h02/2, 0.50, y03 + h03/2)

# ── Split → HTS + Potency ────────────────────────────────────────────────────
y_sp = 0.425; h_sp = 0.080
box(0.19, y_sp, 0.32, h_sp,
    "HTS Space\n327,336 structural references\n(reference only)",
    PAL['gray_light'], tc=PAL['gray_dark'], ec='#CCCCCC', fs=6.5)
box(0.75, y_sp, 0.36, h_sp,
    "Potency Space\n7,815 → 7,319 rows (replicate aggregation)",
    '#E8F4F8', tc=PAL['gray_dark'], ec='#BBDDEE', fs=6.5)
# diverging arrows from Step 03 bottom center
arrow(0.50, y03 - h03/2, 0.22, y_sp + h_sp/2)
arrow(0.50, y03 - h03/2, 0.68, y_sp + h_sp/2)
ax.text(0.19, y_sp - h_sp/2 - 0.012,
        "(no further processing)", ha='center', va='top',
        fontsize=5.5, color=PAL['gray_dark'],
        style='italic', transform=ax.transAxes)

# ── Stage 04 ─────────────────────────────────────────────────────────────────
y04 = 0.29; h04 = 0.068
box(0.50, y04, 0.60, h04,
    "Step 04 — Deduplication + Assembly\n450 replicate groups  ·  independence scoring  ·  3,093 compounds",
    '#F5F5F5', tc='black', ec='#CCCCCC')
arrow(0.75, y_sp - h_sp/2, 0.65, y04 + h04/2)

# ── Stage 05 ─────────────────────────────────────────────────────────────────
y05 = 0.165; h05 = 0.068
box(0.50, y05, 0.60, h05,
    "Step 05 — Scaffold + Activity Cliff Analysis\n1,244 scaffolds  ·  358K pairs  ·  94 severe cliffs",
    '#F5F5F5', tc='black', ec='#CCCCCC')
arrow(0.50, y04 - h04/2, 0.50, y05 + h05/2)

# ── Final output ──────────────────────────────────────────────────────────────
y_fin = 0.048; h_fin = 0.075
box(0.50, y_fin, 0.64, h_fin,
    "PAD4-DB v2\n3,093 curated PAD4 inhibitors  ·  25 columns",
    PAL['navy'], bold=True, fs=8.0, lw=1.0, ec=PAL['navy'])
arrow(0.50, y05 - h05/2, 0.50, y_fin + h_fin/2)

save_fig(fig, 'fig1_pipeline')
plt.close()
print("Fig 1 DONE")
