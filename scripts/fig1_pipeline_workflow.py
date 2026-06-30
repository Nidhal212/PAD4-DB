#!/usr/bin/env python
"""Figure 1 — Pipeline Workflow Diagram (programmatic matplotlib flowchart)."""
import os, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import scienceplots
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
plt.style.use(['science', 'nature', 'no-latex'])

os.makedirs('outputs/figures', exist_ok=True)

fig, ax = plt.subplots(figsize=(8, 10))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')


def draw_box(ax, x, y, w, h, text, facecolor, textcolor='black',
             fontsize=8.5, bold=False):
    box = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle='round,pad=0.05',
        facecolor=facecolor, edgecolor='#AAAAAA', linewidth=0.8,
        transform=ax.transAxes, clip_on=False, zorder=2,
    )
    ax.add_patch(box)
    ax.text(x, y, text, transform=ax.transAxes,
            ha='center', va='center', fontsize=fontsize,
            color=textcolor, fontweight='bold' if bold else 'normal',
            multialignment='center', zorder=3)


def arrow(ax, x1, y1, x2, y2):
    ax.annotate(
        '', xy=(x2, y2), xytext=(x1, y1),
        xycoords='axes fraction', textcoords='axes fraction',
        arrowprops=dict(arrowstyle='->', color='#555555', lw=1.2),
        zorder=1,
    )


# Row 1: Source boxes
y_src = 0.91;  h_src = 0.076
draw_box(ax, 0.18, y_src, 0.28, h_src,
         "PubChem\n95 AIDs\n341,328 records",   '#4A90D9', 'white', fontsize=8)
draw_box(ax, 0.50, y_src, 0.24, h_src,
         "ChEMBL\nCHEMBL6111\n4,925 rows",      '#2ECC71', 'white', fontsize=8)
draw_box(ax, 0.82, y_src, 0.28, h_src,
         "BindingDB\nQ9UM07\n3,087 rows",        '#F39C12', 'white', fontsize=8)

# Stage 01
y_s01 = 0.76;  h_s01 = 0.075
draw_box(ax, 0.50, y_s01, 0.56, h_s01,
         "Step 01 — SMILES Standardization\n341,282 standardized entries",
         '#EEEEEE', 'black', fontsize=8.5)
arrow(ax, 0.18, y_src - h_src/2, 0.37, y_s01 + h_s01/2)
arrow(ax, 0.50, y_src - h_src/2, 0.50, y_s01 + h_s01/2)
arrow(ax, 0.82, y_src - h_src/2, 0.63, y_s01 + h_s01/2)

# Stage 02
y_s02 = 0.61;  h_s02 = 0.085
draw_box(ax, 0.50, y_s02, 0.56, h_s02,
         "Step 02 — Activity Normalization\n6-layer architecture (A–F)\npIC50 conversion",
         '#EEEEEE', 'black', fontsize=8.5)
arrow(ax, 0.50, y_s01 - h_s01/2, 0.50, y_s02 + h_s02/2)

# Stage 03 split
y_s03 = 0.46;  h_s03 = 0.08
draw_box(ax, 0.20, y_s03, 0.32, h_s03,
         "HTS Space\n327,336 structural\nreferences",
         '#DDDDDD', 'black', fontsize=7.5)
draw_box(ax, 0.76, y_s03, 0.36, h_s03,
         "Potency Space\n7,815 measurements\n→ 7,319 rows",
         '#DDDDDD', 'black', fontsize=7.5)
arrow(ax, 0.37, y_s02 - h_s02/2, 0.23, y_s03 + h_s03/2)
arrow(ax, 0.63, y_s02 - h_s02/2, 0.73, y_s03 + h_s03/2)
ax.text(0.20, y_s03 - h_s03/2 - 0.012, "(structural reference only)",
        transform=ax.transAxes, ha='center', va='top',
        fontsize=6, color='#888888', style='italic')

# Stage 04
y_s04 = 0.315; h_s04 = 0.085
draw_box(ax, 0.50, y_s04, 0.56, h_s04,
         "Step 04 — Deduplication + Assembly\nReplicate aggregation (450 groups)\nSource independence scoring",
         '#EEEEEE', 'black', fontsize=8.5)
arrow(ax, 0.76, y_s03 - h_s03/2, 0.63, y_s04 + h_s04/2)

# Stage 05
y_s05 = 0.165; h_s05 = 0.08
draw_box(ax, 0.50, y_s05, 0.56, h_s05,
         "Step 05 — Scaffold + Cliff Analysis\n1,244 Bemis-Murcko scaffolds · 94 severe cliff pairs",
         '#EEEEEE', 'black', fontsize=8.5)
arrow(ax, 0.50, y_s04 - h_s04/2, 0.50, y_s05 + h_s05/2)

# Final output
y_fin = 0.038; h_fin = 0.075
draw_box(ax, 0.50, y_fin, 0.62, h_fin,
         "PAD4-DB v2\n3,093 curated PAD4 inhibitors · 25 columns",
         '#E74C3C', 'white', fontsize=9.5, bold=True)
arrow(ax, 0.50, y_s05 - h_s05/2, 0.50, y_fin + h_fin/2)

# Audit annotation (right side, rotated)
ax.text(0.965, 0.53,
        "Audit Trail\n10-phase master audit\n58/58 checks passed",
        transform=ax.transAxes,
        ha='center', va='center', fontsize=7, color='#666666',
        style='italic', multialignment='center', rotation=90)

ax.set_title("PAD4-DB v2 Data Processing Pipeline", fontsize=11, pad=6)

PNG_PATH = 'outputs/figures/fig1_pipeline_workflow.png'
SVG_PATH = 'outputs/figures/fig1_pipeline_workflow.svg'
plt.savefig(PNG_PATH, dpi=300, bbox_inches='tight')
plt.savefig(SVG_PATH, bbox_inches='tight')
plt.close()

print(f"Saved: {PNG_PATH}  ({os.path.getsize(PNG_PATH)/1024:.1f} KB)")
print(f"Saved: {SVG_PATH}  ({os.path.getsize(SVG_PATH)/1024:.1f} KB)")
print("TASK A: DONE")
