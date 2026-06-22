"""
supp_s01_pipeline.py — S1: PAD4-DB v2 Curation Pipeline (complete rewrite)
Outputs: publication/figures/supplementary/fig_s01_pipeline.{png,pdf}
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, save_fig, C

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

set_style()

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT  = ROOT / 'publication/figures/supplementary'
OUT.mkdir(parents=True, exist_ok=True)

# ── Design constants (data coordinates: x in [0,1], y in [0,1]) ───────────────
FIG_W, FIG_H = 6.0, 9.5      # inches — single tall column
BOX_X0  = 0.13               # left edge of boxes
BOX_W   = 0.78               # box width
BOX_H   = 0.095              # default box height
ACCENT  = 0.022              # left accent bar width
BADGE_R = 0.022              # badge circle radius
GAP     = 0.028              # vertical gap between boxes
ARROW_X = BOX_X0 + BOX_W / 2

# Accent colors per step
CLR_SRC   = C['blue']
CLR_STEP  = C['navy']
CLR_SAR   = C['teal']
CLR_OUT   = C['orange']

# ── Step definitions ──────────────────────────────────────────────────────────
steps = [
    {
        'badge': '00', 'accent': CLR_SRC,
        'title': 'Source Databases',
        'body': 'PubChem: 95 AIDs (HTS · confirmatory · literature · secondary)\n'
                'ChEMBL: CHEMBL6111 assay  ·  BindingDB: UniProt Q9UM07',
    },
    {
        'badge': '01', 'accent': CLR_STEP,
        'title': 'Step 01 — Standardize SMILES',
        'body': 'RDKit salt strip, neutralize, canonicalize → InChIKey generation\n'
                '341,282 rows processed  ·  328,976 unique InChIKeys  ·  6 NO_SMILES',
    },
    {
        'badge': '02', 'accent': CLR_STEP,
        'title': 'Step 02 — Normalize Activities',
        'body': 'IC50 → pIC50 conversion (nM)  ·  endpoint typing (Pct_inh excluded)\n'
                '338,021 OK (99.0%)  ·  7,815 use_in_potency_model rows',
    },
    {
        'badge': '03', 'accent': CLR_STEP,
        'title': 'Step 03 — Replicate Aggregation',
        'body': 'Group by InChIKey × source × AID × endpoint_type  ·  geometric mean\n'
                '450 multi-replicate groups consolidated  ·  7,319 potency rows retained',
    },
    {
        'badge': '04', 'accent': CLR_STEP,
        'title': 'Step 04 — Deduplication & Assembly',
        'body': 'InChIKey cross-source dedup  ·  AID 2202576/77 preferred (n=55)\n'
                '3,093 curated PAD4 inhibitors  ·  source_independence_score computed',
    },
    {
        'badge': '05', 'accent': CLR_SAR,
        'title': 'Step 05 — SAR Analysis',
        'body': 'ECFP4 fingerprints (r=2, 2048 bits)  ·  Murcko scaffold analysis\n'
                '1,244 scaffolds  ·  94 severe cliffs  ·  MMP analysis  ·  ECFP6 sensitivity',
    },
    {
        'badge': 'v2', 'accent': CLR_OUT,
        'title': 'PAD4-DB v2 — Final Dataset',
        'body': '3,093 curated PAD4 inhibitors  ·  pIC50 range 2.00–8.52 (median 6.84)\n'
                '1,244 Murcko scaffolds  ·  94 severe activity cliffs  ·  MMP-validated 85.1%',
        'highlight': True,
    },
]

# ── Compute y positions from top down ─────────────────────────────────────────
n = len(steps)
total_boxes  = n * BOX_H
total_gaps   = (n - 1) * GAP
total_height = total_boxes + total_gaps
y_start      = (1.0 - total_height) / 2 + total_height   # top of first box

ys = []
y = y_start
for i, step in enumerate(steps):
    ys.append(y - BOX_H)
    y -= BOX_H + GAP

# ── Draw ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

for i, (step, y0) in enumerate(zip(steps, ys)):
    accent = step['accent']
    is_highlight = step.get('highlight', False)

    # Main box
    face = '#FFF8F0' if is_highlight else 'white'
    edge = accent if is_highlight else '#DDDDDD'
    lw   = 1.5 if is_highlight else 0.8
    rect = mpatches.FancyBboxPatch(
        (BOX_X0, y0), BOX_W, BOX_H,
        boxstyle='round,pad=0.008',
        facecolor=face, edgecolor=edge, linewidth=lw,
        transform=ax.transAxes, clip_on=False, zorder=2
    )
    ax.add_patch(rect)

    # Left accent bar
    accent_bar = mpatches.FancyBboxPatch(
        (BOX_X0, y0), ACCENT, BOX_H,
        boxstyle='round,pad=0.004',
        facecolor=accent, edgecolor='none',
        transform=ax.transAxes, clip_on=False, zorder=3
    )
    ax.add_patch(accent_bar)

    # Badge circle
    badge_x = BOX_X0 + ACCENT + 0.038
    badge_y = y0 + BOX_H / 2
    circ = plt.Circle(
        (badge_x, badge_y), BADGE_R,
        color=accent, transform=ax.transAxes,
        clip_on=False, zorder=4
    )
    ax.add_patch(circ)
    ax.text(badge_x, badge_y, step['badge'],
            transform=ax.transAxes, fontsize=6.5, fontweight='bold',
            ha='center', va='center', color='white', zorder=5, clip_on=False)

    # Title
    text_x = BOX_X0 + ACCENT + 0.08
    title_color = accent if is_highlight else '#111111'
    ax.text(text_x, y0 + BOX_H * 0.68, step['title'],
            transform=ax.transAxes, fontsize=7.5, fontweight='bold',
            ha='left', va='center', color=title_color, clip_on=False, zorder=4)

    # Body
    ax.text(text_x, y0 + BOX_H * 0.28, step['body'],
            transform=ax.transAxes, fontsize=6, ha='left', va='center',
            color='#555555', clip_on=False, zorder=4,
            linespacing=1.4)

    # Downward arrow to next box
    if i < len(steps) - 1:
        arrow_y_top = y0 - 0.002
        arrow_y_bot = ys[i + 1] + BOX_H + 0.002
        ax.annotate('',
            xy=(ARROW_X, arrow_y_bot),
            xytext=(ARROW_X, arrow_y_top),
            xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', color='#AAAAAA', lw=1.2),
        )

# Legend row at very bottom
legend_y = min(ys) - 0.055
ax.add_patch(plt.Circle((0.22, legend_y + 0.012), 0.009,
    color=CLR_SRC, transform=ax.transAxes, clip_on=False))
ax.text(0.235, legend_y + 0.012, 'Source', transform=ax.transAxes,
        fontsize=6, va='center', color='#555555', clip_on=False)

ax.add_patch(plt.Circle((0.38, legend_y + 0.012), 0.009,
    color=CLR_STEP, transform=ax.transAxes, clip_on=False))
ax.text(0.395, legend_y + 0.012, 'Curation step', transform=ax.transAxes,
        fontsize=6, va='center', color='#555555', clip_on=False)

ax.add_patch(plt.Circle((0.57, legend_y + 0.012), 0.009,
    color=CLR_SAR, transform=ax.transAxes, clip_on=False))
ax.text(0.585, legend_y + 0.012, 'SAR analysis', transform=ax.transAxes,
        fontsize=6, va='center', color='#555555', clip_on=False)

ax.add_patch(plt.Circle((0.74, legend_y + 0.012), 0.009,
    color=CLR_OUT, transform=ax.transAxes, clip_on=False))
ax.text(0.755, legend_y + 0.012, 'Output', transform=ax.transAxes,
        fontsize=6, va='center', color='#555555', clip_on=False)

save_fig(fig, str(OUT / 'fig_s01_pipeline'))
plt.close(fig)
print("S1 complete.")
