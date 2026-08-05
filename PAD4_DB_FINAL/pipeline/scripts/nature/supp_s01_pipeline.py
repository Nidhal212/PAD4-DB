"""
supp_s01_pipeline.py — S1: Pipeline workflow diagram
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Liberation Sans', 'DejaVu Sans', 'Arial'],
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

COLORS = {
    'blue': '#0077BB', 'navy': '#004488', 'grey': '#BBBBBB',
    'dark_grey': '#555555', 'light_grey': '#E8E8E8',
    'teal': '#009988', 'orange': '#EE7733',
}

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/figures/nature_v2'
OUT.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(6, 8))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

boxes = [
    {
        'title': 'Source Databases',
        'body': 'PubChem (95 AIDs)  ·  ChEMBL (CHEMBL6111)\nBindingDB (Q9UM07)',
        'badge': '00',
        'border': COLORS['blue'],
        'badge_color': COLORS['blue'],
        'y': 0.84,
    },
    {
        'title': 'Step 01 — Standardize',
        'body': 'RDKit salt strip, neutralize, canonicalize\nSMILES standardization',
        'badge': '01',
        'border': COLORS['grey'],
        'badge_color': COLORS['teal'],
        'y': 0.65,
    },
    {
        'title': 'Step 02 — Deduplicate',
        'body': 'pIC50 consensus  ·  3,093 compounds\nSource independence scoring',
        'badge': '02',
        'border': COLORS['grey'],
        'badge_color': COLORS['teal'],
        'y': 0.46,
    },
    {
        'title': 'Step 03 — Analysis',
        'body': '358,416 pairs  ·  94 severe cliffs\nMMP analysis  ·  Scaffold series',
        'badge': '03',
        'border': COLORS['grey'],
        'badge_color': COLORS['teal'],
        'y': 0.27,
    },
    {
        'title': 'PAD4-DB v2',
        'body': '3,093 curated PAD4 inhibitors\n1,244 scaffolds  ·  94 activity cliffs',
        'badge': 'v2',
        'border': COLORS['navy'],
        'badge_color': COLORS['navy'],
        'y': 0.06,
        'navy_text': True,
    },
]

box_w = 0.72
box_h = 0.13
badge_r = 0.035
badge_x = 0.07

for i, box in enumerate(boxes):
    x0 = 0.15
    y0 = box['y']

    # Box
    rect = mpatches.FancyBboxPatch(
        (x0, y0), box_w, box_h,
        boxstyle='round,pad=0.01',
        facecolor=COLORS['light_grey'],
        edgecolor=box['border'], lw=2.5,
        transform=ax.transAxes, clip_on=False
    )
    ax.add_patch(rect)

    # Badge circle
    badge_circ = plt.Circle(
        (badge_x, y0 + box_h / 2), badge_r,
        color=box['badge_color'], transform=ax.transAxes, clip_on=False
    )
    ax.add_patch(badge_circ)
    ax.text(badge_x, y0 + box_h / 2, box['badge'],
            transform=ax.transAxes, fontsize=8, fontweight='bold',
            ha='center', va='center', color='white', clip_on=False)

    # Title
    title_color = COLORS['navy'] if box.get('navy_text') else COLORS['dark_grey']
    ax.text(x0 + 0.04, y0 + box_h * 0.72, box['title'],
            transform=ax.transAxes, fontsize=9, fontweight='bold',
            ha='left', va='center', color=title_color)

    # Body
    ax.text(x0 + 0.04, y0 + box_h * 0.28, box['body'],
            transform=ax.transAxes, fontsize=7.5,
            ha='left', va='center', color=COLORS['dark_grey'])

    # Arrow to next box
    if i < len(boxes) - 1:
        next_y = boxes[i + 1]['y'] + box_h
        cur_y = y0
        arrow_x = x0 + box_w / 2
        ax.annotate('', xy=(arrow_x, next_y + 0.005),
                     xytext=(arrow_x, cur_y - 0.005),
                     arrowprops=dict(arrowstyle='->', color=COLORS['grey'],
                                     lw=2.0),
                     xycoords='axes fraction', textcoords='axes fraction')

ax.set_title('PAD4-DB v2: Curation Pipeline', fontsize=12, fontweight='bold',
              color=COLORS['navy'], pad=10)

for ext in ['png', 'pdf']:
    outpath = OUT / f'supp_s01_pipeline.{ext}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"  Saved: {outpath}")

plt.close(fig)
print("S1 complete.")
