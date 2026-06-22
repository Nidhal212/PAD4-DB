"""
supp_s04_reference_recovery.py — S4: Reference compound concordance scatter
Seven recovered PAD4 inhibitors: published pIC50 vs PAD4-DB v2 consensus pIC50.
Outputs: publication/figures/supplementary/fig_s04_reference_recovery.{png,pdf}
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, save_fig, C, SINGLE

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

set_style()

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT_SUPP = ROOT / 'publication/figures/supplementary'
OUT_SUPP.mkdir(parents=True, exist_ok=True)

COMPOUNDS = [
    # name,          published_pic50, db_pic50,  jbi_outlier, gsk_offset
    ('Streptonigrin', 5.60,           5.602,     False,        False),
    ('Cl-amidine',    5.23,           5.219,     False,        False),
    ('F-Amidine',     4.67,           4.571,     False,        False),
    ('GSK484',        7.30,           7.049,     False,        True),
    ('TDFA',          5.64,           5.638,     False,        False),
    ('BMS-P5',        7.01,           7.009,     False,        False),
    ('JBI-589',       6.914,          6.000,     True,         False),
]

names        = [c[0] for c in COMPOUNDS]
x_published  = np.array([c[1] for c in COMPOUNDS])
y_db         = np.array([c[2] for c in COMPOUNDS])
is_jbi       = np.array([c[3] for c in COMPOUNDS])
is_gsk       = np.array([c[4] for c in COMPOUNDS])
is_concordant = ~is_jbi & ~is_gsk

deltas = np.abs(x_published - y_db)
mean_delta_excl_jbi = deltas[~is_jbi].mean()

print("=" * 60)
print("SUPPLEMENTARY S4 — REFERENCE COMPOUND RECOVERY")
print("=" * 60)
for name, xp, yd, jbi, gsk in COMPOUNDS:
    delta = abs(xp - yd)
    flag = ''
    if jbi: flag = ' *** JBI-589 (Ca2+ assay outlier, |d|=0.91)'
    elif gsk: flag = ' ** GSK484 (inter-assay, |d|=0.25)'
    print(f"  {name:<18} published={xp:.3f}  DB={yd:.3f}  |d|={delta:.3f}{flag}")
print(f"  Concordant (|d|<=0.15): {is_concordant.sum()} / 7")
print(f"  Mean |dpIC50| (6 non-JBI): {mean_delta_excl_jbi:.4f}")

fig, ax = plt.subplots(figsize=(SINGLE, SINGLE))

ax_min, ax_max = 4.0, 8.0

ax.plot([ax_min, ax_max], [ax_min, ax_max], color=C['grey'], lw=0.8, ls='--', zorder=1)
ax.fill_between([ax_min, ax_max],
                [ax_min - 0.15, ax_max - 0.15],
                [ax_min + 0.15, ax_max + 0.15],
                color=C['blue'], alpha=0.08, zorder=1)

# 5 concordant
ax.scatter(x_published[is_concordant], y_db[is_concordant],
           s=40, c=C['blue'], zorder=4, edgecolors='white', lw=0.4)

# GSK484 open circle
ax.scatter(x_published[is_gsk], y_db[is_gsk],
           s=45, facecolors='none', edgecolors=C['blue'], lw=1.2, zorder=4)

# JBI-589 red diamond
ax.scatter(x_published[is_jbi], y_db[is_jbi],
           s=60, c=C['vermillion'], marker='D', zorder=5, edgecolors='white', lw=0.4)

label_offsets = {
    'Streptonigrin': (-0.05, -0.17),
    'Cl-amidine':    (+0.03, +0.07),
    'F-Amidine':     (+0.03, +0.07),
    'GSK484':        (+0.03, -0.16),
    'TDFA':          (+0.03, +0.07),
    'BMS-P5':        (+0.03, +0.07),
    'JBI-589':       (+0.05, +0.07),
}
for name, xp, yd, jbi, gsk in COMPOUNDS:
    dx, dy = label_offsets.get(name, (+0.03, +0.07))
    color = C['vermillion'] if jbi else C['grey']
    ax.text(xp + dx, yd + dy, name, fontsize=5.5, color=color, va='bottom', ha='left')

x_gsk = float(x_published[is_gsk])
y_gsk = float(y_db[is_gsk])
ax.annotate('|ΔpIC50|=0.25\n(inter-assay)',
    xy=(x_gsk, y_gsk), xytext=(x_gsk - 0.6, y_gsk - 0.45),
    fontsize=5, color=C['blue'],
    arrowprops=dict(arrowstyle='->', color=C['blue'], lw=0.7),
    ha='right', va='top')

x_jbi = float(x_published[is_jbi])
y_jbi = float(y_db[is_jbi])
ax.annotate('|ΔpIC50|=0.91\n(assay Ca2+)',
    xy=(x_jbi, y_jbi), xytext=(x_jbi + 0.25, y_jbi - 0.55),
    fontsize=5, color=C['vermillion'],
    arrowprops=dict(arrowstyle='->', color=C['vermillion'], lw=0.7),
    ha='left', va='top')

ax.text(0.04, 0.97,
        f'5/7: |ΔpIC50| ≤ 0.15  ●\n'
        f'GSK484: |Δ|=0.25 (inter-assay)  ○\n'
        f'mean (6, excl. JBI-589) = {mean_delta_excl_jbi:.3f}',
        transform=ax.transAxes, fontsize=5.0, va='top', ha='left', color=C['grey'])

ax.set_xlabel('Published pIC50')
ax.set_ylabel('PAD4-DB v2 consensus pIC50')
ax.set_xlim(ax_min, ax_max)
ax.set_ylim(ax_min, ax_max)
ax.set_aspect('equal')

legend_handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=C['blue'],
           markersize=5, label='Concordant, |Δ|≤0.15 (n=5)'),
    Line2D([0], [0], marker='o', color=C['blue'], markerfacecolor='none',
           markersize=5, markeredgewidth=1.2, label='GSK484, |Δ|=0.25'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor=C['vermillion'],
           markersize=5, label='JBI-589, |Δ|=0.91 (Ca2+)'),
    Line2D([0], [0], color=C['grey'], lw=0.8, ls='--', label='y = x'),
]
ax.legend(handles=legend_handles, fontsize=6, frameon=False, loc='lower right', borderaxespad=1)

fig.tight_layout(pad=0.5)
save_fig(fig, str(OUT_SUPP / 'fig_s04_reference_recovery'))
plt.close(fig)
print("S4 reference recovery complete.")
