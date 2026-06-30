#!/usr/bin/env python
"""Figure 4 — pIC50 Distribution Analysis (2×2 panels) + Great Tables."""
import os, sys, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

# ── Library standards ────────────────────────────────────────────────────────
import scienceplots
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
plt.style.use(['science', 'nature', 'no-latex'])

from great_tables import GT, loc, style as gt_style
import great_tables

import numpy as np
import pandas as pd
import scipy.stats as stats
from adjustText import adjust_text
import matplotlib
matplotlib.use('Agg')

os.makedirs('outputs/figures', exist_ok=True)
os.makedirs('outputs/tables', exist_ok=True)

# ── Load data ────────────────────────────────────────────────────────────────
df = pd.read_parquet('data/processed/pad4_compounds.parquet')
assert len(df) == 3093

nodes = pd.read_csv('outputs/figures/fig7_nodes.csv')
assert len(nodes) == 99

HUB_A = {'SMADULGDNOCLOP-GISFHXKWSA-N', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N'}
HUB_B = {'UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'DVCKJOQIVOGXEI-XMMPIXPASA-N'}

# ── Source category helper ───────────────────────────────────────────────────
def get_src_cat(sl):
    if sl == 'pubchem_confirmatory':  return 'PubChem only (patent)'
    if sl == 'bindingdb':             return 'BindingDB only'
    if sl == 'chembl':                return 'ChEMBL only'
    return 'Multi-source'

df['_src_cat'] = df['source_list'].map(get_src_cat)

MECH_ORDER = ['enzymatic', 'enzymatic_confirmed', 'fp_ic50', 'covalent']
MECH_COLORS = {
    'enzymatic':           '#AAAAAA',
    'enzymatic_confirmed': '#4A90D9',
    'fp_ic50':             '#2ECC71',
    'covalent':            '#E05A2B',
}
MECH_LABELS = {
    'enzymatic':           'Enzymatic\n(n=2,079)',
    'enzymatic_confirmed': 'Enzymatic\nConfirmed\n(n=878)',
    'fp_ic50':             'FP IC50\n(n=115)',
    'covalent':            'Covalent\n(n=21)',
}
SRC_COLORS = {
    'PubChem only (patent)': '#E05A2B',
    'BindingDB only':        '#4A90D9',
    'ChEMBL only':           '#2ECC71',
    'Multi-source':          '#AAAAAA',
}

# ── Verification before plotting ─────────────────────────────────────────────
pic50 = df['pic50_consensus'].values
mean_val   = float(np.mean(pic50))
median_val = float(np.median(pic50))
std_val    = float(np.std(pic50, ddof=1))

src_counts = df['_src_cat'].value_counts().to_dict()
mech_counts = df['mechanism_class'].value_counts().to_dict()

hub_a_mask  = nodes['hub_class'] == 'A'
hub_b_mask  = nodes['hub_class'] == 'B'
nonhub_mask = nodes['hub_class'] == 'none'
max_degree  = int(nodes['severe_cliff_degree'].max())
max_degree_ik = nodes.loc[nodes['severe_cliff_degree'].idxmax(), 'inchi_key']

print("=== Pre-plot verification ===")
print(f"Panel A: n={len(df):,}, mean={mean_val:.3f}, median={median_val:.3f}, std={std_val:.3f}")
print(f"Panel B counts: " +
      " | ".join(f"{k}={src_counts.get(k, 0)}" for k in SRC_COLORS))
print(f"Panel C counts: " +
      " | ".join(f"{k}={mech_counts.get(k, 0)}" for k in MECH_ORDER))
print(f"Panel D: total={len(nodes)}, hub_A={hub_a_mask.sum()}, "
      f"hub_B={hub_b_mask.sum()}, non-hub={nonhub_mask.sum()}")
print(f"Max degree: {max_degree} — {max_degree_ik}")
print()

EXPECTED_SRC  = {'PubChem only (patent)': 233, 'BindingDB only': 95,
                 'ChEMBL only': 10, 'Multi-source': 2755}
EXPECTED_MECH = {'enzymatic': 2079, 'enzymatic_confirmed': 878,
                 'fp_ic50': 115, 'covalent': 21}
errs = []
for k, v in EXPECTED_SRC.items():
    if src_counts.get(k) != v:
        errs.append(f"SRC {k}: expected {v}, got {src_counts.get(k)}")
for k, v in EXPECTED_MECH.items():
    if mech_counts.get(k) != v:
        errs.append(f"MECH {k}: expected {v}, got {mech_counts.get(k)}")
if len(nodes) != 99 or hub_a_mask.sum() != 2 or hub_b_mask.sum() != 2:
    errs.append(f"Panel D node counts wrong")
if errs:
    for e in errs: print("FAIL:", e)
    sys.exit(1)
print("All pre-plot checks PASS")
print()

# ── Figure ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(10, 9))
ax_a, ax_b = axes[0, 0], axes[0, 1]
ax_c, ax_d = axes[1, 0], axes[1, 1]

def panel_label(ax, letter):
    ax.text(0.02, 0.96, letter, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top')

# ── Panel A: Overall distribution ────────────────────────────────────────────
counts, bin_edges, _ = ax_a.hist(
    pic50, bins=40, range=(2.0, 8.52),
    color='#4A90D9', alpha=0.6, edgecolor='white', linewidth=0.4,
    label='_nolegend_',
)

ax_a_r = ax_a.twinx()
kde_x = np.linspace(1.8, 8.8, 500)
kde = stats.gaussian_kde(pic50)
ax_a_r.plot(kde_x, kde(kde_x), color='#1A237E', linewidth=1.5, zorder=5)
ax_a_r.set_ylabel('Density', fontsize=9)
ax_a_r.set_ylim(bottom=0)
ax_a_r.spines['right'].set_visible(True)

ax_a.axvline(mean_val,   color='#E74C3C', linestyle='--', linewidth=1.0,
             label=f'Mean = {mean_val:.2f}')
ax_a.axvline(median_val, color='#F39C12', linestyle='--', linewidth=1.0,
             label=f'Median = {median_val:.2f}')
ax_a.set_xlim(1.8, 8.8)
ax_a.set_xlabel('pIC50', fontsize=9)
ax_a.set_ylabel('Count', fontsize=9)
ax_a.legend(loc='upper left', fontsize=7.5, framealpha=0.8)
ax_a.text(0.97, 0.95, 'n = 3,093', transform=ax_a.transAxes,
          ha='right', va='top', fontsize=8)
panel_label(ax_a, 'A')

# ── Panel B: KDE by source category ──────────────────────────────────────────
src_order = ['Multi-source', 'PubChem only (patent)', 'BindingDB only', 'ChEMBL only']
for cat in src_order:
    mask = df['_src_cat'] == cat
    vals = df.loc[mask, 'pic50_consensus'].values
    n = mask.sum()
    color = SRC_COLORS[cat]
    if len(vals) < 3:
        continue
    kde_cat = stats.gaussian_kde(vals, bw_method='scott')
    y_vals = kde_cat(kde_x)
    ax_b.plot(kde_x, y_vals, color=color, linewidth=1.5,
              label=f'{cat} (n={n})')
    ax_b.fill_between(kde_x, y_vals, alpha=0.15, color=color)

ax_b.set_xlim(1.8, 8.8)
ax_b.set_xlabel('pIC50', fontsize=9)
ax_b.set_ylabel('Density', fontsize=9)
ax_b.legend(loc='upper left', fontsize=6.5, framealpha=0.8)
panel_label(ax_b, 'B')

# ── Panel C: Violin by mechanism_class ───────────────────────────────────────
violin_data = [df.loc[df['mechanism_class'] == m, 'pic50_consensus'].values
               for m in MECH_ORDER]

vp = ax_c.violinplot(violin_data, positions=range(len(MECH_ORDER)),
                     showmedians=False, showextrema=True, widths=0.7)

for i, (body, mech) in enumerate(zip(vp['bodies'], MECH_ORDER)):
    body.set_facecolor(MECH_COLORS[mech])
    body.set_alpha(0.7)
    body.set_linewidth(0.8)

vp['cbars'].set_linewidth(0.8)
vp['cmaxes'].set_linewidth(0.8)
vp['cmins'].set_linewidth(0.8)

# Median lines inside each violin
for i, vals in enumerate(violin_data):
    med = np.median(vals)
    ax_c.hlines(med, i - 0.2, i + 0.2, color='white', linewidth=1.5, zorder=5)

# Strip points for fp_ic50 and covalent only
rng = np.random.default_rng(42)
for i, mech in enumerate(MECH_ORDER):
    if mech in ('fp_ic50', 'covalent'):
        vals = df.loc[df['mechanism_class'] == mech, 'pic50_consensus'].values
        jitter = rng.uniform(-0.08, 0.08, size=len(vals))
        ax_c.scatter(np.full_like(vals, i) + jitter, vals,
                     s=12, color=MECH_COLORS[mech], alpha=0.8,
                     edgecolor='white', linewidth=0.3, zorder=6)

ax_c.axhline(6.0, color='gray', linestyle='--', linewidth=1.0, alpha=0.4)
ax_c.set_xticks(range(len(MECH_ORDER)))
ax_c.set_xticklabels([MECH_LABELS[m] for m in MECH_ORDER],
                     rotation=15, ha='right', fontsize=8)
ax_c.set_ylabel('pIC50', fontsize=9)
panel_label(ax_c, 'C')

# ── Panel D: pIC50 vs cliff degree for 99 severe cliff compounds ──────────────
nodes['_hub'] = nodes['hub_class']

mask_nh  = nodes['_hub'] == 'none'
mask_ha  = nodes['_hub'] == 'A'
mask_hb  = nodes['_hub'] == 'B'

ax_d.scatter(nodes.loc[mask_nh, 'severe_cliff_degree'],
             nodes.loc[mask_nh, 'pic50_consensus'],
             c='#AAAAAA', s=40, alpha=0.7, zorder=2, label='_nolegend_')
ax_d.scatter(nodes.loc[mask_hb, 'severe_cliff_degree'],
             nodes.loc[mask_hb, 'pic50_consensus'],
             c='#1A237E', s=200, marker='*', zorder=10,
             label='Hub Class B (n=2)')
ax_d.scatter(nodes.loc[mask_ha, 'severe_cliff_degree'],
             nodes.loc[mask_ha, 'pic50_consensus'],
             c='#E74C3C', s=200, marker='*', zorder=11,
             label='Hub Class A (n=2)')

# Text annotations for all 4 hubs
texts = []
hub_nodes = nodes[nodes['_hub'] != 'none']
for _, row in hub_nodes.iterrows():
    short_ik = row['inchi_key'][:14]
    txt = ax_d.text(row['severe_cliff_degree'], row['pic50_consensus'],
                    short_ik, fontsize=7)
    texts.append(txt)
adjust_text(texts, ax=ax_d, arrowprops=dict(arrowstyle='-', color='gray',
                                             lw=0.5))

ax_d.axvline(5, color='gray', linestyle='--', linewidth=1.0, alpha=0.4,
             label='degree = 5')
ax_d.set_xlabel('Severe Cliff Degree (number of pairs)', fontsize=9)
ax_d.set_ylabel('pIC50', fontsize=9)
ax_d.set_xlim(0, max_degree + 1)
ax_d.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
ax_d.legend(loc='upper right', fontsize=7.5, framealpha=0.8)
panel_label(ax_d, 'D')

plt.tight_layout(pad=0.8)

PNG_PATH = 'outputs/figures/fig4_pic50_distribution.png'
SVG_PATH = 'outputs/figures/fig4_pic50_distribution.svg'
fig.savefig(PNG_PATH, dpi=300, bbox_inches='tight')
fig.savefig(SVG_PATH, bbox_inches='tight')
plt.close()
print(f"Saved: {PNG_PATH}  ({os.path.getsize(PNG_PATH)/1024:.1f} KB)")
print(f"Saved: {SVG_PATH}  ({os.path.getsize(SVG_PATH)/1024:.1f} KB)")
print()

# ── Great Tables summary ─────────────────────────────────────────────────────
rows = []
for mech in MECH_ORDER:
    vals = df.loc[df['mechanism_class'] == mech, 'pic50_consensus']
    rows.append({
        'mechanism_class': mech,
        'n':               int(len(vals)),
        'mean_pic50':      float(vals.mean()),
        'median_pic50':    float(vals.median()),
        'std_pic50':       float(vals.std(ddof=1)),
        'min':             float(vals.min()),
        'max':             float(vals.max()),
    })
gt_df = pd.DataFrame(rows)

gt = (
    GT(gt_df)
    .tab_header(
        title="PAD4-DB pIC50 Distribution by Mechanism Class",
        subtitle="3,093 curated PAD4 inhibitors",
    )
    .cols_label(
        mechanism_class="Mechanism Class",
        n="N",
        mean_pic50="Mean pIC50",
        median_pic50="Median pIC50",
        std_pic50="SD",
        min="Min",
        max="Max",
    )
    .fmt_number(
        columns=["mean_pic50", "median_pic50", "std_pic50", "min", "max"],
        decimals=2,
    )
    .cols_align(
        align="center",
        columns=["n", "mean_pic50", "median_pic50", "std_pic50", "min", "max"],
    )
    .tab_style(
        style=gt_style.fill(color="#E8F4F8"),
        locations=loc.body(rows=[1]),
    )
    .tab_source_note(
        "pIC50 = −log10(IC50 in M). "
        "Mechanism class derived from assay metadata and SMARTS warhead detection."
    )
)

HTML_PATH = 'outputs/tables/fig4_distribution_stats.html'
TEX_PATH  = 'outputs/tables/fig4_distribution_stats.tex'
with open(HTML_PATH, 'w') as f:
    f.write(gt.as_raw_html())
with open(TEX_PATH, 'w') as f:
    f.write(gt.as_latex())
print(f"Great Tables HTML: {HTML_PATH}  ({os.path.getsize(HTML_PATH)/1024:.1f} KB)")
print(f"Great Tables LaTeX: {TEX_PATH}  ({os.path.getsize(TEX_PATH)/1024:.1f} KB)")
print()

# ── Completion report ────────────────────────────────────────────────────────
print("=== Completion Report ===")
print(f"Panel A median: {median_val:.3f}")
print(f"Panel D max degree: {max_degree_ik}, degree={max_degree}")
print()
print("Files written:")
for p in [PNG_PATH, SVG_PATH, HTML_PATH, TEX_PATH]:
    print(f"  {p}  ({os.path.getsize(p)/1024:.1f} KB)")
