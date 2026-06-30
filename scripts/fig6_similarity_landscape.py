#!/usr/bin/env python
"""Figure 6 — Similarity / Activity Landscape (2×2 panels) + Great Tables."""
import os, sys, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

# ── Library standards ─────────────────────────────────────────────────────────
import scienceplots
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.style.use(['science', 'nature', 'no-latex'])

from great_tables import GT, loc, style as gt_style
import great_tables

import numpy as np
import pandas as pd

os.makedirs('outputs/figures', exist_ok=True)
os.makedirs('outputs/tables', exist_ok=True)

print(f"SciencePlots: importable")
print(f"Great Tables:  {great_tables.__version__}")
print()

# ── Load data ─────────────────────────────────────────────────────────────────
pairs  = pd.read_parquet('data/processed/activity_pairs_with_sali.parquet')
cliffs = pd.read_parquet('data/processed/activity_cliffs.parquet')

# ── Locked number verification ────────────────────────────────────────────────
n_total       = len(pairs)
n_sali_valid  = int(pairs['sali'].notna().sum())
n_tan1        = int((pairs['tanimoto'] == 1.0).sum())
sali_max      = float(pairs['sali'].max())
n_severe      = int((cliffs['cliff_tier'] == 'severe').sum())
n_moderate    = int((cliffs['cliff_tier'] == 'moderate').sum())
n_broad       = int((cliffs['cliff_tier'] == 'broad').sum())

n_sim_ge08    = int((pairs['tanimoto'] >= 0.8).sum())
n_sali_gt10   = int((pairs['sali'] > 10).sum())
n_sali_gt20   = int((pairs['sali'] > 20).sum())

print("=== Locked number verification ===")
LOCKED = [
    ("Total pairs (sim≥0.6)",     n_total,      358416),
    ("Pairs with SALI (non-NaN)", n_sali_valid, 358218),
    ("SALI max (≈65.88)",         round(sali_max, 2), 65.88),
    ("Severe cliff pairs",        n_severe,     94),
    ("Moderate cliff pairs",      n_moderate,   193),
    ("Broad cliff pairs",         n_broad,      580),
]
all_pass = True
for label, actual, expected in LOCKED:
    ok = actual == expected
    print(f"  {label:35s}: {actual}  {'PASS' if ok else f'FAIL (expected {expected})'}")
    if not ok: all_pass = False
if not all_pass:
    print("\nVERIFICATION FAILED — stopping.")
    sys.exit(1)
print()
print(f"Additional counts (derived):")
print(f"  Pairs sim≥0.8:  {n_sim_ge08:,}")
print(f"  Pairs SALI>10:  {n_sali_gt10:,}")
print(f"  Pairs SALI>20:  {n_sali_gt20:,}")
print(f"  Pairs sim=1.0:  {n_tan1}")
print()

# ── Subset references ─────────────────────────────────────────────────────────
valid = pairs[pairs['sali'].notna()]          # 358,218 rows
severe_cliffs   = cliffs[cliffs['cliff_tier'] == 'severe']
moderate_cliffs = cliffs[cliffs['cliff_tier'] == 'moderate']

# Severe/moderate pairs from main pairs frame (has tanimoto column needed)
severe_pairs   = pairs[pairs['cliff_tier'] == 'severe']
moderate_pairs = pairs[pairs['cliff_tier'] == 'moderate']

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
ax_a, ax_b = axes[0, 0], axes[0, 1]
ax_c, ax_d = axes[1, 0], axes[1, 1]

def panel_label(ax, letter):
    ax.text(0.02, 0.96, letter, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top')

# ── Panel A: Tanimoto histogram ───────────────────────────────────────────────
tan = pairs['tanimoto'].values.astype(float)
ax_a.hist(tan, bins=50, range=(0.6, 1.0),
          color='#4A90D9', alpha=0.7, edgecolor='white', linewidth=0.3)
ax_a.axvline(0.8, color='#E74C3C', linestyle='--', linewidth=1.0,
             label='Cliff threshold (0.8)')
# Shade region ≥0.8
ylims_a = ax_a.get_ylim()
ax_a.axvspan(0.8, 1.0, alpha=0.08, color='#E74C3C', zorder=0)
# Annotate shaded region
ax_a.text(0.88, ax_a.get_ylim()[1] * 0.65,
          f"n={n_sim_ge08:,}\npairs",
          ha='center', fontsize=8, color='#E74C3C')
ax_a.set_xlabel('Tanimoto Similarity (ECFP4)', fontsize=9)
ax_a.set_ylabel('Number of Pairs', fontsize=9)
ax_a.legend(loc='upper right', fontsize=8, framealpha=0.8)
panel_label(ax_a, 'A')

# ── Panel B: ΔpIC50 vs Tanimoto scatter ───────────────────────────────────────
rng = np.random.default_rng(42)
sample_idx = rng.choice(len(pairs), size=50000, replace=False)
samp = pairs.iloc[sample_idx]

ax_b.scatter(samp['tanimoto'], samp['delta_pic50'],
             c='#CCCCCC', s=1, alpha=0.3, rasterized=True, zorder=1)
ax_b.scatter(moderate_pairs['tanimoto'], moderate_pairs['delta_pic50'],
             c='#F39C12', s=8, alpha=0.7, zorder=4,
             label=f'Moderate cliff (n={len(moderate_pairs)})')
ax_b.scatter(severe_pairs['tanimoto'], severe_pairs['delta_pic50'],
             c='#E74C3C', s=15, alpha=0.9, zorder=5,
             label=f'Severe cliff (n={len(severe_pairs)})')
ax_b.axhline(2.0, color='#E74C3C', linestyle='--', linewidth=0.8, zorder=3)
ax_b.axhline(1.5, color='#F39C12', linestyle='--', linewidth=0.8, zorder=3)
ax_b.axvline(0.8, color='gray',    linestyle='--', linewidth=0.8, zorder=3)
ax_b.set_xlim(0.55, 1.05)
ax_b.set_ylim(-0.5, 3.5)
ax_b.set_xlabel('Tanimoto Similarity (ECFP4)', fontsize=9)
ax_b.set_ylabel('ΔpIC50', fontsize=9)
ax_b.legend(loc='upper left', fontsize=7.5, framealpha=0.8, markerscale=2)
panel_label(ax_b, 'B')

# ── Panel C: SALI histogram (log scale) ───────────────────────────────────────
sali_vals = valid['sali'].values
ax_c.hist(sali_vals, bins=80, range=(0, 70),
          color='#4A90D9', alpha=0.7, edgecolor='white', linewidth=0.3)
ax_c.set_yscale('log')
ax_c.axvline(10, color='#F39C12', linestyle='--', linewidth=1.0, label='SALI=10')
ax_c.axvline(20, color='#E74C3C', linestyle='--', linewidth=1.0, label='SALI=20')
# Annotate counts (after log scale is set, get y position from data range)
ymax_c = ax_c.get_ylim()[1]
ax_c.text(0.97, 0.92, f'SALI>10: n={n_sali_gt10:,}',
          transform=ax_c.transAxes, ha='right', fontsize=8, color='#F39C12')
ax_c.text(0.97, 0.82, f'SALI>20: n={n_sali_gt20:,}',
          transform=ax_c.transAxes, ha='right', fontsize=8, color='#E74C3C')
ax_c.text(0.97, 0.72, f'Max SALI={sali_max:.2f}',
          transform=ax_c.transAxes, ha='right', fontsize=8, color='#1A237E')
ax_c.set_xlabel('SALI', fontsize=9)
ax_c.set_ylabel('Number of Pairs (log scale)', fontsize=9)
ax_c.legend(loc='upper right', fontsize=8, framealpha=0.8,
            bbox_to_anchor=(0.97, 0.65))
panel_label(ax_c, 'C')

# ── Panel D: SAS 2D density map ───────────────────────────────────────────────
h2d = ax_d.hist2d(
    valid['tanimoto'].values.astype(float),
    valid['delta_pic50'].values,
    bins=60,
    range=[[0.6, 1.0], [0, 3.5]],
    cmap='Blues',
    norm=LogNorm(),
)
cb_d = plt.colorbar(h2d[3], ax=ax_d, pad=0.02)
cb_d.set_label('Pair density (log)', fontsize=8)
cb_d.ax.tick_params(labelsize=7)

# Reference lines
ax_d.axvline(0.8, color='white', linestyle='--', linewidth=0.8, alpha=0.8)
ax_d.axhline(2.0, color='white', linestyle='--', linewidth=0.8, alpha=0.8)
ax_d.axhline(1.5, color='white', linestyle='--', linewidth=0.6, alpha=0.6)

# Quadrant labels
ax_d.text(0.905, 2.65, "Severe\nActivity\nCliffs",
          fontsize=7, color='white', ha='center', fontweight='bold')
ax_d.text(0.68,  2.65, "Distant\nActivity\nCliffs",
          fontsize=7, color='white', ha='center', fontweight='bold')
ax_d.text(0.905, 0.35, "Similar,\nConcordant",
          fontsize=7, color='white', ha='center')
ax_d.text(0.68,  0.35, "Diverse,\nConcordant",
          fontsize=7, color='gray',  ha='center')

ax_d.set_xlim(0.6, 1.0)
ax_d.set_ylim(0,   3.5)
ax_d.set_xlabel('Tanimoto Similarity (ECFP4)', fontsize=9)
ax_d.set_ylabel('ΔpIC50', fontsize=9)
panel_label(ax_d, 'D')

plt.tight_layout(pad=0.8)

PNG_PATH = 'outputs/figures/fig6_similarity_landscape.png'
SVG_PATH = 'outputs/figures/fig6_similarity_landscape.svg'
fig.savefig(PNG_PATH, dpi=300, bbox_inches='tight')
fig.savefig(SVG_PATH, bbox_inches='tight')
plt.close()
print(f"Saved: {PNG_PATH}  ({os.path.getsize(PNG_PATH)/1024:.1f} KB)")
print(f"Saved: {SVG_PATH}  ({os.path.getsize(SVG_PATH)/1024:.1f} KB)")
print()

# ── Great Tables: SALI tier summary ──────────────────────────────────────────
def cliff_tier_overlap(sali_lo, sali_hi, df):
    sub = df[(df['sali'] > sali_lo) & (df['sali'] <= sali_hi) & df['cliff_tier'].notna()]
    tiers = sub['cliff_tier'].value_counts()
    if len(tiers) == 0:
        return 'none'
    return ', '.join(f'{t}={n}' for t, n in tiers.items())

tiers_def = [
    ('Low',     0,  5,   '0–5'),
    ('Medium',  5,  10,  '5–10'),
    ('High',    10, 20,  '10–20'),
    ('Extreme', 20, 999, '>20'),
]
rows = []
for name, lo, hi, rng_str in tiers_def:
    mask = (valid['sali'] > lo) & (valid['sali'] <= hi) if hi < 999 else (valid['sali'] > lo)
    sub = valid[mask]
    ct = cliff_tier_overlap(lo, hi, valid)
    rows.append({
        'tier':             name,
        'sali_range':       rng_str,
        'n_pairs':          len(sub),
        'pct_total':        len(sub) / len(valid) * 100,
        'mean_tanimoto':    float(sub['tanimoto'].mean()),
        'mean_delta_pic50': float(sub['delta_pic50'].mean()),
        'example_cliff_tier': ct,
    })
gt_df = pd.DataFrame(rows)

gt = (
    GT(gt_df)
    .tab_header(
        title="PAD4-DB SALI Distribution Summary",
        subtitle=f"358,218 compound pairs with Tanimoto ≥ 0.6 and Tanimoto < 1.0",
    )
    .cols_label(
        tier="SALI Tier",
        sali_range="Range",
        n_pairs="N Pairs",
        pct_total="% Total",
        mean_tanimoto="Mean Tanimoto",
        mean_delta_pic50="Mean ΔpIC50",
        example_cliff_tier="Cliff tier overlap",
    )
    .fmt_number(
        columns=["pct_total", "mean_tanimoto", "mean_delta_pic50"],
        decimals=2,
    )
    .fmt_integer(columns=["n_pairs"])
    .tab_style(
        style=gt_style.fill(color="#FFEBEE"),
        locations=loc.body(rows=[3]),
    )
    .tab_source_note(
        "SALI = |ΔpIC50| / (1 − Tanimoto). "
        "Pairs with Tanimoto = 1.0 (n=198) excluded from SALI computation."
    )
)

HTML_PATH = 'outputs/tables/fig6_sali_stats.html'
TEX_PATH  = 'outputs/tables/fig6_sali_stats.tex'
with open(HTML_PATH, 'w') as f:
    f.write(gt.as_raw_html())
with open(TEX_PATH, 'w') as f:
    f.write(gt.as_latex())
print(f"Great Tables HTML:  {HTML_PATH}  ({os.path.getsize(HTML_PATH)/1024:.1f} KB)")
print(f"Great Tables LaTeX: {TEX_PATH}  ({os.path.getsize(TEX_PATH)/1024:.1f} KB)")
print()

# ── Completion report ─────────────────────────────────────────────────────────
print("=== Completion Report ===")
print(f"Pairs sim≥0.8:   {n_sim_ge08:,}")
print(f"Pairs SALI>10:   {n_sali_gt10:,}")
print(f"Pairs SALI>20:   {n_sali_gt20:,}")
print(f"SALI max:        {sali_max:.2f}  {'✓' if round(sali_max,2)==65.88 else 'FAIL'}")
print()
print("Files written:")
for p in [PNG_PATH, SVG_PATH, HTML_PATH, TEX_PATH]:
    print(f"  {p}  ({os.path.getsize(p)/1024:.1f} KB)")
