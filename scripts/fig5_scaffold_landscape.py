#!/usr/bin/env python
"""Figure 5 — Scaffold Landscape (2×2 panels) + Great Tables."""
import os, sys, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

# ── Library standards ────────────────────────────────────────────────────────
import scienceplots
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.style.use(['science', 'nature', 'no-latex'])

from great_tables import GT, loc, style as gt_style
import great_tables

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

os.makedirs('outputs/figures', exist_ok=True)
os.makedirs('outputs/tables', exist_ok=True)

print(f"SciencePlots: importable")
print(f"Great Tables: {great_tables.__version__}")
print()

# ── Load data ────────────────────────────────────────────────────────────────
df = pd.read_parquet('data/processed/pad4_compounds.parquet')
sc = pd.read_csv('outputs/tables/05_scaffold_summary.csv')
coords = np.load('data/interim/tsne_coords_3093.npy')
assert len(df) == 3093 and coords.shape == (3093, 2)

# ── Derive scaffold SMILES per compound ──────────────────────────────────────
print("Deriving Murcko scaffolds for 3,093 compounds...")
scaffold_smi_list = []
for smi in df['smiles_std']:
    try:
        mol = Chem.MolFromSmiles(smi) if pd.notna(smi) else None
        if mol is None:
            scaffold_smi_list.append(None)
        else:
            scaffold_smi_list.append(
                MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
            )
    except Exception:
        scaffold_smi_list.append(None)
df['scaffold_smi'] = scaffold_smi_list
print(f"  Done. Null scaffolds: {df['scaffold_smi'].isna().sum()}")

# ── Rank scaffolds descending by n_compounds ─────────────────────────────────
sc_sorted = sc.sort_values('n_compounds', ascending=False).reset_index(drop=True)
sc_sorted['rank'] = sc_sorted.index + 1

# Map compound → scaffold rank
smi_to_rank = sc_sorted.set_index('scaffold_smiles')['rank'].to_dict()
smi_to_n    = sc_sorted.set_index('scaffold_smiles')['n_compounds'].to_dict()
df['scaffold_rank'] = df['scaffold_smi'].map(smi_to_rank)
df['scaffold_n']    = df['scaffold_smi'].map(smi_to_n).fillna(0).astype(int)

# ── Compute patent-exclusive count per scaffold ───────────────────────────────
# source_list contains 'pubchem_confirmatory' only for patent-exclusive
df['_patent_only'] = df['source_list'] == 'pubchem_confirmatory'
patent_counts = (
    df[df['_patent_only']]
    .groupby('scaffold_smi')['inchi_key']
    .count()
    .rename('n_patent_exclusive')
)
sc_sorted['n_patent_exclusive'] = (
    sc_sorted['scaffold_smiles'].map(
        df.groupby('scaffold_smi')['_patent_only'].sum()
    ).fillna(0).astype(int)
)

# ── Locked number verification ───────────────────────────────────────────────
series_sc = sc_sorted[sc_sorted['n_compounds'] >= 2]
single_sc = sc_sorted[sc_sorted['n_compounds'] == 1]
total_unique = len(sc_sorted)
n_series     = len(series_sc)
n_singletons = len(single_sc)
largest      = int(sc_sorted['n_compounds'].iloc[0])
coverage_n   = int(series_sc['n_compounds'].sum())
coverage_pct = coverage_n / 3093 * 100

print("\n=== Locked number verification ===")
checks = [
    ("Total unique scaffolds", total_unique, 1244),
    ("Series scaffolds (≥2)",  n_series,     375),
    ("Singleton scaffolds",    n_singletons, 869),
    ("Largest series",         largest,      174),
]
all_pass = True
for label, actual, expected in checks:
    ok = actual == expected
    print(f"  {label:35s}: {actual}  {'PASS' if ok else f'FAIL (expected {expected})'}")
    if not ok: all_pass = False
# Coverage: allow ±0.1%
cov_ok = abs(coverage_pct - 71.9) < 0.2
print(f"  {'Scaffold series coverage':35s}: {coverage_n}/3093 = {coverage_pct:.1f}%  "
      f"{'PASS' if cov_ok else 'FAIL (expected 71.8–71.9%)'}")
if not cov_ok: all_pass = False
if not all_pass:
    print("\nVERIFICATION FAILED — stopping.")
    sys.exit(1)
print()

# ── Panel D scaffold groups ───────────────────────────────────────────────────
HUB_A_IKS = {'SMADULGDNOCLOP-GISFHXKWSA-N', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N'}
mask_rank1    = df['scaffold_rank'] == 1          # azaindole-benzimidazole hub scaffold
mask_top2_10  = df['scaffold_rank'].between(2, 10)
mask_hub_a    = df['inchi_key'].isin(HUB_A_IKS)
mask_other    = ~(mask_rank1 | mask_top2_10)

print(f"Panel D groups:")
print(f"  Hub scaffold series (rank 1): {mask_rank1.sum()}  (expected 174)")
print(f"  Top 2–10 series:              {mask_top2_10.sum()}")
print(f"  Other/singleton:              {mask_other.sum()}")
print(f"  Hub Class A compounds:        {mask_hub_a.sum()}  (expected 2)")
print()

# ── Build figure ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
ax_a, ax_b = axes[0, 0], axes[0, 1]
ax_c, ax_d = axes[1, 0], axes[1, 1]

def panel_label(ax, letter):
    ax.text(0.02, 0.96, letter, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top')

# ── Panel A: Scaffold frequency (top 30) ─────────────────────────────────────
top30 = sc_sorted.head(30)
ranks   = top30['rank'].values
counts  = top30['n_compounds'].values
colors  = ['#E74C3C' if r == 1 else '#4A90D9' for r in ranks]

ax_a.bar(range(1, 31), counts, color=colors, width=0.8,
         edgecolor='white', linewidth=0.3)
# Annotate rank-1 bar
ax_a.text(1, counts[0] + 3,
          f"174 cpds\n(Hub Class A\nscaffold)",
          ha='center', va='bottom', fontsize=7, color='#E74C3C',
          fontweight='bold')
ax_a.set_xlabel('Scaffold Rank', fontsize=9)
ax_a.set_ylabel('Number of Compounds', fontsize=9)
ax_a.set_xlim(0.4, 30.6)
ax_a.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax_a.text(0.97, 0.95, 'Top 30 of 375 series scaffolds',
          transform=ax_a.transAxes, ha='right', va='top', fontsize=8)
panel_label(ax_a, 'A')

# ── Panel B: Scaffold size distribution histogram ─────────────────────────────
series_sizes = series_sc['n_compounds'].values
med_size = float(np.median(series_sizes))
# Log-spaced bins from 2 to max
bins = np.logspace(np.log10(2), np.log10(largest + 1), 30)
ax_b.hist(series_sizes, bins=bins, color='#4A90D9', alpha=0.7,
          edgecolor='white', linewidth=0.5)
ax_b.set_xscale('log')
ax_b.axvline(med_size, color='#E05A2B', linestyle='--', linewidth=1.2,
             label=f'Median = {med_size:.0f}')
ax_b.text(med_size * 1.15, ax_b.get_ylim()[1] * 0.85,
          f'Median = {med_size:.0f}', color='#E05A2B', fontsize=8)
ax_b.set_xlabel('Scaffold Series Size (log scale)', fontsize=9)
ax_b.set_ylabel('Number of Scaffolds', fontsize=9)
ax_b.text(0.97, 0.95, f'375 series scaffolds\n869 singletons',
          transform=ax_b.transAxes, ha='right', va='top', fontsize=8)
panel_label(ax_b, 'B')

# ── Panel C: Gini + Lorenz curve ─────────────────────────────────────────────
all_counts = np.sort(sc_sorted['n_compounds'].values)  # ascending
n_sc = len(all_counts)
cumsum_c = np.cumsum(all_counts) / all_counts.sum()
lorenz_y = np.concatenate([[0], cumsum_c])
lorenz_x = np.linspace(0, 1, len(lorenz_y))
gini = 1.0 - 2.0 * np.trapz(lorenz_y, lorenz_x)

ax_c.plot(lorenz_x, lorenz_y, color='#4A90D9', linewidth=1.5,
          label=f'PAD4-DB (Gini = {gini:.3f})')
ax_c.plot([0, 1], [0, 1], 'k--', linewidth=0.8, label='Perfect equality')
ax_c.fill_between(lorenz_x, lorenz_y, lorenz_x, alpha=0.15, color='#4A90D9')
ax_c.set_xlabel('Cumulative Fraction of Scaffolds', fontsize=9)
ax_c.set_ylabel('Cumulative Fraction of Compounds', fontsize=9)
ax_c.legend(loc='lower right', fontsize=8, framealpha=0.8)
ax_c.text(0.05, 0.85, f'Gini = {gini:.3f}',
          transform=ax_c.transAxes, fontsize=11, fontweight='bold',
          color='#4A90D9', va='top')
panel_label(ax_c, 'C')

# ── Panel D: Scaffold groups on t-SNE ────────────────────────────────────────
x, y = coords[:, 0], coords[:, 1]

# Layer order: other → top2-10 → rank1 → hub A stars
ax_d.scatter(x[mask_other],   y[mask_other],   c='#DDDDDD', s=4,  alpha=0.3,
             zorder=1, rasterized=True, label='Other / singleton scaffolds')
ax_d.scatter(x[mask_top2_10], y[mask_top2_10], c='#4A90D9', s=8,  alpha=0.5,
             zorder=3, rasterized=True, label='Top 2–10 scaffold series')
ax_d.scatter(x[mask_rank1],   y[mask_rank1],   c='#E74C3C', s=15, alpha=0.8,
             zorder=5, rasterized=True, label='Hub scaffold series (n=174)')
ax_d.scatter(x[mask_hub_a],   y[mask_hub_a],   c='black',   s=150, marker='*',
             zorder=10, label='Hub Class A compounds')

ax_d.set_xlabel('t-SNE 1', fontsize=9)
ax_d.set_ylabel('t-SNE 2', fontsize=9)
ax_d.set_xticks([])
ax_d.set_yticks([])
ax_d.legend(loc='upper right', fontsize=7, framealpha=0.8, markerscale=1.3)
panel_label(ax_d, 'D')

plt.tight_layout(pad=0.8)

PNG_PATH = 'outputs/figures/fig5_scaffold_landscape.png'
SVG_PATH = 'outputs/figures/fig5_scaffold_landscape.svg'
fig.savefig(PNG_PATH, dpi=300, bbox_inches='tight')
fig.savefig(SVG_PATH, bbox_inches='tight')
plt.close()
print(f"Saved: {PNG_PATH}  ({os.path.getsize(PNG_PATH)/1024:.1f} KB)")
print(f"Saved: {SVG_PATH}  ({os.path.getsize(SVG_PATH)/1024:.1f} KB)")
print()

# ── Great Tables: top 20 scaffolds ───────────────────────────────────────────
top20 = sc_sorted.head(20).copy()
top20['rank']             = top20.index + 1
top20['pct_dataset']      = top20['n_compounds'] / 3093 * 100
top20['scaffold_smiles_t']= top20['scaffold_smiles'].str[:40]
top20['contains_hub']     = top20['rank'] == 1

gt_df = top20[[
    'rank', 'scaffold_smiles_t', 'n_compounds', 'pct_dataset',
    'mean_pic50', 'std_pic50', 'n_patent_exclusive', 'contains_hub',
]].rename(columns={'scaffold_smiles_t': 'scaffold_smiles'})

gt = (
    GT(gt_df)
    .tab_header(
        title="PAD4-DB Top 20 Scaffold Series",
        subtitle="Ranked by compound count; 375 total series (≥2 compounds)",
    )
    .cols_label(
        rank="Rank",
        scaffold_smiles="Scaffold (truncated)",
        n_compounds="N",
        pct_dataset="% Dataset",
        mean_pic50="Mean pIC50",
        std_pic50="SD pIC50",
        n_patent_exclusive="Patent-exclusive",
        contains_hub="Hub scaffold",
    )
    .fmt_number(columns=["mean_pic50", "std_pic50"], decimals=2)
    .fmt_number(columns=["pct_dataset"], decimals=1)
    .tab_style(
        style=gt_style.fill(color="#FFEBEE"),
        locations=loc.body(rows=[0]),
    )
    .tab_style(
        style=gt_style.text(weight="bold"),
        locations=loc.body(columns=["n_compounds"], rows=[0]),
    )
    .tab_source_note(
        "Bemis-Murcko scaffolds computed with RDKit. "
        "Hub scaffold = azaindole-benzimidazole series containing "
        "Hub Class A compounds."
    )
)

HTML_PATH = 'outputs/tables/fig5_scaffold_stats.html'
TEX_PATH  = 'outputs/tables/fig5_scaffold_stats.tex'
with open(HTML_PATH, 'w') as f:
    f.write(gt.as_raw_html())
with open(TEX_PATH, 'w') as f:
    f.write(gt.as_latex())
print(f"Great Tables HTML:  {HTML_PATH}  ({os.path.getsize(HTML_PATH)/1024:.1f} KB)")
print(f"Great Tables LaTeX: {TEX_PATH}  ({os.path.getsize(TEX_PATH)/1024:.1f} KB)")
print()

# ── Completion report ─────────────────────────────────────────────────────────
print("=== Completion Report ===")
print(f"Gini coefficient:             {gini:.3f}")
print(f"Median scaffold series size:  {med_size:.0f}")
rank1_n = int(sc_sorted.iloc[0]['n_compounds'])
print(f"Rank 1 compound count:        {rank1_n}  {'✓' if rank1_n == 174 else 'FAIL'}")
print(f"Scaffold series coverage:     {coverage_pct:.1f}%  {'✓' if cov_ok else 'FAIL'}")
print()
print("Files written:")
for p in [PNG_PATH, SVG_PATH, HTML_PATH, TEX_PATH]:
    print(f"  {p}  ({os.path.getsize(p)/1024:.1f} KB)")
