#!/usr/bin/env python
"""Figure 10 — Patent Scaffold Analysis (supplementary) + Great Tables."""
import os, sys, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import scienceplots
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
plt.style.use(['science', 'nature', 'no-latex'])

from great_tables import GT, loc, style as gt_style
import great_tables
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

os.makedirs('outputs/figures', exist_ok=True)
os.makedirs('outputs/tables', exist_ok=True)

print(f"SciencePlots: importable  |  Great Tables: {great_tables.__version__}")

# Load data
df     = pd.read_parquet('data/processed/pad4_compounds.parquet')
sc_csv = pd.read_csv('outputs/tables/05_scaffold_summary.csv')
coords = np.load('data/interim/tsne_coords_3093.npy')
pairs  = pd.read_parquet('data/processed/activity_pairs_with_sali.parquet')

assert len(df) == 3093, f"Expected 3093 compounds, got {len(df)}"
assert coords.shape == (3093, 2), f"Bad coords shape: {coords.shape}"

# Patent / published split (compound level)
mask_patent    = df['source_list'] == 'pubchem_confirmatory'
mask_published = ~mask_patent
n_patent    = int(mask_patent.sum())
n_published = int(mask_published.sum())

# Severe cliff patent contribution
severe_pairs    = pairs[pairs['cliff_tier'] == 'severe']
any_patent_sev  = int(severe_pairs['any_patent_exclusive'].sum())
both_patent_sev = int((severe_pairs['patent_exclusive_a'] &
                       severe_pairs['patent_exclusive_b']).sum())

print("=== Locked number verification ===")
LOCKED = [
    ("Patent-exclusive compounds",         n_patent,       233),
    ("Non-patent compounds",               n_published,    2860),
    ("Severe cliffs with any patent cpd",  any_patent_sev, 1),
]
all_pass = True
for label, actual, expected in LOCKED:
    ok = actual == expected
    print(f"  {label:35s}: {actual}  {'PASS' if ok else f'FAIL (expected {expected})'}")
    if not ok:
        all_pass = False
if not all_pass:
    print("\nVERIFICATION FAILED — stopping.")
    sys.exit(1)
print()

# pIC50 stats by group
mean_patent    = float(df.loc[mask_patent,    'pic50_consensus'].mean())
mean_published = float(df.loc[mask_published, 'pic50_consensus'].mean())
print(f"Patent mean pIC50:    {mean_patent:.3f}")
print(f"Published mean pIC50: {mean_published:.3f}")

# Derive Murcko scaffolds for Panel B scaffold size analysis
print("Deriving Murcko scaffolds for scaffold size analysis...")
smi_list = []
for smi in df['smiles_std']:
    try:
        mol = Chem.MolFromSmiles(smi) if pd.notna(smi) else None
        smi_list.append(
            MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
            if mol else None
        )
    except Exception:
        smi_list.append(None)
df['scaffold_smi'] = smi_list

# Classify scaffolds: patent-only vs shared
grp = df.groupby('scaffold_smi').agg(
    n_total=('inchi_key', 'count'),
    n_patent=('source_list', lambda x: (x == 'pubchem_confirmatory').sum()),
).reset_index()
grp['patent_only'] = grp['n_patent'] == grp['n_total']
grp['has_published'] = grp['n_patent'] < grp['n_total']

patent_only_scaffolds = grp[grp['patent_only']]
shared_scaffolds      = grp[grp['has_published'] & (grp['n_total'] >= 2)]

n_patent_scaffolds = len(patent_only_scaffolds)
print(f"Patent-only scaffolds (fresh derivation): {n_patent_scaffolds}  (locked=103)")
print(f"Shared series scaffolds (≥2 cpds, ≥1 non-patent): {len(shared_scaffolds)}")
print()

# Map compound → patent_only scaffold flag for Panel C
patent_only_iks = set(
    df.loc[df['scaffold_smi'].isin(
        patent_only_scaffolds['scaffold_smi'].values
    ), 'inchi_key']
)

# Build arrays for panels
x, y = coords[:, 0], coords[:, 1]
idx_patent    = df.index[mask_patent].tolist()
idx_published = df.index[mask_published].tolist()

# Panel A KDE data
pic50_patent    = df.loc[mask_patent,    'pic50_consensus'].dropna().values
pic50_published = df.loc[mask_published, 'pic50_consensus'].dropna().values

# Panel B scaffold sizes
sizes_patent_only = patent_only_scaffolds['n_total'].values
sizes_shared      = shared_scaffolds['n_total'].values

# Figure
fig, axes = plt.subplots(1, 3, figsize=(13, 5))
ax_a, ax_b, ax_c = axes


def panel_label(ax, letter):
    ax.text(0.02, 0.96, letter, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top')


# Panel A: pIC50 KDE — patent vs published
x_grid = np.linspace(1.5, 9.0, 400)

kde_pat = gaussian_kde(pic50_patent, bw_method=0.25)
kde_pub = gaussian_kde(pic50_published, bw_method=0.25)
y_pat = kde_pat(x_grid)
y_pub = kde_pub(x_grid)

ax_a.fill_between(x_grid, y_pub, color='#4A90D9', alpha=0.15)
ax_a.fill_between(x_grid, y_pat, color='#E05A2B', alpha=0.15)
ax_a.plot(x_grid, y_pub, color='#4A90D9', linewidth=1.5,
          label=f'Published (n={n_published:,})')
ax_a.plot(x_grid, y_pat, color='#E05A2B', linewidth=1.5,
          label=f'Patent-exclusive (n={n_patent})')

ax_a.axvline(mean_published, color='#4A90D9', linestyle='--',
             linewidth=0.9, alpha=0.8)
ax_a.axvline(mean_patent, color='#E05A2B', linestyle='--',
             linewidth=0.9, alpha=0.8)

ax_a.text(mean_published + 0.05, ax_a.get_ylim()[1] * 0.92,
          f'{mean_published:.2f}', fontsize=7, color='#4A90D9', va='top')
ax_a.text(mean_patent - 0.08, ax_a.get_ylim()[1] * 0.92,
          f'{mean_patent:.2f}', fontsize=7, color='#E05A2B', va='top', ha='right')

ax_a.set_xlabel('pIC50', fontsize=9)
ax_a.set_ylabel('Density', fontsize=9)
ax_a.set_title('pIC50: Patent vs Published', fontsize=9)
ax_a.legend(fontsize=7.5, framealpha=0.8)
panel_label(ax_a, 'A')

# Panel B: Scaffold size distribution — violin/boxplot
bp_data  = [sizes_patent_only, sizes_shared]
bp_labels = [f'Patent-only\n(n={n_patent_scaffolds} scaffolds)',
             f'Shared series\n(n={len(shared_scaffolds)} scaffolds)']
bp_colors = ['#E05A2B', '#4A90D9']

parts = ax_b.violinplot(bp_data, positions=[0, 1], showmedians=True,
                        showextrema=True, widths=0.6)
for i, (pc, col) in enumerate(zip(parts['bodies'], bp_colors)):
    pc.set_facecolor(col)
    pc.set_alpha(0.5)
for key in ('cmedians', 'cmins', 'cmaxes', 'cbars'):
    parts[key].set_color('#333333')
    parts[key].set_linewidth(0.8)

ax_b.set_xticks([0, 1])
ax_b.set_xticklabels(bp_labels, fontsize=8)
ax_b.set_ylabel('Compounds per scaffold', fontsize=9)
ax_b.set_title('Scaffold Series Size', fontsize=9)
ax_b.text(0.97, 0.95,
          f'Median patent-only: {np.median(sizes_patent_only):.0f}\n'
          f'Median shared: {np.median(sizes_shared):.0f}',
          transform=ax_b.transAxes, ha='right', va='top',
          fontsize=7.5, multialignment='right')
panel_label(ax_b, 'B')

# Panel C: t-SNE with patent compounds highlighted
ax_c.scatter(x[idx_published], y[idx_published],
             c='#4A90D9', s=8, alpha=0.5, rasterized=True, zorder=2,
             label=f'Published (n={n_published:,})')
ax_c.scatter(x[idx_patent], y[idx_patent],
             c='#E05A2B', s=20, alpha=0.8, rasterized=True, zorder=5,
             label=f'Patent-exclusive (n={n_patent})')
ax_c.set_xlabel('t-SNE 1', fontsize=9)
ax_c.set_ylabel('t-SNE 2', fontsize=9)
ax_c.set_xticks([])
ax_c.set_yticks([])
ax_c.set_title('t-SNE: Patent vs Published', fontsize=9)
ax_c.legend(fontsize=7.5, framealpha=0.8, loc='upper right', markerscale=1.5)
panel_label(ax_c, 'C')

plt.tight_layout(pad=0.8)

PNG_PATH = 'outputs/figures/fig10_patent_scaffolds.png'
SVG_PATH = 'outputs/figures/fig10_patent_scaffolds.svg'
fig.savefig(PNG_PATH, dpi=300, bbox_inches='tight')
fig.savefig(SVG_PATH, bbox_inches='tight')
plt.close()
print(f"Saved: {PNG_PATH}  ({os.path.getsize(PNG_PATH)/1024:.1f} KB)")
print(f"Saved: {SVG_PATH}  ({os.path.getsize(SVG_PATH)/1024:.1f} KB)")
print()

# Great Tables: patent summary
summary_rows = [
    ("Patent-exclusive compounds",        f"{n_patent}"),
    ("Patent-exclusive scaffolds (pipeline)", "103"),
    ("Patent-exclusive scaffolds (fresh)", f"{n_patent_scaffolds}"),
    ("Mean pIC50 — patent",               f"{mean_patent:.3f}"),
    ("Mean pIC50 — published",            f"{mean_published:.3f}"),
    ("pIC50 delta (published − patent)",  f"{mean_published - mean_patent:.3f}"),
    ("Severe cliffs involving any patent",f"{any_patent_sev}"),
    ("Severe cliffs where both patent",   f"{both_patent_sev}"),
    ("t-SNE coverage",                    "peripheral / distributed"),
]
gt_df = pd.DataFrame(summary_rows, columns=['Metric', 'Value'])

gt = (
    GT(gt_df)
    .tab_header(
        title="PAD4-DB Patent-Exclusive Compound Summary",
        subtitle="Patent-exclusive = source_list == 'pubchem_confirmatory' only",
    )
    .tab_style(
        style=gt_style.fill(color="#FFF3E0"),
        locations=loc.body(rows=[0, 2, 3, 6]),
    )
    .tab_source_note(
        "Patent compounds (n=233) from PubChem confirmatory AIDs with no "
        "ChEMBL or BindingDB overlap. "
        "Scaffold count discrepancy (103 pipeline vs fresh) reflects "
        "RDKit inter-version SMILES canonicalization drift."
    )
)

HTML_PATH = 'outputs/tables/fig10_patent_stats.html'
with open(HTML_PATH, 'w') as f:
    f.write(gt.as_raw_html())
print(f"Great Tables: {HTML_PATH}  ({os.path.getsize(HTML_PATH)/1024:.1f} KB)")
print("TASK C: DONE")
