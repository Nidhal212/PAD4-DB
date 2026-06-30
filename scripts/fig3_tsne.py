#!/usr/bin/env python
"""Figure 3 — t-SNE Chemical Space (4-panel) + Great Tables summary."""
import os, sys, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

# ── Step 0: Verify libraries ─────────────────────────────────────────────────
import scienceplots
from great_tables import GT, loc, style as gt_style
import great_tables
print(f"SciencePlots: importable (version attribute unavailable in this build)")
print(f"Great Tables: {great_tables.__version__}")
print()

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.manifold import TSNE

os.makedirs('outputs/figures', exist_ok=True)
os.makedirs('outputs/tables', exist_ok=True)

# ── Step 1: Morgan fingerprints ──────────────────────────────────────────────
FP_CACHE = 'data/interim/morgan_fps_3093.npy'
df = pd.read_parquet('data/processed/pad4_compounds.parquet')
assert len(df) == 3093, f"Expected 3093 compounds, got {len(df)}"

if os.path.exists(FP_CACHE):
    print(f"Loading fingerprints from cache: {FP_CACHE}")
    fps_matrix = np.load(FP_CACHE)
else:
    print("Computing Morgan fingerprints (ECFP4, radius=2, nBits=2048)...")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps_list = []
    for smi in df['smiles_std']:
        mol = Chem.MolFromSmiles(smi) if pd.notna(smi) else None
        if mol is None:
            fps_list.append(np.zeros(2048, dtype=np.uint8))
        else:
            arr = np.zeros(2048, dtype=np.uint8)
            fp = gen.GetFingerprint(mol)
            for bit in fp.GetOnBits():
                arr[bit] = 1
            fps_list.append(arr)
    fps_matrix = np.vstack(fps_list).astype(np.uint8)
    np.save(FP_CACHE, fps_matrix)

print(f"Fingerprints computed: {fps_matrix.shape[0]} × {fps_matrix.shape[1]}")
print()

# ── Step 2: t-SNE ────────────────────────────────────────────────────────────
TSNE_CACHE = 'data/interim/tsne_coords_3093.npy'

if os.path.exists(TSNE_CACHE):
    print(f"Loading t-SNE coordinates from cache: {TSNE_CACHE}")
    coords = np.load(TSNE_CACHE)
else:
    print("Running t-SNE — this may take up to 15 minutes...")
    tsne = TSNE(
        n_components=2,
        perplexity=40,
        max_iter=1000,
        random_state=42,
        metric='jaccard',
        init='pca',
    )
    coords = tsne.fit_transform(fps_matrix.astype(float))
    np.save(TSNE_CACHE, coords)

x, y = coords[:, 0], coords[:, 1]
print(f"t-SNE complete. X range: [{x.min():.2f}, {x.max():.2f}]  "
      f"Y range: [{y.min():.2f}, {y.max():.2f}]")
print()

# ── Step 3: Four-panel figure ─────────────────────────────────────────────────
plt.style.use(['science', 'nature', 'no-latex'])
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
ax_a, ax_b = axes[0, 0], axes[0, 1]
ax_c, ax_d = axes[1, 0], axes[1, 1]

def clean_ax(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('t-SNE 1', fontsize=9)
    ax.set_ylabel('t-SNE 2', fontsize=9)

def panel_label(ax, letter):
    ax.text(0.02, 0.96, letter, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top')

# ── Panel A: Source category ─────────────────────────────────────────────────
SOURCE_CATS = {
    'pubchem_confirmatory': ('PubChem only (patent)', '#E05A2B'),
    'bindingdb':            ('BindingDB only',        '#4A90D9'),
    'chembl':               ('ChEMBL only',           '#2ECC71'),
}

def get_source_cat(sl):
    if sl == 'pubchem_confirmatory':
        return 'PubChem only (patent)'
    elif sl == 'bindingdb':
        return 'BindingDB only'
    elif sl == 'chembl':
        return 'ChEMBL only'
    else:
        return 'Multi-source'

df['_src_cat'] = df['source_list'].map(get_source_cat)

# Plot order: Multi-source first (background), singles on top
plot_order = [
    ('Multi-source',        '#AAAAAA'),
    ('PubChem only (patent)', '#E05A2B'),
    ('BindingDB only',       '#4A90D9'),
    ('ChEMBL only',          '#2ECC71'),
]

panel_a_counts = {}
for cat, color in plot_order:
    mask = df['_src_cat'] == cat
    panel_a_counts[cat] = int(mask.sum())
    ax_a.scatter(x[mask], y[mask], c=color, s=8, alpha=0.6,
                 label=f'{cat} (n={mask.sum()})', rasterized=True)

ax_a.legend(loc='upper right', framealpha=0.8, fontsize=7, markerscale=1.5)
clean_ax(ax_a)
panel_label(ax_a, 'A')

# ── Panel B: pIC50 continuous ─────────────────────────────────────────────────
sc_b = ax_b.scatter(x, y, c=df['pic50_consensus'], cmap='viridis',
                    vmin=2.0, vmax=8.52, s=8, alpha=0.6, rasterized=True)
cb = plt.colorbar(sc_b, ax=ax_b, pad=0.02)
cb.set_label('pIC50', fontsize=9)
cb.ax.tick_params(labelsize=8)
clean_ax(ax_b)
panel_label(ax_b, 'B')

# ── Panel C: mechanism_class ─────────────────────────────────────────────────
MECH_COLORS = {
    'enzymatic':           '#AAAAAA',
    'enzymatic_confirmed': '#4A90D9',
    'fp_ic50':             '#2ECC71',
    'covalent':            '#E05A2B',
}
mech_order = ['enzymatic', 'enzymatic_confirmed', 'fp_ic50', 'covalent']

for mech in mech_order:
    mask = df['mechanism_class'] == mech
    ax_c.scatter(x[mask], y[mask], c=MECH_COLORS[mech], s=8, alpha=0.7,
                 label=f'{mech} (n={mask.sum()})', rasterized=True)

ax_c.legend(loc='upper right', framealpha=0.8, fontsize=7, markerscale=1.5)
clean_ax(ax_c)
panel_label(ax_c, 'C')

# ── Panel D: Hub highlight ────────────────────────────────────────────────────
HUB_A = {'SMADULGDNOCLOP-GISFHXKWSA-N', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N'}
HUB_B = {'UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'DVCKJOQIVOGXEI-XMMPIXPASA-N'}

severe_csv = pd.read_csv('outputs/figures/fig7_nodes.csv')
severe_iks = set(severe_csv['inchi_key'])

mask_hub_a = df['inchi_key'].isin(HUB_A)
mask_hub_b = df['inchi_key'].isin(HUB_B)
mask_severe_nonhub = df['inchi_key'].isin(severe_iks - HUB_A - HUB_B)

panel_d_counts = {
    'Hub Class A': int(mask_hub_a.sum()),
    'Hub Class B': int(mask_hub_b.sum()),
    'Severe cliff non-hub': int(mask_severe_nonhub.sum()),
    'Background': len(df),
}

# Layer 1: all compounds (background)
ax_d.scatter(x, y, c='#DDDDDD', s=6, alpha=0.4, rasterized=True, zorder=1)
# Layer 2: severe cliff non-hubs
ax_d.scatter(x[mask_severe_nonhub], y[mask_severe_nonhub],
             c='#F39C12', s=20, alpha=0.8, zorder=2,
             label=f'Severe cliff compound (n={mask_severe_nonhub.sum()})')
# Layer 3: Class B hubs
ax_d.scatter(x[mask_hub_b], y[mask_hub_b],
             c='#1A237E', s=150, marker='*', alpha=1.0, zorder=10,
             label='Hub Class B — singleton (n=2)')
# Layer 4: Class A hubs
ax_d.scatter(x[mask_hub_a], y[mask_hub_a],
             c='#E74C3C', s=150, marker='*', alpha=1.0, zorder=11,
             label='Hub Class A — series-embedded (n=2)')

ax_d.legend(loc='upper right', framealpha=0.8, fontsize=7, markerscale=1.2)
clean_ax(ax_d)
panel_label(ax_d, 'D')

plt.tight_layout(pad=0.8)

PNG_PATH = 'outputs/figures/fig3_tsne_chemical_space.png'
SVG_PATH = 'outputs/figures/fig3_tsne_chemical_space.svg'
fig.savefig(PNG_PATH, dpi=300, bbox_inches='tight')
fig.savefig(SVG_PATH, bbox_inches='tight')
plt.close()
print(f"Saved: {PNG_PATH}  ({os.path.getsize(PNG_PATH)/1024:.1f} KB)")
print(f"Saved: {SVG_PATH}  ({os.path.getsize(SVG_PATH)/1024:.1f} KB)")
print()

# ── Step 4: Great Tables summary ─────────────────────────────────────────────
summary_rows = []
for mech in mech_order:
    mask = df['mechanism_class'] == mech
    idx = np.where(mask.values)[0]
    summary_rows.append({
        'mechanism_class': mech,
        'n_compounds':     int(mask.sum()),
        'mean_pic50':      float(df.loc[mask, 'pic50_consensus'].mean()),
        'std_pic50':       float(df.loc[mask, 'pic50_consensus'].std()),
        'mean_tsne1':      float(x[idx].mean()),
        'mean_tsne2':      float(y[idx].mean()),
    })
summary_df = pd.DataFrame(summary_rows)

gt = (
    GT(summary_df)
    .tab_header(
        title="PAD4-DB Chemical Space Summary",
        subtitle="t-SNE embedding of 3,093 compounds by ECFP4 fingerprints",
    )
    .cols_label(
        mechanism_class="Mechanism Class",
        n_compounds="N",
        mean_pic50="Mean pIC50",
        std_pic50="SD pIC50",
        mean_tsne1="t-SNE 1 (mean)",
        mean_tsne2="t-SNE 2 (mean)",
    )
    .fmt_number(
        columns=["mean_pic50", "std_pic50", "mean_tsne1", "mean_tsne2"],
        decimals=2,
    )
    .tab_style(
        style=gt_style.fill(color="#FFF3E0"),
        locations=loc.body(rows=[0]),
    )
    .tab_source_note(
        "ECFP4 fingerprints (radius=2, 2048 bits); "
        "t-SNE: perplexity=40, random_state=42, metric=jaccard"
    )
)

HTML_PATH = 'outputs/tables/fig3_tsne_summary.html'
with open(HTML_PATH, 'w') as f:
    f.write(gt.as_raw_html())
print(f"Great Tables HTML written: {HTML_PATH}  ({os.path.getsize(HTML_PATH)/1024:.1f} KB)")
print()

# ── Step 5: Verification ─────────────────────────────────────────────────────
print("=" * 60)
print("VERIFICATION REPORT")
print("=" * 60)
print(f"Library check:")
print(f"  SciencePlots: importable")
print(f"  Great Tables: {great_tables.__version__}")
print()
print(f"Data:")
print(f"  Compounds plotted:   {len(df):,} (all panels)")
print(f"  Morgan FP shape:     {fps_matrix.shape}")
print(f"  t-SNE coords shape:  {coords.shape}")
print()
print(f"Panel A counts:")
for cat, n in panel_a_counts.items():
    print(f"  {cat}: {n}")
total_a = sum(panel_a_counts.values())
print(f"  TOTAL: {total_a}  {'PASS' if total_a == 3093 else 'FAIL'}")

EXPECTED_A = {
    'PubChem only (patent)': 233,
    'BindingDB only': 95,
    'ChEMBL only': 10,
    'Multi-source': 2755,
}
a_pass = all(panel_a_counts.get(k) == v for k, v in EXPECTED_A.items())
print(f"  Locked count match:  {'PASS' if a_pass else 'FAIL'}")

print()
print(f"Panel D counts:")
for k, v in panel_d_counts.items():
    print(f"  {k}: {v}")
d_hub_a_ok  = panel_d_counts['Hub Class A'] == 2
d_hub_b_ok  = panel_d_counts['Hub Class B'] == 2
d_severe_ok = panel_d_counts['Severe cliff non-hub'] == 95
d_bg_ok     = panel_d_counts['Background'] == 3093
print(f"  Hub A=2:  {'PASS' if d_hub_a_ok else 'FAIL'}")
print(f"  Hub B=2:  {'PASS' if d_hub_b_ok else 'FAIL'}")
print(f"  Severe non-hub=95: {'PASS' if d_severe_ok else 'FAIL'}")
print(f"  Background=3093:   {'PASS' if d_bg_ok else 'FAIL'}")

print()
print(f"Files written:")
for p in [PNG_PATH, SVG_PATH, HTML_PATH]:
    kb = os.path.getsize(p) / 1024
    print(f"  {p}  ({kb:.1f} KB)")

all_pass = a_pass and d_hub_a_ok and d_hub_b_ok and d_severe_ok and d_bg_ok and total_a == 3093
print()
print("OVERALL:", "PASS" if all_pass else "FAIL — review above")
if not all_pass:
    sys.exit(1)
