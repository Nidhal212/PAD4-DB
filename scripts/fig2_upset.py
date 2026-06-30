#!/usr/bin/env python
import os
import sys
import warnings
warnings.filterwarnings('ignore')

# Ensure working from project root
os.chdir('/home/nidhal/PAD4-db_V2')

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from upsetplot import UpSet, from_memberships

os.makedirs('outputs/figures', exist_ok=True)

# Font setup
available = {f.name for f in fm.fontManager.ttflist}
font_family = 'Arial' if 'Arial' in available else 'DejaVu Sans'
plt.rcParams.update({
    'font.family': font_family,
    'font.size': 10,
    'axes.titlesize': 10,
    'axes.labelsize': 10,
})

# Load data
df = pd.read_parquet('data/processed/pad4_compounds.parquet')
SOURCE_MAP = {
    'pubchem_confirmatory': 'PubChem',
    'chembl': 'ChEMBL',
    'bindingdb': 'BindingDB',
}

memberships = []
for sl in df['source_list']:
    memberships.append([SOURCE_MAP[p] for p in sl.split('|') if p in SOURCE_MAP])
data = from_memberships(memberships)

# Verify locked counts
LOCKED = {
    (False, False, True): 233,
    (True,  False, False): 95,
    (False, True,  False): 10,
    (True,  False, True): 1199,
    (True,  True,  False): 167,
    (False, True,  True): 23,
    (True,  True,  True): 1366,
}
counts_agg = data.groupby(level=['BindingDB', 'ChEMBL', 'PubChem']).count()
all_pass = True
for k, v in LOCKED.items():
    actual = counts_agg.get(k, 0)
    status = 'PASS' if actual == v else f'FAIL (got {actual})'
    print(f'{k}: expected={v}, {status}')
    if actual != v:
        all_pass = False
total = len(data)
print(f'Total: {total} (expected 3093) {"PASS" if total == 3093 else "FAIL"}')
if not all_pass or total != 3093:
    sys.exit(1)
print()

# Build figure
fig = plt.figure(figsize=(8, 5))
upset = UpSet(
    data,
    sort_by='degree',
    sort_categories_by=None,
    show_counts=True,
    totals_plot_elements=0,
    subset_size='count',
)
# PubChem-only bar colored orange
upset.style_subsets(
    present='PubChem',
    absent=['BindingDB', 'ChEMBL'],
    facecolor='#E05A2B',
    label='_pubchem_only',
)
dict_axes = upset.plot(fig)

ax_bars = dict_axes['intersections']

# Style bars: grey for all, orange retained for PubChem-only
for bar in ax_bars.patches:
    r, g, b, _ = bar.get_facecolor()
    is_orange = (r > 0.8 and g < 0.5 and b < 0.3)
    if not is_orange:
        bar.set_facecolor('#AAAAAA')
    bar.set_edgecolor('none')

# Remove default count text, add custom labels
for txt in ax_bars.texts:
    txt.set_visible(False)

for bar in ax_bars.patches:
    x = bar.get_x() + bar.get_width() / 2
    h = bar.get_height()
    r, g, b, _ = bar.get_facecolor()
    is_orange = (r > 0.8 and g < 0.5 and b < 0.3)
    if is_orange:
        ax_bars.text(
            x, h + 20,
            "233\n(patent-\nexclusive)",
            ha='center', va='bottom',
            fontsize=8, fontweight='bold',
            color='#E05A2B',
        )
    else:
        ax_bars.text(
            x, h + 15,
            str(int(h)),
            ha='center', va='bottom',
            fontsize=8, color='#333333',
        )

ax_bars.set_ylabel('Compounds', fontsize=10)
ax_bars.spines['top'].set_visible(False)
ax_bars.spines['right'].set_visible(False)

# Footnote
footnote = (
    "Source abbreviations: PubChem = PubChem Confirmatory Assay data (includes patent-derived records); "
    "ChEMBL = ChEMBL v34; BindingDB = BindingDB (UniProt Q9UM07).\n"
    "Counts reflect deduplicated compounds (unique InChIKey). "
    "Intersection membership = measured in ≥1 assay from that source."
)
fig.text(
    0.5, -0.02,
    footnote,
    ha='center', va='top',
    fontsize=8, color='#666666',
    transform=fig.transFigure,
)

fig.suptitle(
    'Figure 2. Source overlap of PAD4 inhibitor database (n = 3,093)',
    fontsize=10, fontweight='bold', y=1.01,
)

plt.tight_layout(rect=[0, 0.08, 1, 1])

png_path = 'outputs/figures/fig2_source_overlap_upset.png'
svg_path = 'outputs/figures/fig2_source_overlap_upset.svg'
fig.savefig(png_path, dpi=300, bbox_inches='tight')
fig.savefig(svg_path, format='svg', bbox_inches='tight')
plt.close()

print(f'Saved: {png_path}  ({os.path.getsize(png_path):,} bytes)')
print(f'Saved: {svg_path}  ({os.path.getsize(svg_path):,} bytes)')
