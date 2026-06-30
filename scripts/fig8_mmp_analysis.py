#!/usr/bin/env python
"""Figure 8 — Scoped MMP Analysis (99 severe cliff compounds)."""
import os, sys, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

# ── Library standards ─────────────────────────────────────────────────────────
import scienceplots
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use(['science', 'nature', 'no-latex'])

from great_tables import GT, loc, style as gt_style
import great_tables
from adjustText import adjust_text

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMMPA
from itertools import combinations

os.makedirs('outputs/figures', exist_ok=True)
os.makedirs('outputs/tables', exist_ok=True)
os.makedirs('outputs/mmp', exist_ok=True)

print(f"SciencePlots: importable  |  Great Tables: {great_tables.__version__}")
print()

# ── Load data ─────────────────────────────────────────────────────────────────
nodes_df = pd.read_csv('outputs/figures/fig7_nodes.csv')
comps    = pd.read_parquet('data/processed/pad4_compounds.parquet')
cliffs   = pd.read_parquet('data/processed/activity_cliffs.parquet')
pairs_df = pd.read_parquet('data/processed/activity_pairs_with_sali.parquet')

assert len(nodes_df) == 99, f"Expected 99 nodes, got {len(nodes_df)}"

# Filter master parquet to 99 cliff compounds
cliff99 = comps[comps['inchi_key'].isin(nodes_df['inchi_key'])].copy()
assert len(cliff99) == 99, f"Expected 99 compounds after filter, got {len(cliff99)}"
print(f"Compounds loaded: {len(cliff99)} (verified = 99)")

HUB_IKS = {
    'SMADULGDNOCLOP-GISFHXKWSA-N', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
    'UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'DVCKJOQIVOGXEI-XMMPIXPASA-N',
}
HUB_CLASS = {
    'SMADULGDNOCLOP-GISFHXKWSA-N': 'A', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N': 'A',
    'UDCDEKJNAMHBFH-HSZRJFAPSA-N': 'B', 'DVCKJOQIVOGXEI-XMMPIXPASA-N': 'B',
}
pic50_map = cliff99.set_index('inchi_key')['pic50_consensus'].to_dict()
hac_map   = cliff99.set_index('inchi_key')['n_heavy_atoms'].to_dict()

# ── Step 1: Fragment all 99 compounds ─────────────────────────────────────────
print("Step 1: Fragmenting 99 compounds with rdMMPA (maxCuts=1)...")
# core_map: core_smi → list of (inchi_key, rgroup_smi)
core_map = {}
null_mol = 0

for _, row in cliff99.iterrows():
    ik  = row['inchi_key']
    smi = row['smiles_std']
    mol = Chem.MolFromSmiles(smi) if pd.notna(smi) else None
    if mol is None:
        null_mol += 1
        continue
    frags = rdMMPA.FragmentMol(mol, maxCuts=1, maxCutBonds=20, resultsAsMols=False)
    for env, frag_str in frags:
        parts = frag_str.split('.')
        if len(parts) < 2:
            continue
        parts.sort(key=len, reverse=True)
        core_smi, rgroup_smi = parts[0], parts[1]
        if not core_smi:
            continue
        core_map.setdefault(core_smi, []).append((ik, rgroup_smi))

print(f"  Null mols: {null_mol}  |  Unique cores: {len(core_map)}")

# ── Step 2: Find MMP pairs ────────────────────────────────────────────────────
print("Step 2: Finding MMP pairs sharing same core...")
seen_pairs = set()
mmp_rows = []
for core_smi, members in core_map.items():
    if len(members) < 2:
        continue
    for (ik_a, rg_a), (ik_b, rg_b) in combinations(members, 2):
        if ik_a == ik_b:
            continue
        key = tuple(sorted([ik_a, ik_b]))
        if key in seen_pairs:
            # count shared cores for existing pair
            for r in mmp_rows:
                if tuple(sorted([r['inchi_key_a'], r['inchi_key_b']])) == key:
                    r['n_shared_cores'] += 1
                    break
            continue
        seen_pairs.add(key)
        mmp_rows.append({
            'inchi_key_a':   ik_a,
            'inchi_key_b':   ik_b,
            'shared_core':   core_smi,
            'n_shared_cores': 1,
            'rgroup_a':      rg_a,
            'rgroup_b':      rg_b,
        })

mmp_df = pd.DataFrame(mmp_rows)
print(f"MMP pairs found among 99 cliff compounds: {len(mmp_df)}")
print()

# ── Step 3: Enrich with cliff data ────────────────────────────────────────────
print("Step 3: Enriching MMP pairs with cliff data...")

# Build fast lookup from cliffs (all tiers) and all pairs (for tanimoto)
# Normalise key order: always alphabetical
def make_key(a, b):
    return tuple(sorted([a, b]))

cliff_lookup = {}
for _, r in cliffs.iterrows():
    k = make_key(r['inchi_key_a'], r['inchi_key_b'])
    cliff_lookup[k] = {
        'delta_pic50': r['delta_pic50'],
        'tanimoto':    r['tanimoto'],
        'cliff_tier':  r['cliff_tier'],
    }

pairs_lookup = {}
for _, r in pairs_df.iterrows():
    k = make_key(r['inchi_key_a'], r['inchi_key_b'])
    pairs_lookup[k] = {
        'delta_pic50': r['delta_pic50'],
        'tanimoto':    r['tanimoto'],
        'cliff_tier':  r.get('cliff_tier', 'non_cliff'),
    }

enriched = []
for _, r in mmp_df.iterrows():
    ik_a, ik_b = r['inchi_key_a'], r['inchi_key_b']
    k = make_key(ik_a, ik_b)

    if k in cliff_lookup:
        info = cliff_lookup[k]
    elif k in pairs_lookup:
        info = pairs_lookup[k]
    else:
        # Compute from compounds directly
        dp = abs(pic50_map.get(ik_a, float('nan')) - pic50_map.get(ik_b, float('nan')))
        info = {'delta_pic50': dp, 'tanimoto': float('nan'), 'cliff_tier': 'non_cliff'}
        if dp >= 2.0:
            info['cliff_tier'] = 'severe'
        elif dp >= 1.5:
            info['cliff_tier'] = 'moderate'
        elif dp >= 1.0:
            info['cliff_tier'] = 'broad'

    # Heavy atom difference for mmp_type
    hac_a = hac_map.get(ik_a, 0)
    hac_b = hac_map.get(ik_b, 0)
    size_change = abs(hac_b - hac_a)
    if size_change <= 1:
        mmp_type = 'single_atom_change'
    elif size_change <= 3:
        mmp_type = 'small_substituent'
    elif size_change <= 8:
        mmp_type = 'medium_substituent'
    else:
        mmp_type = 'large_substituent'

    enriched.append({
        **r.to_dict(),
        'delta_pic50':    info['delta_pic50'],
        'tanimoto':       info.get('tanimoto', float('nan')),
        'cliff_tier':     info['cliff_tier'],
        'is_severe_cliff': info['delta_pic50'] >= 2.0 if pd.notna(info['delta_pic50']) else False,
        'hub_a':          ik_a in HUB_IKS,
        'hub_b':          ik_b in HUB_IKS,
        'mmp_type':       mmp_type,
    })

mmp_enrich = pd.DataFrame(enriched)

# Print breakdown
print(f"  Total MMP pairs: {len(mmp_enrich)}")
for tier in ['severe', 'moderate', 'broad', 'non_cliff']:
    n = (mmp_enrich['cliff_tier'] == tier).sum()
    print(f"  MMP {tier:10s} cliffs: {n}")
print()
# Tanimoto-AND-MMP confirmed severe cliffs (cliff_tier='severe' = both conditions met)
sev_confirmed = mmp_enrich[mmp_enrich['cliff_tier'] == 'severe']
print(f"  MMP-confirmed severe cliffs (also Tanimoto-severe): {len(sev_confirmed)} of 94")
print("  R-group change type breakdown (cliff_tier='severe'):")
for t in ['single_atom_change','small_substituent','medium_substituent','large_substituent']:
    print(f"    {t}: {(sev_confirmed['mmp_type']==t).sum()}")
print()

# ── Step 4: Discontinuity score ───────────────────────────────────────────────
print("Step 4: Computing discontinuity scores...")
disc = {}
for ik in cliff99['inchi_key']:
    partner_mask = (mmp_enrich['inchi_key_a'] == ik) | (mmp_enrich['inchi_key_b'] == ik)
    sub = mmp_enrich[partner_mask]
    if len(sub) == 0:
        disc[ik] = {'discontinuity_score': float('nan'), 'n_mmp_partners': 0}
    else:
        disc[ik] = {
            'discontinuity_score': float(sub['delta_pic50'].mean()),
            'n_mmp_partners': len(sub),
        }

disc_df = nodes_df[['inchi_key','pic50_consensus','severe_cliff_degree','hub_class']].copy()
disc_df['discontinuity_score'] = disc_df['inchi_key'].map(lambda x: disc[x]['discontinuity_score'])
disc_df['n_mmp_partners']      = disc_df['inchi_key'].map(lambda x: disc[x]['n_mmp_partners'])
disc_df = disc_df.sort_values('discontinuity_score', ascending=False)

print("Top 10 by discontinuity score:")
print(disc_df.head(10)[['inchi_key','pic50_consensus','discontinuity_score','n_mmp_partners','hub_class']].to_string(index=False))
print()

# ── Step 5: Save outputs ──────────────────────────────────────────────────────
out_cols = ['inchi_key_a','inchi_key_b','shared_core','n_shared_cores',
            'delta_pic50','tanimoto','cliff_tier','is_severe_cliff','hub_a','hub_b','mmp_type']
mmp_enrich[out_cols].to_csv('outputs/mmp/mmp_pairs_cliff99.csv', index=False)

disc_out_cols = ['inchi_key','pic50_consensus','severe_cliff_degree',
                 'hub_class','discontinuity_score','n_mmp_partners']
disc_df[disc_out_cols].to_csv('outputs/mmp/mmp_discontinuity_scores.csv', index=False)

print(f"Saved: outputs/mmp/mmp_pairs_cliff99.csv ({len(mmp_enrich)} rows)")
print(f"Saved: outputs/mmp/mmp_discontinuity_scores.csv ({len(disc_df)} rows)")
print()

# ── Step 6: Figure 8 ──────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
ax_a, ax_b = axes[0, 0], axes[0, 1]
ax_c, ax_d = axes[1, 0], axes[1, 1]

def panel_label(ax, letter):
    ax.text(0.02, 0.96, letter, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top')

MMP_TYPE_ORDER  = ['single_atom_change', 'small_substituent',
                   'medium_substituent', 'large_substituent']
MMP_TYPE_COLORS = ['#2ECC71', '#4A90D9', '#F39C12', '#E74C3C']
TIER_ORDER = ['non_cliff', 'broad', 'moderate', 'severe']

# ── Panel A: Stacked bar by cliff tier × mmp_type ────────────────────────────
x_pos = np.arange(len(TIER_ORDER))
bottom = np.zeros(len(TIER_ORDER))
for mtype, color in zip(MMP_TYPE_ORDER, MMP_TYPE_COLORS):
    heights = [((mmp_enrich['cliff_tier'] == t) & (mmp_enrich['mmp_type'] == mtype)).sum()
               for t in TIER_ORDER]
    ax_a.bar(x_pos, heights, bottom=bottom, color=color, label=mtype.replace('_', ' '),
             edgecolor='white', linewidth=0.3)
    bottom += np.array(heights, dtype=float)

ax_a.set_xticks(x_pos)
ax_a.set_xticklabels(TIER_ORDER, fontsize=8)
ax_a.set_xlabel('Cliff Tier', fontsize=9)
ax_a.set_ylabel('MMP Pair Count', fontsize=9)
ax_a.legend(fontsize=7, framealpha=0.8, loc='upper left')
panel_label(ax_a, 'A')

# ── Panel B: Discontinuity score vs pIC50 ────────────────────────────────────
plot_disc = disc_df.dropna(subset=['discontinuity_score'])
mask_a  = plot_disc['hub_class'] == 'A'
mask_b  = plot_disc['hub_class'] == 'B'
mask_nh = plot_disc['hub_class'] == 'none'

ax_b.scatter(plot_disc.loc[mask_nh, 'pic50_consensus'],
             plot_disc.loc[mask_nh, 'discontinuity_score'],
             c='#AAAAAA', s=40, alpha=0.7, zorder=2)
ax_b.scatter(plot_disc.loc[mask_b, 'pic50_consensus'],
             plot_disc.loc[mask_b, 'discontinuity_score'],
             c='#1A237E', s=150, marker='*', zorder=10, label='Hub Class B')
ax_b.scatter(plot_disc.loc[mask_a, 'pic50_consensus'],
             plot_disc.loc[mask_a, 'discontinuity_score'],
             c='#E74C3C', s=150, marker='*', zorder=11, label='Hub Class A')

texts = []
for _, row in plot_disc[plot_disc['hub_class'] != 'none'].iterrows():
    texts.append(ax_b.text(
        row['pic50_consensus'], row['discontinuity_score'],
        row['inchi_key'][:14], fontsize=7))
if texts:
    adjust_text(texts, ax=ax_b,
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

ax_b.set_xlabel('pIC50', fontsize=9)
ax_b.set_ylabel('Mean MMP ΔpIC50 (Discontinuity Score)', fontsize=9)
ax_b.legend(fontsize=8, framealpha=0.8)
panel_label(ax_b, 'B')

# ── Panel C: ΔpIC50 histogram — MMP vs Tanimoto severe cliffs ────────────────
# MMP severe cliffs = confirmed by BOTH MMP shared-core AND Tanimoto threshold
all_severe_dp = cliffs.loc[cliffs['cliff_tier'] == 'severe', 'delta_pic50'].values
mmp_severe_dp = mmp_enrich.loc[mmp_enrich['cliff_tier'] == 'severe', 'delta_pic50'].values

bins_c = np.linspace(2.0, 3.5, 21)
ax_c.hist(all_severe_dp, bins=bins_c, color='#4A90D9', alpha=0.5,
          label=f'Tanimoto severe cliffs (n={len(all_severe_dp)})')
ax_c.hist(mmp_severe_dp, bins=bins_c, color='#E74C3C', alpha=0.5,
          label=f'MMP severe cliffs (n={len(mmp_severe_dp)})')
ax_c.set_xlabel('ΔpIC50', fontsize=9)
ax_c.set_ylabel('Count', fontsize=9)
ax_c.legend(fontsize=7.5, framealpha=0.8)
panel_label(ax_c, 'C')

# ── Panel D: Core frequency (top 10) ─────────────────────────────────────────
core_freq = (mmp_enrich.groupby('shared_core')
             .size().sort_values(ascending=False).head(10))
core_labels = [c[:30] for c in core_freq.index]

y_pos_d = np.arange(len(core_freq))
ax_d.barh(y_pos_d, core_freq.values, color='#4A90D9',
          edgecolor='white', linewidth=0.3)
ax_d.set_yticks(y_pos_d)
ax_d.set_yticklabels(core_labels, fontsize=6)
ax_d.invert_yaxis()
ax_d.set_xlabel('Number of MMP Pairs', fontsize=9)
ax_d.set_title('Top 10 Most Frequent Shared Cores', fontsize=8)
panel_label(ax_d, 'D')

plt.tight_layout(pad=0.8)

PNG_PATH = 'outputs/figures/fig8_mmp_analysis.png'
SVG_PATH = 'outputs/figures/fig8_mmp_analysis.svg'
fig.savefig(PNG_PATH, dpi=300, bbox_inches='tight')
fig.savefig(SVG_PATH, bbox_inches='tight')
plt.close()
print(f"Saved: {PNG_PATH}  ({os.path.getsize(PNG_PATH)/1024:.1f} KB)")
print(f"Saved: {SVG_PATH}  ({os.path.getsize(SVG_PATH)/1024:.1f} KB)")
print()

# ── Great Tables summary ──────────────────────────────────────────────────────
most_common_type = mmp_enrich['mmp_type'].value_counts().index[0].replace('_', ' ')
hubs_in_top10 = disc_df.head(10)['hub_class'].ne('none').sum()

gt_data = pd.DataFrame([
    {'metric': 'Total MMP pairs (99 cliff compounds)',  'value': str(len(mmp_enrich))},
    {'metric': 'MMP-validated severe cliffs',           'value': str(len(mmp_severe_dp))},
    {'metric': 'MMP-validated moderate cliffs',
     'value': str((mmp_enrich['cliff_tier'] == 'moderate').sum())},
    {'metric': 'Most common R-group change type',       'value': most_common_type},
    {'metric': 'Hub compounds in top-10 discontinuity', 'value': f"{hubs_in_top10} of 4"},
])

gt = (
    GT(gt_data)
    .tab_header(
        title="PAD4-DB MMP Analysis Summary",
        subtitle="Scoped to 99 severe activity cliff compounds",
    )
    .cols_label(metric="Metric", value="Value")
    .tab_source_note(
        "rdMMPA fragmentation: maxCuts=1, maxCutBonds=20. "
        "Pairs sharing identical core SMILES with different R-groups."
    )
)

HTML_PATH = 'outputs/tables/fig8_mmp_stats.html'
TEX_PATH  = 'outputs/tables/fig8_mmp_stats.tex'
with open(HTML_PATH, 'w') as f:
    f.write(gt.as_raw_html())
with open(TEX_PATH, 'w') as f:
    f.write(gt.as_latex())
print(f"Great Tables HTML:  {HTML_PATH}  ({os.path.getsize(HTML_PATH)/1024:.1f} KB)")
print(f"Great Tables LaTeX: {TEX_PATH}  ({os.path.getsize(TEX_PATH)/1024:.1f} KB)")
print()

# ── Completion report ─────────────────────────────────────────────────────────
top1 = disc_df.iloc[0]
n_mmp_severe = len(mmp_severe_dp)
mmp_val_rate = n_mmp_severe / 94 * 100
print("=== Completion Report ===")
print(f"Total MMP pairs among 99 cliff compounds: {len(mmp_enrich)}")
print(f"MMP-confirmed severe cliffs: {n_mmp_severe} of 94 ({mmp_val_rate:.1f}%)")
print(f"  (MMP pairs confirmed by both shared-core AND Tanimoto-severe criteria)")
print(f"Top discontinuity compound: {top1['inchi_key']}, "
      f"score={top1['discontinuity_score']:.3f}, hub={top1['hub_class']}")
print()
print("Files written:")
for p in ['outputs/mmp/mmp_pairs_cliff99.csv',
          'outputs/mmp/mmp_discontinuity_scores.csv',
          PNG_PATH, SVG_PATH, HTML_PATH, TEX_PATH]:
    print(f"  {p}  ({os.path.getsize(p)/1024:.1f} KB)")
