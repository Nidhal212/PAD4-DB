#!/usr/bin/env python
"""Nature Tables 1–9 — Great Tables redesign."""
import os, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import pandas as pd
import numpy as np
from great_tables import GT, loc, style, html

OUT = 'outputs/tables/nature'
os.makedirs(OUT, exist_ok=True)

# ── Shared styling helper ─────────────────────────────────────────────────────
def nature_gt(df, title, subtitle='', source_note=None):
    """Apply Nature-standard GT styling: no vertical lines, clean header."""
    gt = GT(df)
    if title:
        gt = gt.tab_header(title=title, subtitle=subtitle or None)
    gt = (gt
        .tab_style(
            style=style.borders(sides='top',    weight='1.5px', color='#000000'),
            locations=loc.column_labels())
        .tab_style(
            style=style.borders(sides='bottom', weight='1px',   color='#000000'),
            locations=loc.column_labels())
        .tab_style(
            style=style.borders(sides='bottom', weight='1px',   color='#000000'),
            locations=loc.body(rows=[-1]))
        .tab_style(
            style=style.text(weight='bold', size='9px'),
            locations=loc.column_labels())
        .tab_style(
            style=style.text(size='8px'),
            locations=loc.body())
    )
    if source_note:
        gt = gt.tab_source_note(source_note)
    return gt

def save_gt(gt, name):
    html_path = f'{OUT}/{name}.html'
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(gt.as_raw_html())
    sz = os.path.getsize(html_path) // 1024
    print(f'  {name}.html  {sz} KB')

# ── Load core data once ────────────────────────────────────────────────────────
df     = pd.read_parquet('data/processed/pad4_compounds.parquet')
cliffs = pd.read_parquet('data/processed/activity_cliffs.parquet')
pairs  = pd.read_parquet('data/processed/activity_pairs_with_sali.parquet')
mmp    = pd.read_csv('outputs/mmp/mmp_pairs_cliff99.csv')
disc   = pd.read_csv('outputs/mmp/mmp_discontinuity_scores.csv')

print("=== Building Nature Tables ===")

# ── TABLE 1: Source Database Overview ─────────────────────────────────────────
print("\nTable 1: Source database coverage")
def src_stats(mask, name):
    sub = df.loc[mask, 'pic50_consensus']
    n   = mask.sum()
    return {
        'Source': name,
        'N compounds': f'{n:,}',
        '% Dataset': f'{n/3093*100:.1f}%',
        'Mean pIC50': round(float(sub.mean()), 2),
        'Median pIC50': round(float(sub.median()), 2),
        'SD pIC50': round(float(sub.std()), 2),
    }

sl = df['source_list']
rows_t1 = [
    src_stats(sl.str.contains('pubchem_confirmatory'), 'PubChem (≥1 AID)'),
    src_stats(sl.str.contains('chembl'),               'ChEMBL'),
    src_stats(sl.str.contains('bindingdb'),            'BindingDB'),
    src_stats(sl == 'bindingdb|chembl|pubchem_confirmatory', 'All three sources'),
    src_stats(sl == 'pubchem_confirmatory',            'PubChem only'),
]
t1 = pd.DataFrame(rows_t1)
gt1 = (nature_gt(t1, 'Source database coverage', 'PAD4-DB v2 (N=3,093)',
                  source_note='Values computed from deduplicated master dataset.')
       .cols_align(align='right',
                   columns=['N compounds','% Dataset','Mean pIC50','Median pIC50','SD pIC50'])
       .cols_align(align='left', columns=['Source']))
save_gt(gt1, 'table1_source_coverage')

# ── TABLE 2: pIC50 by Mechanism ───────────────────────────────────────────────
print("Table 2: pIC50 by assay mechanism")
mech_map = {
    'enzymatic':           'Enzymatic (BAEE colorimetric)',
    'enzymatic_confirmed': 'Enzymatic, RFMS-confirmed',
    'fp_ic50':             'FP-based IC50',
    'covalent':            'Covalent (assay-flagged)',
}
rows_t2 = []
for mech, label in mech_map.items():
    sub = df.loc[df['mechanism_class']==mech, 'pic50_consensus']
    rows_t2.append({
        'Assay mechanism': label,
        'N': f'{len(sub):,}',
        'Mean pIC50': round(float(sub.mean()), 2),
        'Median pIC50': round(float(sub.median()), 2),
        'SD': round(float(sub.std()), 2),
        'Min pIC50': round(float(sub.min()), 2),
        'Max pIC50': round(float(sub.max()), 2),
    })
t2 = pd.DataFrame(rows_t2)
gt2 = (nature_gt(t2, 'pIC50 by assay mechanism class',
                  source_note='pIC50 from consensus of all measurements per compound.')
       .cols_align(align='right',
                   columns=['N','Mean pIC50','Median pIC50','SD','Min pIC50','Max pIC50'])
       .cols_align(align='left', columns=['Assay mechanism']))
save_gt(gt2, 'table2_mechanism')

# ── TABLE 3: Top 20 Scaffold Series ──────────────────────────────────────────
print("Table 3: Top 20 scaffold series")
try:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    def get_sc(smi):
        try:
            mol = Chem.MolFromSmiles(smi)
            s   = MurckoScaffold.GetScaffoldForMol(mol)
            return Chem.MolToSmiles(s)
        except: return '__FAIL__'
    df['_scaffold'] = df['smiles_std'].map(get_sc)
    pat_mask = df['source_list'] == 'pubchem_confirmatory'
    df['_is_patent'] = pat_mask.values
    sc_grp = df[df['_scaffold'] != '__FAIL__'].groupby('_scaffold')
    sc_tab  = sc_grp.agg(
        n_compounds=('inchi_key','count'),
        mean_pic50=('pic50_consensus','mean'),
        pct_patent=('_is_patent','mean'),
    ).reset_index().sort_values('n_compounds', ascending=False).head(20)

    # Hub scaffold = rank 1 (n=174 azaindole-benzimidazole series)
    HUB_A_IK = {'SMADULGDNOCLOP-GISFHXKWSA-N','RAVBZQAQTVGKIV-XBPDSQQVSA-N'}
    hub_sc = df.loc[df['inchi_key'].isin(HUB_A_IK), '_scaffold'].iloc[0]
    sc_tab['hub_scaffold'] = sc_tab['_scaffold'].apply(lambda s: '✓' if s == hub_sc else '')
    sc_tab['pct_patent'] = (sc_tab['pct_patent']*100).round(1).astype(str) + '%'
    sc_tab['mean_pic50'] = sc_tab['mean_pic50'].round(2)
    sc_tab['Scaffold SMILES'] = sc_tab['_scaffold'].apply(
        lambda s: s[:35]+'…' if len(s) > 35 else s)
    sc_tab = sc_tab.reset_index(drop=True)
    sc_tab.insert(0, 'Rank', range(1, len(sc_tab)+1))
    t3 = sc_tab[['Rank','Scaffold SMILES','n_compounds','mean_pic50',
                  'pct_patent','hub_scaffold']].rename(columns={
        'n_compounds':'N compounds','mean_pic50':'Mean pIC50',
        'pct_patent':'% Patent','hub_scaffold':'Hub scaffold',
    })
    gt3 = (nature_gt(t3, 'Top 20 scaffold series (Murcko)',
                      source_note='✓ = dominant 174-compound azaindole-benzimidazole hub scaffold.')
           .cols_align(align='right', columns=['Rank','N compounds','Mean pIC50','% Patent'])
           .cols_align(align='left',  columns=['Scaffold SMILES','Hub scaffold']))
    save_gt(gt3, 'table3_scaffolds')
except Exception as e:
    print(f"  Table 3 RDKit error: {e}")

# ── TABLE 4: SALI Distribution Summary ───────────────────────────────────────
print("Table 4: SALI distribution")
tier_order = ['severe','moderate','broad','non_cliff']
rows_t4 = []
for tier in tier_order:
    sub = pairs[pairs['cliff_tier']==tier] if tier != 'non_cliff' \
          else pairs[~pairs['cliff_tier'].isin(['severe','moderate','broad'])]
    sali_sub = sub['sali'].dropna()
    n_cliff_pairs = len(cliffs[cliffs['cliff_tier']==tier]) if tier != 'non_cliff' else \
                    len(pairs) - len(cliffs)
    rows_t4.append({
        'Cliff tier': tier.replace('_',' ').capitalize(),
        'SALI pairs (n)': f'{len(sali_sub):,}',
        'Mean SALI': round(float(sali_sub.mean()), 2) if len(sali_sub) else '—',
        'Median SALI': round(float(sali_sub.median()), 2) if len(sali_sub) else '—',
        'Max SALI': round(float(sali_sub.max()), 2) if len(sali_sub) else '—',
        'Example cliff pairs': n_cliff_pairs,
    })
t4 = pd.DataFrame(rows_t4)
gt4 = (nature_gt(t4, 'SALI distribution by cliff tier',
                  source_note='SALI = |ΔpIC50| / (1 − Tanimoto). Pairs with Tanimoto=1.0 excluded (NaN).')
       .cols_align(align='right',
                   columns=['SALI pairs (n)','Mean SALI','Median SALI','Max SALI','Example cliff pairs'])
       .cols_align(align='left', columns=['Cliff tier']))
save_gt(gt4, 'table4_sali')

# ── TABLE 5: Hub Compound Summary ─────────────────────────────────────────────
print("Table 5: Hub compound summary")
hub_data = [
    {'ID':'A1', 'InChIKey (14 chars)':'SMADULGDNOCLOP',
     'Class':'A', 'pIC50':5.39, 'Severe pairs':15,
     '% of 94': f'{15/94*100:.1f}%', 'Scaffold type':'Series member',
     'Compounds in scaffold':174, 'Sources':'PC+BDB+ChEMBL'},
    {'ID':'A2', 'InChIKey (14 chars)':'RAVBZQAQTVGKIV',
     'Class':'A', 'pIC50':5.34, 'Severe pairs':12,
     '% of 94': f'{12/94*100:.1f}%', 'Scaffold type':'Series member',
     'Compounds in scaffold':174, 'Sources':'PC+BDB'},
    {'ID':'B1', 'InChIKey (14 chars)':'UDCDEKJNAMHBFH',
     'Class':'B', 'pIC50':4.30, 'Severe pairs':12,
     '% of 94': f'{12/94*100:.1f}%', 'Scaffold type':'Singleton',
     'Compounds in scaffold':1, 'Sources':'PC+BDB'},
    {'ID':'B2', 'InChIKey (14 chars)':'DVCKJOQIVOGXEI',
     'Class':'B', 'pIC50':4.30, 'Severe pairs':11,
     '% of 94': f'{11/94*100:.1f}%', 'Scaffold type':'Singleton',
     'Compounds in scaffold':1, 'Sources':'PC+BDB'},
]
t5 = pd.DataFrame(hub_data)
t5.insert(0, 'Hub class', ['Class A — series-embedded hub']*2 +
                          ['Class B — scaffold-singleton hub']*2)
gt5 = (nature_gt(t5, 'Cliff hub compound summary',
                  source_note=(
    'Severe cliff threshold: Tanimoto ≥ 0.8, ΔpIC50 ≥ 2.0 log units. '
    'Class A: members of the dominant 174-compound azaindole-benzimidazole series. '
    'Class B: unique Murcko scaffolds (scaffold singletons).'))
       .tab_style(
           style=style.borders(sides='left', weight='3px', color='#1A237E'),
           locations=loc.body(rows=[0,1]))
       .tab_style(
           style=style.borders(sides='left', weight='3px', color='#CC3311'),
           locations=loc.body(rows=[2,3]))
       .cols_align(align='right',
                   columns=['pIC50','Severe pairs','% of 94','Compounds in scaffold'])
       .cols_align(align='left',
                   columns=['Hub class','ID','InChIKey (14 chars)','Class',
                            'Scaffold type','Sources']))
save_gt(gt5, 'table5_hubs')

# ── TABLE 6: MMP Analysis Summary ─────────────────────────────────────────────
print("Table 6: MMP analysis")
t6 = pd.DataFrame([
    {'Section':'Overview', 'Metric':'Total MMP pairs',         'Value':'707',       'Note':'All 99 cliff compounds'},
    {'Section':'Overview', 'Metric':'MMP-validated severe cliffs','Value':'85 of 94 (90.4%)','Note':'Tanimoto ≥ 0.8, ΔpIC50 ≥ 2.0'},
    {'Section':'Overview', 'Metric':'Non-MMP severe cliffs',   'Value':'9',         'Note':'Scaffold hops / fragment merges'},
    {'Section':'MMP type breakdown', 'Metric':'Single R-group change', 'Value':'49','Note':'57.6% of MMP severe'},
    {'Section':'MMP type breakdown', 'Metric':'Small substituent',     'Value':'28','Note':'32.9% of MMP severe'},
    {'Section':'MMP type breakdown', 'Metric':'Medium substituent',    'Value':'8', 'Note':'9.4% of MMP severe'},
])
gt6 = (nature_gt(t6, 'MMP analysis summary',
                  source_note='MMP pairs identified by rdkit.Chem.MolStandardize. '
                  'Cliff tier assigned by Tanimoto similarity and ΔpIC50 thresholds.')
       .cols_align(align='right', columns=['Value'])
       .cols_align(align='left',  columns=['Section','Metric','Note']))
save_gt(gt6, 'table6_mmp')

# ── TABLE 7: Top 20 SALI Pairs (deduplicated) ─────────────────────────────────
print("Table 7: Top 20 SALI pairs")
pairs_sali = pairs[pairs['sali'].notna()].copy()
# Deduplicate: ensure each pair appears only once (normalize A/B order)
pairs_sali['pair_key'] = pairs_sali.apply(
    lambda r: '|'.join(sorted([r['inchi_key_a'], r['inchi_key_b']])), axis=1)
pairs_dedup = pairs_sali.drop_duplicates('pair_key').sort_values('sali', ascending=False).head(20)
t7 = pd.DataFrame({
    'Rank': range(1, len(pairs_dedup)+1),
    'Compound A': pairs_dedup['inchi_key_a'].str[:14].values,
    'Compound B': pairs_dedup['inchi_key_b'].str[:14].values,
    'Tanimoto': pairs_dedup['tanimoto'].round(3).values,
    'ΔpIC50':   pairs_dedup['delta_pic50'].round(3).values,
    'SALI':     pairs_dedup['sali'].round(2).values,
    'Cliff tier': pairs_dedup['cliff_tier'].fillna('non_cliff').str.capitalize().values,
})
gt7 = (nature_gt(t7, 'Top 20 SALI pairs (deduplicated)',
                  source_note='SALI = |ΔpIC50| / (1 − Tanimoto). Pairs are unique.')
       .tab_style(
           style=style.fill(color='#FFE5E0'),
           locations=loc.body(rows=[i for i,t in enumerate(t7['Cliff tier']) if t=='Severe']))
       .tab_style(
           style=style.fill(color='#FFF3E0'),
           locations=loc.body(rows=[i for i,t in enumerate(t7['Cliff tier']) if t=='Moderate']))
       .tab_style(
           style=style.fill(color='#E3F2FD'),
           locations=loc.body(rows=[i for i,t in enumerate(t7['Cliff tier']) if t=='Broad']))
       .cols_align(align='right', columns=['Rank','Tanimoto','ΔpIC50','SALI'])
       .cols_align(align='left',  columns=['Compound A','Compound B','Cliff tier']))
save_gt(gt7, 'table7_sali_pairs')

# ── TABLE 8: Patent Compound Summary ─────────────────────────────────────────
print("Table 8: Patent compound summary")
from scipy.stats import mannwhitneyu
pat_mask = df['source_list']=='pubchem_confirmatory'
pat_pic  = df.loc[pat_mask,  'pic50_consensus']
pub_pic  = df.loc[~pat_mask, 'pic50_consensus']
_, p_mwu = mannwhitneyu(pat_pic, pub_pic, alternative='two-sided')
p_str = 'p < 0.001' if p_mwu < 0.001 else f'p = {p_mwu:.3f}'

t8 = pd.DataFrame([
    {'Metric': 'Patent-exclusive compounds (PubChem-only)', 'Value': '233'},
    {'Metric': 'Published compounds (multi-source)',         'Value': '2,860'},
    {'Metric': 'Patent mean pIC50',                         'Value': '6.082'},
    {'Metric': 'Published mean pIC50',                      'Value': '6.588'},
    {'Metric': 'pIC50 difference (published − patent)',      'Value': '0.506'},
    {'Metric': 'Mann-Whitney U p-value',                    'Value': p_str},
    {'Metric': 'Patent-exclusive scaffolds',                'Value': '103'},
    {'Metric': 'Patent severe cliff contribution',           'Value': '1 of 94 pairs (1.1%)'},
])
gt8 = (nature_gt(t8, 'Patent compound analysis summary',
                  source_note='Patent-exclusive = compounds appearing only in PubChem bioassay AIDs. '
                  'pIC50 difference reflects lower potency of patent screening libraries vs '
                  'ChEMBL/BindingDB primary literature compounds.')
       .cols_align(align='right', columns=['Value'])
       .cols_align(align='left',  columns=['Metric']))
save_gt(gt8, 'table8_patent')

# ── TABLE 9: Source Independence Summary ──────────────────────────────────────
print("Table 9: Source independence")
combo_info = [
    ('bindingdb',                           1.0,  95,  'Source-exclusive'),
    ('pubchem_confirmatory',                1.0, 233,  'Source-exclusive'),
    ('chembl',                              1.0,  10,  'Source-exclusive'),
    ('chembl|pubchem_confirmatory',         0.7,  23,  'Likely independent'),
    ('bindingdb|chembl',                    0.6, 167,  'Likely independent'),
    ('bindingdb|pubchem_confirmatory',      0.5,1199,  'Partial redundancy'),
    ('bindingdb|chembl|pubchem_confirmatory',0.3,1366, 'High redundancy'),
]
source_labels = {
    'bindingdb':                             'BindingDB only',
    'pubchem_confirmatory':                  'PubChem only',
    'chembl':                                'ChEMBL only',
    'chembl|pubchem_confirmatory':           'ChEMBL + PubChem',
    'bindingdb|chembl':                      'BindingDB + ChEMBL',
    'bindingdb|pubchem_confirmatory':        'BindingDB + PubChem',
    'bindingdb|chembl|pubchem_confirmatory': 'BindingDB + ChEMBL + PubChem',
}
rows_t9 = []
for sl_key, score, n, interp in combo_info:
    sub = df[df['source_list']==sl_key]['pic50_consensus']
    rows_t9.append({
        'Source combination': source_labels[sl_key],
        'Independence score': score,
        'N compounds': f'{n:,}',
        'Mean pIC50': round(float(sub.mean()), 2) if len(sub) else '—',
        'is_true_multi-source (≥0.6)': '✓' if score >= 0.6 else '✗',
        'Interpretation': interp,
    })
t9 = pd.DataFrame(rows_t9).sort_values('Independence score', ascending=False)
gt9 = (nature_gt(t9, 'Source independence score by combination',
                  source_note=(
    'Single-source compounds (BDB only, PC only, ChEMBL only) receive independence score = 1.0 '
    'by definition, as no cross-source comparison is available. The is_true_multi_source flag '
    'for these compounds reflects source exclusivity, not experimental replication. '
    'Threshold = 0.6 (528 true multi-source, 2,565 pipeline redundancy).'))
       .cols_align(align='right', columns=['Independence score','N compounds','Mean pIC50'])
       .cols_align(align='left',
                   columns=['Source combination','is_true_multi-source (≥0.6)','Interpretation']))
save_gt(gt9, 'table9_independence')

print("\n=== All 9 tables saved ===")
for f in sorted(os.listdir(OUT)):
    if f.endswith('.html'):
        sz = os.path.getsize(f'{OUT}/{f}') // 1024
        print(f'  {f}  {sz} KB')
