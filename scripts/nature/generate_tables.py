"""
generate_tables.py — All main + supplementary tables for PAD4-DB v2
Outputs:
  outputs/tables/nature/csv/*.csv
  outputs/tables/nature/latex/*.tex
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

ROOT = Path('/home/nidhal/PAD4-db_V2')
CSV_OUT = ROOT / 'outputs/tables/nature_v2/csv'
TEX_OUT = ROOT / 'outputs/tables/nature_v2/latex'
CSV_OUT.mkdir(parents=True, exist_ok=True)
TEX_OUT.mkdir(parents=True, exist_ok=True)

HUB_IKS = {
    'A1': 'SMADULGDNOCLOP-GISFHXKWSA-N',
    'A2': 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
    'B1': 'UDCDEKJNAMHBFH-HSZRJFAPSA-N',
    'B2': 'DVCKJOQIVOGXEI-XMMPIXPASA-N',
}

CANON = {
    'n_compounds': 3093,
    'n_severe': 94, 'n_moderate': 193, 'n_broad': 580,
    'n_in_severe': 99, 'n_scaffolds': 1244, 'n_series': 375,
    'n_singletons': 869, 'largest_series': 174, 'n_patent': 233,
    'n_multi_06': 528, 'n_multi_07': 361,
}

print("=" * 60)
print("GENERATE TABLES")
print("=" * 60)

# ── Load data ──────────────────────────────────────────────────────────────────
df = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
cliffs = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
pairs = pd.read_parquet(ROOT / 'data/processed/activity_pairs_with_sali.parquet')
mmp = pd.read_csv(ROOT / 'outputs/mmp/mmp_pairs_cliff99.csv')
disc = pd.read_csv(ROOT / 'outputs/mmp/mmp_discontinuity_scores.csv')
scaffold_sum = pd.read_csv(ROOT / 'outputs/tables/05_scaffold_summary.csv')
scaffold_sum = scaffold_sum.sort_values('n_compounds', ascending=False).reset_index(drop=True)
scaffold_sum['scaffold_rank'] = scaffold_sum.index + 1

print(f"  Loaded: {len(df)} cpds, {len(cliffs)} cliffs, {len(pairs)} pairs, {len(mmp)} mmp pairs")

assert len(df) == CANON['n_compounds']

# ── Helper: LaTeX table ───────────────────────────────────────────────────────
def to_latex(df_tbl, caption, label, col_fmt=None, float_fmt='%.2f'):
    n_cols = len(df_tbl.columns)
    if col_fmt is None:
        col_fmt = '@{}l' + 'r' * (n_cols - 1) + '@{}'
    header = ' & '.join([str(c).replace('_', ' ') for c in df_tbl.columns]) + ' \\\\'

    rows = []
    for _, row in df_tbl.iterrows():
        cells = []
        for val in row:
            if isinstance(val, float):
                if np.isnan(val):
                    cells.append('—')
                elif abs(val) < 0.001:
                    cells.append(f'{val:.2e}')
                else:
                    cells.append(f'{val:.2f}')
            else:
                cells.append(str(val).replace('_', '\\_').replace('&', '\\&').replace('%', '\\%'))
        rows.append(' & '.join(cells) + ' \\\\')

    tex = f"""\\begin{{table}}[h]
\\centering
\\caption{{{caption}}}
\\label{{{label}}}
\\begin{{tabular}}{{{col_fmt}}}
\\toprule
{header}
\\midrule
""" + '\n'.join(rows) + """
\\bottomrule
\\end{tabular}
\\end{table}
"""
    return tex

# ── Table 1 — Cliff Hub Summary ───────────────────────────────────────────────
print("\n[Table 1] Cliff hub summary ...")
disc_joined = disc.merge(df[['inchi_key', 'pIC50', 'murcko_smiles', 'scaffold_rank',
                               'scaffold_series_size']], on='inchi_key', how='left')

hub_rows = []
for hub_id, ik in HUB_IKS.items():
    row = disc_joined[disc_joined['inchi_key'] == ik]
    if len(row) == 0:
        print(f"  WARNING: {hub_id} ({ik}) not found in disc_joined")
        hub_rows.append({
            'ID': hub_id,
            'InChIKey (14 chars)': ik[:14],
            'Class': 'A' if hub_id.startswith('A') else 'B',
            'pIC50': float('nan'),
            'Severe pairs': float('nan'),
            '% of 94': float('nan'),
            'In scaffold': float('nan'),
            'Series size': float('nan'),
        })
        continue
    r = row.iloc[0]
    n_pairs = int(r['severe_cliff_degree']) if not pd.isna(r['severe_cliff_degree']) else 0
    pct = n_pairs / CANON['n_severe'] * 100
    hub_rows.append({
        'ID': hub_id,
        'InChIKey (14 chars)': ik[:14],
        'Class': 'A' if hub_id.startswith('A') else 'B',
        'pIC50': round(r['pIC50'], 2) if not pd.isna(r['pIC50']) else float('nan'),
        'Severe pairs': n_pairs,
        '% of 94': round(pct, 1),
        'In scaffold': int(r['scaffold_rank']) if not pd.isna(r['scaffold_rank']) else '—',
        'Series size': int(r['scaffold_series_size']) if not pd.isna(r['scaffold_series_size']) else '—',
    })

t1 = pd.DataFrame(hub_rows)
print(t1)
t1.to_csv(CSV_OUT / 'table1_hub_summary.csv', index=False)
tex1 = to_latex(t1, 'Activity cliff hub compound summary\\footnotemark', 'tab:hubs',
                 col_fmt='@{}llrrrrr@{}')
tex1 = tex1.replace('\\end{table}',
    '\\footnotetext{Scaffold rank ordered by series size (rank 1 = 174-compound series). '
    'Hub A compounds belong to the dominant azaindole-benzimidazole scaffold series. '
    'Hub B compounds are in smaller series (ranks 526/509).}\n\\end{table}')
(TEX_OUT / 'table1_hub_summary.tex').write_text(tex1)
print("  Saved Table 1")

# ── Table 2 — Activity Cliff Summary ─────────────────────────────────────────
print("\n[Table 2] Activity cliff summary ...")
cliff_rows = []
total_iks = set()
for tier in ['severe', 'moderate', 'broad']:
    sub = cliffs[cliffs['cliff_tier'] == tier]
    iks = set(sub['inchi_key_a'].tolist() + sub['inchi_key_b'].tolist())
    total_iks.update(iks)
    delta = sub['delta_pic50'].abs()
    cliff_rows.append({
        'Tier': tier.capitalize(),
        'N pairs': len(sub),
        'N compounds': len(iks),
        '% of dataset': round(len(iks) / CANON['n_compounds'] * 100, 1),
        'Median |ΔpIC50|': round(delta.median(), 3),
        'Max |ΔpIC50|': round(delta.max(), 3),
    })

# Total row
delta_all = cliffs['delta_pic50'].abs()
cliff_rows.append({
    'Tier': 'Total',
    'N pairs': len(cliffs),
    'N compounds': len(total_iks),
    '% of dataset': round(len(total_iks) / CANON['n_compounds'] * 100, 1),
    'Median |ΔpIC50|': round(delta_all.median(), 3),
    'Max |ΔpIC50|': round(delta_all.max(), 3),
})

t2 = pd.DataFrame(cliff_rows)
print(t2)
t2.to_csv(CSV_OUT / 'table2_cliff_summary.csv', index=False)
tex2 = to_latex(t2, 'Activity cliff tier summary', 'tab:cliffs')
(TEX_OUT / 'table2_cliff_summary.tex').write_text(tex2)
print("  Saved Table 2")

# ── Table 3 — MMP Analysis Summary ───────────────────────────────────────────
print("\n[Table 3] MMP analysis summary ...")
mmp_types = ['single_atom_change', 'small_substituent', 'medium_substituent', 'large_substituent']
type_labels = ['Single atom change', 'Small substituent', 'Medium substituent', 'Large substituent']

n_severe_mmp = (mmp['is_canonical_severe_cliff'] == True).sum()  # 85 MMP-validated severe cliffs
mmp_rows = []
for mtype, mlabel in zip(mmp_types, type_labels):
    sub = mmp[mmp['mmp_type'] == mtype]
    in_severe = sub[sub['cliff_tier'] == 'severe']
    pct_val = len(in_severe) / n_severe_mmp * 100 if n_severe_mmp > 0 else 0
    mmp_rows.append({
        'MMP type': mlabel,
        'N pairs': len(sub),
        'In severe cliffs': len(in_severe),
        '% of 85 validated': round(pct_val, 1),
    })

t3 = pd.DataFrame(mmp_rows)
print(t3)
t3.to_csv(CSV_OUT / 'table3_mmp_summary.csv', index=False)
tex3 = to_latex(t3, 'MMP analysis summary by transformation type', 'tab:mmp')
(TEX_OUT / 'table3_mmp_summary.tex').write_text(tex3)
print("  Saved Table 3")

# ── Table 4 — Source Independence by Combination ─────────────────────────────
print("\n[Table 4] Source independence by combination ...")
df['has_pubchem'] = df['source_list'].str.contains('pubchem', na=False)
df['has_chembl'] = df['source_list'].str.contains('chembl', na=False)
df['has_binding'] = df['source_list'].str.contains('bindingdb', na=False)

def get_combo(row):
    parts = []
    if row['has_pubchem']: parts.append('PubChem')
    if row['has_chembl']:  parts.append('ChEMBL')
    if row['has_binding']: parts.append('BindingDB')
    return '+'.join(sorted(parts))

df['combo'] = df.apply(get_combo, axis=1)

t4_rows = []
for combo, sub in df.groupby('combo'):
    sis = sub['source_independence_score']
    t4_rows.append({
        'Source combination': combo,
        'N compounds': len(sub),
        'Mean score': round(sis.mean(), 3),
        'Median score': round(sis.median(), 3),
        'Threshold met (>=0.6)': (sis >= 0.6).sum(),
    })

t4 = pd.DataFrame(t4_rows).sort_values('N compounds', ascending=False)
print(t4)
t4.to_csv(CSV_OUT / 'table4_source_independence.csv', index=False)
tex4 = to_latex(t4, 'Source independence scores by source combination', 'tab:source_indep')
(TEX_OUT / 'table4_source_independence.tex').write_text(tex4)
print("  Saved Table 4")

# ── Table S1 — Source Coverage ────────────────────────────────────────────────
print("\n[Table S1] Source coverage ...")
src_rows = []
for src_name, src_col in [('PubChem', 'has_pubchem'), ('ChEMBL', 'has_chembl'),
                            ('BindingDB', 'has_binding')]:
    sub = df[df[src_col]]
    p = sub['pIC50']
    src_rows.append({
        'Source': src_name,
        'N': len(sub),
        '%': round(len(sub) / CANON['n_compounds'] * 100, 1),
        'Mean pIC50': round(p.mean(), 2),
        'Median pIC50': round(p.median(), 2),
        'SD': round(p.std(), 2),
    })

# All three
all_three = df[df['has_pubchem'] & df['has_chembl'] & df['has_binding']]
p_all = all_three['pIC50']
src_rows.append({
    'Source': 'All three',
    'N': len(all_three),
    '%': round(len(all_three) / CANON['n_compounds'] * 100, 1),
    'Mean pIC50': round(p_all.mean(), 2),
    'Median pIC50': round(p_all.median(), 2),
    'SD': round(p_all.std(), 2),
})

ts1 = pd.DataFrame(src_rows)
print(ts1)
ts1.to_csv(CSV_OUT / 'tableS1_source_coverage.csv', index=False)
tex_s1 = to_latex(ts1, 'Source database coverage and pIC50 statistics', 'tab:source_cov')
(TEX_OUT / 'tableS1_source_coverage.tex').write_text(tex_s1)
print("  Saved Table S1")

# ── Table S2 — pIC50 by Assay Mechanism ──────────────────────────────────────
print("\n[Table S2] pIC50 by mechanism ...")
from scipy.stats import kruskal
mechs = ['enzymatic', 'enzymatic_confirmed', 'fp_ic50', 'covalent']
mech_data_list = [df[df['mechanism_class'] == m]['pIC50'].dropna().values for m in mechs]
kw_stat, kw_p = kruskal(*mech_data_list)
print(f"  Kruskal-Wallis H={kw_stat:.2f}, p={kw_p:.2e}")

ts2_rows = []
for mech, data in zip(mechs, mech_data_list):
    ts2_rows.append({
        'Mechanism': mech.replace('_', ' ').title(),
        'N': len(data),
        'Mean pIC50': round(np.mean(data), 2),
        'Median pIC50': round(np.median(data), 2),
        'SD': round(np.std(data), 2),
        'Min': round(np.min(data), 2),
        'Max': round(np.max(data), 2),
    })

ts2 = pd.DataFrame(ts2_rows)
print(ts2)
ts2.to_csv(CSV_OUT / 'tableS2_mechanism_pic50.csv', index=False)
tex_s2 = to_latex(ts2, f'pIC50 by assay mechanism (Kruskal-Wallis H={kw_stat:.2f}, p={kw_p:.2e})',
                   'tab:mech_pic50')
(TEX_OUT / 'tableS2_mechanism_pic50.tex').write_text(tex_s2)
print("  Saved Table S2")

# ── Table S3 — SALI Distribution ─────────────────────────────────────────────
print("\n[Table S3] SALI distribution ...")
ts3_rows = []
for tier in ['severe', 'moderate', 'broad']:
    sub = pairs[pairs['cliff_tier'] == tier]
    ts3_rows.append({
        'Tier': tier.capitalize(),
        'N pairs': len(sub),
        'SALI>10': (sub['sali'] > 10).sum(),
        'SALI>20': (sub['sali'] > 20).sum(),
        'Max SALI': round(sub['sali'].max(), 2),
        'Mean |ΔpIC50|': round(sub['delta_pic50'].abs().mean(), 3),
    })

# All pairs
ts3_rows.append({
    'Tier': 'All',
    'N pairs': len(pairs),
    'SALI>10': (pairs['sali'] > 10).sum(),
    'SALI>20': (pairs['sali'] > 20).sum(),
    'Max SALI': round(pairs['sali'].max(), 2),
    'Mean |ΔpIC50|': round(pairs['delta_pic50'].abs().mean(), 3),
})

ts3 = pd.DataFrame(ts3_rows)
print(ts3)
ts3.to_csv(CSV_OUT / 'tableS3_sali_distribution.csv', index=False)
tex_s3 = to_latex(ts3, 'SALI statistics by activity cliff tier', 'tab:sali')
(TEX_OUT / 'tableS3_sali_distribution.tex').write_text(tex_s3)
print("  Saved Table S3")

# ── Table S4 — Patent Analysis ────────────────────────────────────────────────
print("\n[Table S4] Patent analysis ...")
n_patent = df['patent_flag'].sum()
scaffold_sum_s4 = pd.read_csv(ROOT / 'outputs/tables/05_scaffold_summary.csv')
n_patent_scaffolds = scaffold_sum_s4['patent_exclusive_scaffold'].sum()
cliffs_patent = cliffs[cliffs['any_patent_exclusive'] == True]
n_patent_cliffs = (cliffs_patent['cliff_tier'] == 'severe').sum()

ts4_rows = [
    {'Metric': 'Patent-exclusive compounds', 'Value': n_patent},
    {'Metric': 'Patent-exclusive scaffolds', 'Value': n_patent_scaffolds},
    {'Metric': 'Patent compounds in severe cliffs', 'Value': n_patent_cliffs},
    {'Metric': 'Patent pIC50 mean', 'Value': round(df[df['patent_flag']]['pIC50'].mean(), 2)},
    {'Metric': 'Patent pIC50 median', 'Value': round(df[df['patent_flag']]['pIC50'].median(), 2)},
    {'Metric': 'Published pIC50 mean', 'Value': round(df[~df['patent_flag']]['pIC50'].mean(), 2)},
    {'Metric': 'Published pIC50 median', 'Value': round(df[~df['patent_flag']]['pIC50'].median(), 2)},
]

ts4 = pd.DataFrame(ts4_rows)
print(ts4)
ts4.to_csv(CSV_OUT / 'tableS4_patent_analysis.csv', index=False)
tex_s4 = to_latex(ts4, 'Patent-exclusive compound analysis', 'tab:patent',
                   col_fmt='@{}lr@{}')
(TEX_OUT / 'tableS4_patent_analysis.tex').write_text(tex_s4)
print("  Saved Table S4")

# ── Table S5 — Top 20 SALI Pairs ─────────────────────────────────────────────
print("\n[Table S5] Top 20 SALI pairs ...")
pairs_dedup = pairs.copy()
pairs_dedup['_pair_key'] = pairs_dedup.apply(
    lambda r: tuple(sorted([r['inchi_key_a'], r['inchi_key_b']])), axis=1)
pairs_dedup = pairs_dedup.drop_duplicates('_pair_key')
top20_sali = pairs_dedup.nlargest(20, 'sali').reset_index(drop=True)
top20_sali['rank'] = top20_sali.index + 1
top20_sali['ik_a_14'] = top20_sali['inchi_key_a'].str[:14]
top20_sali['ik_b_14'] = top20_sali['inchi_key_b'].str[:14]

ts5 = top20_sali[['rank', 'ik_a_14', 'ik_b_14', 'sali',
                    'delta_pic50', 'tanimoto', 'cliff_tier']].copy()
ts5.columns = ['Rank', 'InChIKey A (14)', 'InChIKey B (14)', 'SALI',
                '|ΔpIC50|', 'Tanimoto', 'Tier']
ts5['SALI'] = ts5['SALI'].round(2)
ts5['|ΔpIC50|'] = ts5['|ΔpIC50|'].abs().round(3)
ts5['Tanimoto'] = ts5['Tanimoto'].round(3)
ts5['Tier'] = ts5['Tier'].str.capitalize()

print(ts5.head(5))
ts5.to_csv(CSV_OUT / 'tableS5_top20_sali_pairs.csv', index=False)
tex_s5 = to_latex(ts5, 'Top 20 SALI activity cliff pairs', 'tab:top_sali',
                   col_fmt='@{}rllrrrr@{}')
(TEX_OUT / 'tableS5_top20_sali_pairs.tex').write_text(tex_s5)
print("  Saved Table S5")

# ── Table S6 — Full Compound List (CSV only) ──────────────────────────────────
print("\n[Table S6] Full compound list ...")
ts6_cols = ['inchi_key', 'smiles_std', 'pIC50', 'source_list', 'mechanism_class',
             'patent_flag', 'source_independence_score']
ts6 = df[ts6_cols].copy()
ts6['pIC50'] = ts6['pIC50'].round(4)
ts6['source_independence_score'] = ts6['source_independence_score'].round(4)
ts6.to_csv(CSV_OUT / 'tableS6_full_compound_list.csv', index=False)
print(f"  Saved Table S6 ({len(ts6)} rows) — CSV only")

# ── Summary ───────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("TABLE GENERATION COMPLETE")
print("=" * 60)
csv_files = list(CSV_OUT.glob('*.csv'))
tex_files = list(TEX_OUT.glob('*.tex'))
print(f"  CSV files: {len(csv_files)}")
print(f"  LaTeX files: {len(tex_files)}")
for f in sorted(csv_files):
    print(f"    {f.name}")
for f in sorted(tex_files):
    print(f"    {f.name}")
