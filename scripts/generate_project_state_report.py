# PAD4-DB v2 — Full Project State Report
# Run this from /home/nidhal/PAD4-db_V2 with pad4bench active
# Purpose: generate a complete snapshot of everything done so far
# so the human can paste the output into a new chat session

import pandas as pd
import numpy as np
import os, json, datetime
from pathlib import Path

ROOT = Path('/home/nidhal/PAD4-db_V2')
REPORT = []

def section(title):
    line = '=' * 70
    REPORT.append(f'\n{line}\n  {title}\n{line}')
    print(f'\n{line}\n  {title}\n{line}')

def log(msg):
    REPORT.append(str(msg))
    print(msg)

def warn(msg):
    msg = f'  ⚠️  {msg}'
    REPORT.append(msg)
    print(msg)

def ok(msg):
    msg = f'  ✅  {msg}'
    REPORT.append(msg)
    print(msg)

def fail(msg):
    msg = f'  ❌  {msg}'
    REPORT.append(msg)
    print(msg)

# ─────────────────────────────────────────────────────────────────────
section('0. ENVIRONMENT')
# ─────────────────────────────────────────────────────────────────────
import sys
log(f'Python:      {sys.version.split()[0]}')
log(f'Root:        {ROOT}')
log(f'Run date:    {datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}')

try:
    from rdkit import Chem, rdBase
    log(f'RDKit:       {rdBase.rdkitVersion}')
except:
    warn('RDKit not importable')

try:
    import sklearn
    log(f'scikit-learn: {sklearn.__version__}')
except:
    warn('scikit-learn not importable')

# ─────────────────────────────────────────────────────────────────────
section('1. MASTER PARQUET — pad4_compounds.parquet')
# ─────────────────────────────────────────────────────────────────────
master_path = ROOT / 'data/processed/pad4_compounds.parquet'
if master_path.exists():
    df = pd.read_parquet(master_path)
    log(f'Rows:        {len(df):,}')
    log(f'Columns:     {len(df.columns)}')
    log(f'Column list: {df.columns.tolist()}')

    # pIC50 stats — handle both column name variants
    pic50_col = None
    for candidate in ['pic50_consensus', 'pIC50', 'pic50']:
        if candidate in df.columns:
            pic50_col = candidate
            break
    if pic50_col:
        log(f'\npIC50 column: {pic50_col}')
        log(f'  min:    {df[pic50_col].min():.4f}')
        log(f'  max:    {df[pic50_col].max():.4f}')
        log(f'  mean:   {df[pic50_col].mean():.4f}')
        log(f'  median: {df[pic50_col].median():.4f}')
        log(f'  std:    {df[pic50_col].std():.4f}')
    else:
        warn('No pIC50 column found')

    # New columns
    for col in ['is_covalent', 'warhead_class', 'mechanism_class',
                'fragment_flag']:
        if col in df.columns:
            ok(f'{col} present')
            if df[col].dtype == bool:
                log(f'    True: {df[col].sum()}, False: {(~df[col]).sum()}')
            else:
                log(f'    value_counts:\n{df[col].value_counts(dropna=False).to_string()}')
        else:
            warn(f'{col} MISSING from master parquet')

    # Source combination breakdown
    if 'source_list' in df.columns:
        log('\nSource combination breakdown:')
        for combo, count in df['source_list'].value_counts().items():
            log(f'  {combo:<55} {count:>5}')
else:
    fail(f'Master parquet not found: {master_path}')
    df = None

# ─────────────────────────────────────────────────────────────────────
section('2. ALL PROCESSED PARQUET FILES')
# ─────────────────────────────────────────────────────────────────────
processed_dir = ROOT / 'data/processed'
if processed_dir.exists():
    for f in sorted(processed_dir.glob('*.parquet')):
        try:
            tmp = pd.read_parquet(f)
            size_mb = f.stat().st_size / 1e6
            log(f'  {f.name:<50} {len(tmp):>8,} rows  '
                f'{len(tmp.columns):>3} cols  {size_mb:.1f} MB')
            log(f'    columns: {tmp.columns.tolist()}')
        except Exception as e:
            warn(f'  {f.name}: could not read — {e}')
else:
    warn('data/processed/ directory not found')

# ─────────────────────────────────────────────────────────────────────
section('3. INTERIM FILES')
# ─────────────────────────────────────────────────────────────────────
interim_dir = ROOT / 'data/interim'
if interim_dir.exists():
    for f in sorted(interim_dir.rglob('*.parquet')) + \
             sorted(interim_dir.rglob('*.npy')):
        size_mb = f.stat().st_size / 1e6
        log(f'  {str(f.relative_to(ROOT)):<60} {size_mb:.1f} MB')
else:
    warn('data/interim/ not found')

# ─────────────────────────────────────────────────────────────────────
section('4. ACTIVITY CLIFFS')
# ─────────────────────────────────────────────────────────────────────
cliff_path = ROOT / 'data/processed/activity_cliffs.parquet'
if cliff_path.exists():
    cliffs = pd.read_parquet(cliff_path)
    log(f'Total cliff pairs: {len(cliffs):,}')
    log(f'Columns: {cliffs.columns.tolist()}')

    # Detect delta column
    delta_col = None
    for candidate in ['delta_pic50', 'delta_pIC50', 'deltapic50']:
        if candidate in cliffs.columns:
            delta_col = candidate
            break

    # Detect tier column
    tier_col = None
    for candidate in ['cliff_tier', 'tier']:
        if candidate in cliffs.columns:
            tier_col = candidate
            break

    if delta_col:
        severe   = cliffs[cliffs[delta_col] >= 2.0]
        moderate = cliffs[(cliffs[delta_col] >= 1.5) &
                          (cliffs[delta_col] < 2.0)]
        broad    = cliffs[(cliffs[delta_col] >= 1.0) &
                          (cliffs[delta_col] < 1.5)]
        log(f'Severe   (Δ≥2.0): {len(severe)}')
        log(f'Moderate (Δ≥1.5): {len(moderate)}')
        log(f'Broad    (Δ≥1.0): {len(broad)}')
        log(f'Max delta_pic50:   {cliffs[delta_col].max():.4f}')

        # Hub degrees
        ik_cols = [c for c in cliffs.columns
                   if 'inchi_key' in c.lower() or c in ('ik_a','ik_b',
                   'compound_a_ik','compound_b_ik')]
        if len(ik_cols) >= 2:
            col_a, col_b = ik_cols[0], ik_cols[1]
            all_iks = pd.concat([severe[col_a], severe[col_b]])
            degree = all_iks.value_counts()
            log('\nTop 10 compounds by severe cliff degree:')
            for ik, cnt in degree.head(10).items():
                log(f'  {ik}  degree={cnt}')
    if tier_col:
        log(f'\nTier column ({tier_col}) value_counts:')
        log(cliffs[tier_col].value_counts().to_string())
else:
    warn('activity_cliffs.parquet not found')

# ─────────────────────────────────────────────────────────────────────
section('5. SALI PAIRS FILE')
# ─────────────────────────────────────────────────────────────────────
sali_path = ROOT / 'data/processed/activity_pairs_with_sali.parquet'
if sali_path.exists():
    pairs = pd.read_parquet(sali_path)
    log(f'Total pairs:       {len(pairs):,}')
    log(f'Columns:           {pairs.columns.tolist()}')
    if 'sali' in pairs.columns or 'SALI' in pairs.columns:
        sc = 'sali' if 'sali' in pairs.columns else 'SALI'
        log(f'SALI NaN count:    {pairs[sc].isna().sum()}')
        log(f'SALI max:          {pairs[sc].max():.4f}')
        log(f'SALI > 10:         {(pairs[sc] > 10).sum()}')
        log(f'SALI > 20:         {(pairs[sc] > 20).sum()}')
    if 'tanimoto' in pairs.columns:
        log(f'Pairs sim>=0.8:    {(pairs["tanimoto"] >= 0.8).sum():,}')
else:
    warn('activity_pairs_with_sali.parquet not found')
    alt = ROOT / 'data/processed/activity_pairs_sim_ge06.parquet'
    if alt.exists():
        p2 = pd.read_parquet(alt)
        log(f'activity_pairs_sim_ge06.parquet: {len(p2):,} rows')

# ─────────────────────────────────────────────────────────────────────
section('6. MMP OUTPUTS')
# ─────────────────────────────────────────────────────────────────────
mmp_dir = ROOT / 'outputs/mmp'
if mmp_dir.exists():
    for f in sorted(mmp_dir.glob('*.csv')):
        tmp = pd.read_csv(f)
        log(f'  {f.name:<50} {len(tmp):>6,} rows')
        log(f'    columns: {tmp.columns.tolist()}')
        if 'cliff_tier' in tmp.columns:
            log(f'    tier breakdown: '
                f'{tmp["cliff_tier"].value_counts().to_dict()}')
        if 'mmp_type' in tmp.columns:
            log(f'    mmp_type breakdown: '
                f'{tmp["mmp_type"].value_counts().to_dict()}')
        if 'discontinuity_score' in tmp.columns:
            top = tmp.nlargest(1, 'discontinuity_score')
            log(f'    top discontinuity: '
                f'{top.iloc[0]["inchi_key"][:20]} '
                f'score={top.iloc[0]["discontinuity_score"]:.4f}')
else:
    warn('outputs/mmp/ not found')

# ─────────────────────────────────────────────────────────────────────
section('7. SCAFFOLD FILES')
# ─────────────────────────────────────────────────────────────────────
scaffold_csv = ROOT / 'outputs/tables/05_scaffold_summary.csv'
if scaffold_csv.exists():
    sc_df = pd.read_csv(scaffold_csv)
    log(f'Scaffold summary rows: {len(sc_df):,}')
    log(f'Columns: {sc_df.columns.tolist()}')
    # Find size column
    size_col = None
    for c in ['n_compounds','series_size','count','size']:
        if c in sc_df.columns:
            size_col = c
            break
    if size_col:
        log(f'Unique scaffolds:      {len(sc_df)}')
        log(f'Series (>=2):          {(sc_df[size_col] >= 2).sum()}')
        log(f'Singletons:            {(sc_df[size_col] == 1).sum()}')
        log(f'Largest series:        {sc_df[size_col].max()}')
else:
    warn(f'05_scaffold_summary.csv not found at {scaffold_csv}')

# ─────────────────────────────────────────────────────────────────────
section('8. AUDIT OUTPUTS')
# ─────────────────────────────────────────────────────────────────────
audit_dir = ROOT / 'outputs/audit'
if audit_dir.exists():
    for f in sorted(audit_dir.glob('*')):
        size_kb = f.stat().st_size / 1e3
        log(f'  {f.name:<50} {size_kb:.1f} KB')
else:
    warn('outputs/audit/ not found')

# Final audit result
master_verify = ROOT / 'outputs/audit/MASTER_VERIFICATION_TABLE.txt'
if master_verify.exists():
    content = master_verify.read_text()
    pass_count = content.count('PASS')
    fail_count = content.count('FAIL')
    log(f'\nMASTER_VERIFICATION_TABLE: {pass_count} PASS, {fail_count} FAIL')

# ─────────────────────────────────────────────────────────────────────
section('9. FIGURE FILES — ALL LOCATIONS')
# ─────────────────────────────────────────────────────────────────────
fig_dirs = [
    ROOT / 'outputs/figures',
    ROOT / 'outputs/figures/nature',
    ROOT / 'outputs/figures/nature_v2',
]
for fig_dir in fig_dirs:
    if fig_dir.exists():
        pngs = sorted(fig_dir.glob('*.png'))
        svgs = sorted(fig_dir.glob('*.svg'))
        pdfs = sorted(fig_dir.glob('*.pdf'))
        log(f'\n{fig_dir.relative_to(ROOT)}: '
            f'{len(pngs)} PNG, {len(svgs)} SVG, {len(pdfs)} PDF')
        for f in sorted(fig_dir.glob('*')):
            if f.suffix in ('.png', '.svg', '.pdf'):
                size_kb = f.stat().st_size / 1e3
                log(f'  {f.name:<55} {size_kb:>7.1f} KB')

# ─────────────────────────────────────────────────────────────────────
section('10. TABLE FILES — ALL LOCATIONS')
# ─────────────────────────────────────────────────────────────────────
table_dirs = [
    ROOT / 'outputs/tables',
    ROOT / 'outputs/tables/nature',
    ROOT / 'outputs/tables/nature_v2',
    ROOT / 'outputs/tables/latex',
    ROOT / 'outputs/tables/csv',
]
for tdir in table_dirs:
    if tdir.exists():
        files = sorted(tdir.glob('*'))
        files = [f for f in files if f.is_file()]
        log(f'\n{tdir.relative_to(ROOT)}: {len(files)} files')
        for f in files:
            size_kb = f.stat().st_size / 1e3
            log(f'  {f.name:<55} {size_kb:>7.1f} KB')

# ─────────────────────────────────────────────────────────────────────
section('11. SCRIPTS INVENTORY')
# ─────────────────────────────────────────────────────────────────────
script_dirs = [
    ROOT / 'scripts',
    ROOT / 'outputs/figures',
]
for sdir in script_dirs:
    if sdir.exists():
        pyscripts = sorted(sdir.rglob('*.py'))
        if pyscripts:
            log(f'\n{sdir.relative_to(ROOT)}: {len(pyscripts)} Python scripts')
            for f in pyscripts:
                mtime = datetime.datetime.fromtimestamp(
                    f.stat().st_mtime).strftime('%Y-%m-%d %H:%M')
                log(f'  {str(f.relative_to(ROOT)):<65} {mtime}')

# ─────────────────────────────────────────────────────────────────────
section('12. CACHED ASSETS')
# ─────────────────────────────────────────────────────────────────────
cache_files = [
    ROOT / 'data/interim/morgan_fps_3093.npy',
    ROOT / 'data/interim/tsne_coords_3093.npy',
    ROOT / 'data/interim/shared_assets.parquet',
]
for f in cache_files:
    if f.exists():
        size_mb = f.stat().st_size / 1e6
        ok(f'{f.name}: {size_mb:.1f} MB')
    else:
        warn(f'{f.name}: NOT FOUND')

# Also scan interim for any other cached files
interim = ROOT / 'data/interim'
if interim.exists():
    for f in sorted(interim.rglob('*')):
        if f.is_file() and f not in cache_files:
            size_mb = f.stat().st_size / 1e6
            log(f'  {str(f.relative_to(ROOT)):<60} {size_mb:.1f} MB')

# ─────────────────────────────────────────────────────────────────────
section('13. LOCKED NUMBERS SPOT CHECK')
# ─────────────────────────────────────────────────────────────────────
LOCKED = {
    'SAR compounds': 3093,
    'HTS compounds': 327336,
    'SAR∩HTS overlap': 1453,
    'Total unique InChIKeys': 328976,
    'Severe cliff pairs': 94,
    'Moderate cliff pairs': 193,
    'Broad cliff pairs': 580,
    'Severe cliff compounds': 99,
    'Total SALI pairs': 358416,
    'SALI NaN (Tanimoto=1.0)': 198,
    'SALI max': 65.88,
    'Unique scaffolds': 1244,
    'Series scaffolds': 375,
    'Singleton scaffolds': 869,
    'Largest series': 174,
    'Patent-exclusive compounds': 233,
    'is_true_multi_source True': 528,
    'is_covalent True': 107,
    'fragment_flag True': 5,
    'MMP pairs (99 compounds)': 707,
    'MMP-confirmed severe cliffs': 85,
    'Hub A1 degree': 15,
    'Hub A2 degree': 12,
    'Hub B1 degree': 12,
    'Hub B2 degree': 11,
    'Collective hub pairs': 50,
}

log('\nSpot-checking locked numbers against live data:')
log(f'{"Metric":<40} {"Locked":>10} {"Status"}')
log('-' * 65)

if df is not None:
    n = len(df)
    status = '✅' if n == 3093 else '❌'
    log(f'{"SAR compounds":<40} {"3,093":>10} {status} actual={n}')

    if pic50_col:
        mx = df[pic50_col].max()
        status = '✅' if abs(mx - 8.52) < 0.05 else '❌'
        log(f'{"pIC50 max":<40} {"8.52":>10} {status} actual={mx:.4f}')

if cliff_path.exists() and delta_col:
    sv = len(cliffs[cliffs[delta_col] >= 2.0])
    status = '✅' if sv == 94 else '❌'
    log(f'{"Severe cliff pairs":<40} {"94":>10} {status} actual={sv}')

if sali_path.exists() and 'sali' in pairs.columns:
    tot = len(pairs)
    status = '✅' if tot == 358416 else '❌'
    log(f'{"Total SALI pairs":<40} {"358,416":>10} {status} actual={tot}')
    mx = pairs['sali'].max()
    status = '✅' if abs(mx - 65.88) < 0.1 else '❌'
    log(f'{"SALI max":<40} {"65.88":>10} {status} actual={mx:.4f}')

# ─────────────────────────────────────────────────────────────────────
section('14. CLAUDE.md SNAPSHOT')
# ─────────────────────────────────────────────────────────────────────
claude_md = ROOT / 'CLAUDE.md'
if claude_md.exists():
    content = claude_md.read_text()
    size_kb = claude_md.stat().st_size / 1e3
    log(f'CLAUDE.md exists: {size_kb:.1f} KB')
    log(f'First 3,000 characters:')
    log(content[:3000])
    log('...[truncated]...')
else:
    warn('CLAUDE.md not found')

# ─────────────────────────────────────────────────────────────────────
section('15. DIRECTORY TREE SUMMARY')
# ─────────────────────────────────────────────────────────────────────
import subprocess
result = subprocess.run(
    ['find', str(ROOT), '-maxdepth', '3', '-type', 'f',
     '-not', '-path', '*/raw/*',
     '-not', '-path', '*/.git/*',
     '-not', '-path', '*/node_modules/*'],
    capture_output=True, text=True
)
file_list = sorted(result.stdout.strip().split('\n'))
log(f'\nTotal files (excl. raw data): {len(file_list)}')
for f in file_list:
    if f:
        try:
            size = Path(f).stat().st_size
            log(f'  {f.replace(str(ROOT)+"/",""):<70} {size/1e3:>8.1f} KB')
        except:
            pass

# ─────────────────────────────────────────────────────────────────────
section('16. SAVE REPORT')
# ─────────────────────────────────────────────────────────────────────
report_path = ROOT / 'outputs/PROJECT_STATE_REPORT.txt'
report_path.parent.mkdir(parents=True, exist_ok=True)
full_report = '\n'.join(REPORT)
report_path.write_text(full_report)
print(f'\nReport saved to: {report_path}')
print(f'Report length:   {len(full_report):,} characters')
print('\n' + '=' * 70)
print('  COPY EVERYTHING ABOVE THIS LINE AND PASTE INTO THE NEW CHAT')
print('=' * 70)