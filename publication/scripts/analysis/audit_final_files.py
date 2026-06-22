"""
audit_final_files.py — Verifies the final deliverable file set for PAD4-DB.

Checks (1) existence, (2) canonical numbers recomputed from each data file,
(3) figure PNG staleness vs the script that generates it, (4) which manuscript
is current. Read-only; nothing is modified.

Run: conda run -n pad4bench python publication/scripts/analysis/audit_final_files.py
"""
from pathlib import Path
import os, time
import pandas as pd

ROOT = Path('/home/nidhal/PAD4-db_V2')
OK, WARN, FAIL = '[OK]  ', '[WARN]', '[FAIL]'
issues = []

def status(cond, label, detail='', warn=False):
    tag = OK if cond else (WARN if warn else FAIL)
    if not cond:
        issues.append((warn, label, detail))
    print(f"  {tag} {label}" + (f"  — {detail}" if detail else ''))

def mtime(p): return os.path.getmtime(p) if Path(p).exists() else 0

print("=" * 70)
print("PAD4-DB FINAL FILE AUDIT")
print("=" * 70)

# ── 1. CORE DATA FILES + canonical numbers ────────────────────────────────────
print("\n[1] CORE DATA FILES & CANONICAL NUMBERS")
pc = ROOT / 'publication/data/pad4_compounds.parquet'
status(pc.exists(), 'pad4_compounds.parquet exists')
if pc.exists():
    df = pd.read_parquet(pc)
    status(df.inchi_key.nunique() == 3093, '3,093 unique compounds', f'got {df.inchi_key.nunique()}')
    status('hub_class' in df.columns, "hub_class column present")
    if 'hub_class' in df.columns:
        status(int((df.hub_class != 'none').sum()) == 4, '4 hub compounds tagged',
               f"got {int((df.hub_class!='none').sum())}")

acp = ROOT / 'publication/data/activity_cliffs.parquet'
status(acp.exists(), 'activity_cliffs.parquet exists')
if acp.exists():
    ac = pd.read_parquet(acp)
    sev = ac[ac.cliff_tier == 'severe']
    status(len(sev) == 94, '94 severe cliffs', f'got {len(sev)}')
    status('ecfp4_only_cliff' in ac.columns, 'ecfp4_only_cliff column present')
    status('tanimoto_ecfp6' in ac.columns, 'tanimoto_ecfp6 column present')
    if 'ecfp4_only_cliff' in ac.columns:
        status(int(sev.ecfp4_only_cliff.sum()) == 13, '13 ecfp4_only severe pairs',
               f'got {int(sev.ecfp4_only_cliff.sum())}')

mmp = ROOT / 'outputs/mmp/mmp_pairs_cliff99.csv'
status(mmp.exists(), 'mmp_pairs_cliff99.csv exists')
if mmp.exists():
    m = pd.read_csv(mmp)
    status(len(m) == 707, '707 MMP relationships', f'got {len(m)}')
    status(m.shared_core.nunique() == 24, '24 unique shared cores', f'got {m.shared_core.nunique()}')

# publication/data vs data/processed consistency
a1, a2 = ROOT/'publication/data/activity_cliffs.parquet', ROOT/'data/processed/activity_cliffs.parquet'
if a1.exists() and a2.exists():
    status(os.path.getsize(a1) == os.path.getsize(a2),
           'publication/data activity_cliffs == data/processed copy',
           f'{os.path.getsize(a1)} vs {os.path.getsize(a2)}', warn=True)

# ── 2. FIGURES: existence + staleness vs generating script ────────────────────
print("\n[2] FIGURES — existence and freshness (PNG newer than its script?)")
fig_map = [
    ('fig01_headline', 'figures/fig01_headline.py', 'main'),
    ('fig02_source_overlap', 'figures/fig02_source_overlap.py', 'main'),
    ('fig03_potency', 'figures/fig03_potency.py', 'main'),
    ('fig04_scaffold', 'figures/fig04_scaffold.py', 'main'),
    ('fig05_cliff_network', 'figures/fig05_cliff_network.py', 'main'),
    ('fig06b_cliff_pairs', 'figures/fig06b_cliff_pairs.py', 'main'),
    ('fig_s01_scaffold_cliff_density', 'figures/fig_s01_scaffold_cliff_density.py', 'supplementary'),
    ('fig_s02_assay_enrichment', 'analysis/supp_assay_cliff_enrichment.py', 'supplementary'),
    ('fig_s03_sas_map', 'figures/fig_s03_sas_map.py', 'supplementary'),
]
for fig, script, sub in fig_map:
    png = ROOT / f'publication/figures/{sub}/{fig}.png'
    scr = ROOT / f'publication/scripts/{script}'
    if not png.exists():
        status(False, f'{fig}.png exists'); continue
    if scr.exists() and mtime(png) < mtime(scr):
        status(False, f'{fig}.png is STALE',
               f'script edited {time.ctime(mtime(scr))} after PNG {time.ctime(mtime(png))}', warn=True)
    else:
        status(True, f'{fig}.png fresh (PDF: {(png.with_suffix(".pdf")).exists()})')

# Pre-v7 supplementary figures were archived (not part of the v7 deliverable)
arch = ROOT / 'publication/figures/supplementary/_archive_pre_v7'
n_arch = len(list(arch.glob('*.png'))) if arch.exists() else 0
print(f"\n    Pre-v7 legacy supplementary figures archived to _archive_pre_v7/ ({n_arch} PNGs) — not in deliverable.")
# Confirm the supplementary folder holds EXACTLY the 3 current v7 figures
supp_pngs = sorted(p.name for p in (ROOT/'publication/figures/supplementary').glob('*.png'))
status(supp_pngs == ['fig_s01_scaffold_cliff_density.png', 'fig_s02_assay_enrichment.png',
                     'fig_s03_sas_map.png'],
       'supplementary/ holds exactly S1, S2, S3', f'{supp_pngs}')

# ── 3. MANUSCRIPT — which is current ──────────────────────────────────────────
print("\n[3] MANUSCRIPT FILES")
mans = {
    'PAD4_DB_manuscript_DRAFT_v7.md': 'CURRENT source',
    'PAD4_DB_manuscript_DRAFT_v7.docx': 'CURRENT build',
    'PAD4_DB_manuscript_DRAFT_v7.pdf': 'CURRENT pdf',
}
for f, role in mans.items():
    p = ROOT/'publication/manuscript'/f
    status(p.exists(), f'{f} ({role})', time.ctime(mtime(p)) if p.exists() else '')
# build chain freshness: md -> docx -> pdf
md = ROOT/'publication/manuscript/PAD4_DB_manuscript_DRAFT_v7.md'
dx = ROOT/'publication/manuscript/PAD4_DB_manuscript_DRAFT_v7.docx'
pdf = ROOT/'publication/manuscript/PAD4_DB_manuscript_DRAFT_v7.pdf'
status(mtime(dx) >= mtime(md), 'docx newer than md source', warn=True)
status(mtime(pdf) >= mtime(dx), 'pdf newer than docx', warn=True)
# docx embeds expected figures
if dx.exists():
    import zipfile
    n_img = len([n for n in zipfile.ZipFile(dx).namelist() if n.startswith('word/media/')])
    status(n_img >= 9, f'docx embeds >= 9 figure images', f'got {n_img}')

print("\n    Legacy/superseded manuscripts present (not the deliverable):")
for f in ['PAD4_DB_v2_FINAL.docx', 'PAD4_DB_v2_manuscript.md', 'PAD4_DB_v2_manuscript_integrated.docx']:
    p = ROOT/'publication/manuscript'/f
    if p.exists():
        print(f"  {WARN} {f}  — superseded by DRAFT_v7 ({time.ctime(mtime(p))})")

# ── 4. SUPPLEMENTARY TABLES + analysis scripts ────────────────────────────────
print("\n[4] SUPPLEMENTARY TABLES & ANALYSIS SCRIPTS")
for t in ['supp_statistical_tests.csv', 'supp_scaffold_cliff_density.csv', 'supp_hub_properties.csv',
          'supp_assay_enrichment.csv', 'supp_patent_cliff_odds.csv', 'supp_sas_quadrants.csv']:
    status((ROOT/'outputs/tables'/t).exists(), f'{t} exists')
for s in ['supp_statistical_tests.py', 'supp_hub_properties.py', 'supp_assay_cliff_enrichment.py',
          'supp_patent_analysis.py', 'audit_final_files.py']:
    status((ROOT/'publication/scripts/analysis'/s).exists(), f'{s} exists')
for extra in ['outputs/review/PAD4DB_v2_figures_and_tables_review.pdf',
              '00_validate_canonical_numbers.py']:
    status((ROOT/extra).exists(), f'{extra} exists')

# ── SUMMARY ───────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
fails = [i for i in issues if not i[0]]
warns = [i for i in issues if i[0]]
if not fails and not warns:
    print("RESULT: ALL CLEAR — every final file present, current, and consistent.")
else:
    print(f"RESULT: {len(fails)} FAIL, {len(warns)} WARN")
    for _, lbl, det in fails:
        print(f"  FAIL: {lbl} {('— '+det) if det else ''}")
    for _, lbl, det in warns:
        print(f"  WARN: {lbl} {('— '+det) if det else ''}")
print("=" * 70)
