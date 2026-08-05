#!/usr/bin/env python
"""Build combined review PDF — all 11 figures + 9 tables."""
import os, warnings
warnings.filterwarnings('ignore')
os.chdir('/home/nidhal/PAD4-db_V2')

import io
import pandas as pd
import numpy as np
from PIL import Image
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image as RLImage,
    Table, TableStyle, PageBreak, HRFlowable, KeepTogether
)
from reportlab.lib.colors import HexColor

W, H = A4          # 595.27 × 841.89 pts
MARGIN = 18 * mm
CW = W - 2 * MARGIN  # content width

NAV    = HexColor('#1A237E')
RED    = HexColor('#CC3311')
BLUE   = HexColor('#0077BB')
ORANGE = HexColor('#EE7733')
GRAY   = HexColor('#555555')
LGRAY  = HexColor('#DDDDDD')
WHITE  = colors.white
BLACK  = colors.black

styles = getSampleStyleSheet()

def S(name='Normal', **kw):
    base = styles.get(name, styles['Normal'])
    return ParagraphStyle(name+'_custom', parent=base, **kw)

title_s   = S(fontSize=18, textColor=NAV, fontName='Helvetica-Bold',
              spaceAfter=4, alignment=TA_CENTER)
subtitle_s = S(fontSize=11, textColor=GRAY, fontName='Helvetica',
               spaceAfter=2, alignment=TA_CENTER)
body_s    = S(fontSize=8,  textColor=BLACK, fontName='Helvetica',
              leading=11,  spaceAfter=4)
caption_s = S(fontSize=8,  textColor=GRAY,  fontName='Helvetica-Oblique',
              spaceAfter=3, alignment=TA_CENTER)
head_s    = S(fontSize=11, textColor=NAV,   fontName='Helvetica-Bold',
              spaceBefore=6, spaceAfter=3)
small_s   = S(fontSize=6.5,textColor=GRAY,  fontName='Helvetica',
              leading=9, alignment=TA_CENTER)

OUT = 'outputs/figures/nature'

# ── helpers ───────────────────────────────────────────────────────────────────
def png_flowable(path, max_w=CW, max_h=H - 2*MARGIN - 40*mm):
    """Embed a PNG respecting aspect ratio within max bounds."""
    img = Image.open(path)
    pw, ph = img.size
    scale = min(max_w / pw, max_h / ph)
    return RLImage(path, width=pw*scale, height=ph*scale)

def section_divider(label):
    return [
        Spacer(1, 4*mm),
        HRFlowable(width=CW, thickness=1.2, color=NAV),
        Paragraph(label, S(fontSize=13, textColor=NAV, fontName='Helvetica-Bold',
                            spaceBefore=2, spaceAfter=2)),
        HRFlowable(width=CW, thickness=0.4, color=LGRAY),
        Spacer(1, 3*mm),
    ]

def make_table(rows, col_widths, header_row=True):
    """Build a reportlab Table with Nature styling."""
    tbl = Table(rows, colWidths=col_widths, repeatRows=1 if header_row else 0)
    n_rows = len(rows)
    n_cols = len(rows[0])
    ts = TableStyle([
        # Header
        ('FONTNAME',      (0,0), (-1,0),  'Helvetica-Bold'),
        ('FONTSIZE',      (0,0), (-1,0),  8),
        ('TEXTCOLOR',     (0,0), (-1,0),  BLACK),
        ('ALIGN',         (0,0), (-1,0),  'CENTER'),
        ('LINEABOVE',     (0,0), (-1,0),  1.2, BLACK),
        ('LINEBELOW',     (0,0), (-1,0),  0.8, BLACK),
        # Body
        ('FONTNAME',      (0,1), (-1,-1), 'Helvetica'),
        ('FONTSIZE',      (0,1), (-1,-1), 7.5),
        ('LEADING',       (0,1), (-1,-1), 10),
        ('TEXTCOLOR',     (0,1), (-1,-1), BLACK),
        ('ALIGN',         (0,1), (-1,-1), 'LEFT'),
        # Bottom border
        ('LINEBELOW',     (0,-1),(-1,-1), 0.8, BLACK),
        # Row padding
        ('TOPPADDING',    (0,0), (-1,-1), 3),
        ('BOTTOMPADDING', (0,0), (-1,-1), 3),
        ('LEFTPADDING',   (0,0), (-1,-1), 4),
        ('RIGHTPADDING',  (0,0), (-1,-1), 4),
        # No vertical lines
        ('GRID',          (0,0), (-1,-1), 0, WHITE),
    ])
    tbl.setStyle(ts)
    return tbl

# ── Build story ───────────────────────────────────────────────────────────────
story = []

# ── COVER PAGE ────────────────────────────────────────────────────────────────
story += [
    Spacer(1, 30*mm),
    Paragraph("PAD4-DB v2", title_s),
    Paragraph("Comprehensive Review Package", subtitle_s),
    Spacer(1, 3*mm),
    HRFlowable(width=CW*0.6, thickness=1.5, color=NAV),
    Spacer(1, 5*mm),
    Paragraph("N = 3,093 curated PAD4 inhibitors · 25 columns · 11 Figures · 9 Tables",
              S(fontSize=10, textColor=GRAY, fontName='Helvetica', alignment=TA_CENTER)),
    Spacer(1, 8*mm),
    Paragraph(
        "<b>Contents:</b><br/>"
        "Figures 1–7 · Pipeline &amp; chemical space (Session 1)<br/>"
        "Figures 8–11 · MMP, SALI, patent, independence (Session 2)<br/>"
        "Tables 1–9 · Source coverage, mechanisms, scaffolds, hubs, MMP, SALI, patent, independence",
        S(fontSize=9, textColor=GRAY, fontName='Helvetica',
          alignment=TA_CENTER, leading=14)),
    Spacer(1, 10*mm),
    HRFlowable(width=CW*0.6, thickness=0.5, color=LGRAY),
    Spacer(1, 5*mm),
    Paragraph("Nature-standard design · 600 dpi · Colorblind-safe palette",
              S(fontSize=8, textColor=LGRAY, fontName='Helvetica-Oblique',
                alignment=TA_CENTER)),
    PageBreak(),
]

# ── FIGURES SECTION ───────────────────────────────────────────────────────────
fig_meta = [
    ('fig1_pipeline',     'Figure 1',  'Pipeline Workflow Diagram',
     '5-stage curation pipeline: SMILES standardization → activity normalization → '
     'HTS extraction → deduplication → scaffold/cliff analysis. '
     '341,328 input records → 3,093 curated inhibitors.'),
    ('fig2_upset',        'Figure 2',  'Source Database Overlap (UpSet)',
     'Intersection sizes for 7 source combinations across PubChem (95 AIDs), '
     'ChEMBL (CHEMBL6111), and BindingDB (Q9UM07). '
     'Orange = PubChem-only (n=233); blue = multi-source combinations.'),
    ('fig3_tsne',         'Figure 3',  't-SNE Chemical Space',
     '2D t-SNE projection of ECFP4 fingerprints (Morgan r=2, 2048 bits). '
     'Panels: (A) source category, (B) pIC50 gradient, (C) assay mechanism, '
     '(D) cliff hub compounds (★ = hubs A1/A2 navy, B1/B2 red).'),
    ('fig4_pic50',        'Figure 4',  'pIC50 Distribution',
     'Panels: (A) global histogram + KDE, (B) PubChem-only vs multi-source KDE + Mann-Whitney U, '
     '(C) violin by mechanism (strip removed for enzymatic n=2,079), '
     '(D) cliff degree vs pIC50 with hub stars.'),
    ('fig5_scaffold',     'Figure 5',  'Scaffold Landscape',
     'Panels: (A) top-30 series ranked bar, (B) series size log histogram, '
     '(C) Lorenz curve (Gini=0.532), (D) t-SNE colored by scaffold membership '
     '(orange = 174-compound hub scaffold series).'),
    ('fig6_similarity',   'Figure 6',  'Similarity Landscape',
     'Panels: (A) Tanimoto histogram (sim≥0.6), (B) ΔpIC50 vs Tanimoto cliff scatter, '
     '(C) SALI histogram (linear x, log y), (D) SALI landscape (SALI vs Tanimoto, colored by tier).'),
    ('fig7_network',      'Figure 7',  'Severe Cliff Network + Degree Bar',
     '99 nodes, 94 severe cliff edges. Hub A (navy ★), Hub B (red ★). '
     'Cross-mechanism edges shown as dashed. '
     'Degree bar: top 20 compounds by severe cliff partner count.'),
    ('fig8_mmp',          'Figure 8',  'MMP Analysis',
     'Panels: (A) MMP pairs by cliff tier and type (stacked bar), '
     '(B) discontinuity score vs pIC50 (hub stars), '
     '(C) ΔpIC50: Tanimoto (n=94) vs MMP-validated (n=85) severe cliffs, '
     '(D) top-10 shared MMP cores by frequency.'),
    ('fig9_sali',         'Figure 9',  'SALI Analysis (Supplementary)',
     'Panels: (A) top-20 SALI pairs (max=65.88), '
     '(B) SALI vs ΔpIC50 colored by Tanimoto (viridis), SALI>20 highlighted, '
     '(C) cumulative ECDF with reference lines at SALI=10 (335 pairs) and SALI=20 (19 pairs).'),
    ('fig10_patent',      'Figure 10', 'Patent Scaffold Analysis (Supplementary)',
     'Panels: (A) pIC50 KDE patent-exclusive (n=233, mean=6.082) vs published '
     '(n=2,860, mean=6.588) with Mann-Whitney U, '
     '(B) scaffold size categories side-by-side, '
     '(C) t-SNE overlay: patent (orange) on published (blue).'),
    ('fig11_independence','Figure 11', 'Source Independence (Supplementary)',
     'Panels: (A) independence score histogram (log y), threshold=0.6, '
     '(B) dot plot by source combination (size ∝ log n), '
     '(C) pIC50 KDE: true multi-source (n=528) vs pipeline redundancy (n=2,565).'),
]

story += section_divider("FIGURES  (1 – 11)")

for fname, fig_num, fig_title, caption in fig_meta:
    png_path = f'{OUT}/{fname}.png'
    if not os.path.exists(png_path):
        story.append(Paragraph(f"⚠ Missing: {png_path}", body_s))
        story.append(PageBreak())
        continue

    story.append(Paragraph(f"<b>{fig_num}</b> — {fig_title}", head_s))
    story.append(png_flowable(png_path))
    story.append(Spacer(1, 2*mm))
    story.append(Paragraph(caption, caption_s))
    story.append(PageBreak())

# ── TABLES SECTION ────────────────────────────────────────────────────────────
story += section_divider("TABLES  (1 – 9)")

df     = pd.read_parquet('data/processed/pad4_compounds.parquet')

pairs  = pd.read_parquet('data/processed/activity_pairs_with_sali.parquet')
mmp    = pd.read_csv('outputs/mmp/mmp_pairs_cliff99.csv')

# helper
def R(*args):
    """Row of Paragraph cells for a table."""
    return [Paragraph(str(a), body_s) for a in args]
def RB(*args):
    return [Paragraph(str(a), S(fontName='Helvetica-Bold', fontSize=8)) for a in args]

# ── Table 1: Source coverage ──────────────────────────────────────────────────
story.append(Paragraph("<b>Table 1</b> — Source database coverage (PAD4-DB v2, N=3,093)", head_s))
sl = df['source_list']
def src_row(mask, label):
    sub = df.loc[mask, 'pic50_consensus']
    n = mask.sum()
    return R(label, f'{n:,}', f'{n/3093*100:.1f}%',
             f'{sub.mean():.2f}', f'{sub.median():.2f}', f'{sub.std():.2f}')

t1_data = [
    RB('Source','N compounds','% Dataset','Mean pIC50','Median pIC50','SD pIC50'),
    src_row(sl.str.contains('pubchem_confirmatory'), 'PubChem (≥1 AID)'),
    src_row(sl.str.contains('chembl'),               'ChEMBL'),
    src_row(sl.str.contains('bindingdb'),            'BindingDB'),
    src_row(sl=='bindingdb|chembl|pubchem_confirmatory','All three sources'),
    src_row(sl=='pubchem_confirmatory',              'PubChem only'),
]
story.append(make_table(t1_data, [95, 65, 55, 65, 70, 55]))
story.append(Paragraph('Values computed from deduplicated master dataset.', small_s))
story.append(Spacer(1, 6*mm))

# ── Table 2: Mechanism ────────────────────────────────────────────────────────
story.append(Paragraph("<b>Table 2</b> — pIC50 by assay mechanism class", head_s))
mech_map = {
    'enzymatic':           'Enzymatic (BAEE colorimetric)',
    'enzymatic_confirmed': 'Enzymatic, RFMS-confirmed',
    'fp_ic50':             'FP-based IC50',
    'covalent':            'Covalent (assay-flagged)',
}
t2_data = [RB('Assay mechanism','N','Mean pIC50','Median pIC50','SD','Min','Max')]
for mech, label in mech_map.items():
    sub = df.loc[df['mechanism_class']==mech,'pic50_consensus']
    t2_data.append(R(label, f'{len(sub):,}',
                     f'{sub.mean():.2f}', f'{sub.median():.2f}',
                     f'{sub.std():.2f}', f'{sub.min():.2f}', f'{sub.max():.2f}'))
story.append(make_table(t2_data, [150, 40, 60, 65, 40, 40, 40]))
story.append(Spacer(1, 6*mm))

# ── Table 3: Hub Compound Summary ─────────────────────────────────────────────
story.append(Paragraph("<b>Table 3</b> — Cliff hub compound summary", head_s))
t3_data = [
    RB('ID','InChIKey (14 chars)','Class','pIC50','Severe pairs','% of 94','Scaffold type','In scaffold'),
    R('A1','SMADULGDNOCLOP','A','5.39','15','16.0%','Series member','174'),
    R('A2','RAVBZQAQTVGKIV','A','5.34','12','12.8%','Series member','174'),
    R('B1','UDCDEKJNAMHBFH','B','4.30','12','12.8%','Singleton','1'),
    R('B2','DVCKJOQIVOGXEI','B','4.30','11','11.7%','Singleton','1'),
]
tbl3 = make_table(t3_data, [22, 82, 32, 38, 52, 40, 70, 54])
# Left border: Class A rows navy, Class B rows red
tbl3.setStyle(TableStyle([
    ('LEFTPADDING', (0,1),(0,2), 6),
    ('LEFTPADDING', (0,3),(0,4), 6),
    ('LINEAFTER',   (0,1),(0,2), 3, HexColor('#1A237E')),
    ('LINEAFTER',   (0,3),(0,4), 3, HexColor('#CC3311')),
]))
story.append(tbl3)
story.append(Paragraph(
    'Class A: 174-compound azaindole-benzimidazole series (series-embedded hubs). '
    'Class B: unique Murcko scaffolds (scaffold singletons). '
    'Severe threshold: Tanimoto ≥ 0.8, ΔpIC50 ≥ 2.0. '
    'Combined: 50/94 severe pairs (53.2%).', small_s))
story.append(Spacer(1, 6*mm))

# ── Table 4: Cliff Summary ────────────────────────────────────────────────────
story.append(Paragraph("<b>Table 4</b> — Activity cliff summary", head_s))
t4_data = [RB('Tier','Threshold','N pairs','% of 867','Compounds involved')]
tier_def = [
    ('Severe',   'Tan≥0.8, ΔpIC50≥2.0', 94,  '10.8%', 99),
    ('Moderate', 'Tan≥0.8, ΔpIC50≥1.5', 193, '22.3%', '—'),
    ('Broad',    'Tan≥0.8, ΔpIC50≥1.0', 580, '66.9%', '—'),
    ('Total',    'All tiers combined',   867, '100%',  '—'),
]
for t, d, n, p, c in tier_def:
    t4_data.append(R(t, d, str(n), p, str(c)))
story.append(make_table(t4_data, [55, 130, 55, 60, 105]))
story.append(Spacer(1, 6*mm))

# ── Table 5: SALI Distribution ────────────────────────────────────────────────
story.append(Paragraph("<b>Table 5</b> — SALI distribution summary", head_s))
t5_data = [RB('Tier','SALI pairs','Mean SALI','Median SALI','Max SALI')]
for tier in ['severe','moderate','broad']:
    sub_p = pairs[pairs['cliff_tier']==tier]['sali'].dropna()
    t5_data.append(R(tier.capitalize(), f'{len(sub_p):,}',
                     f'{sub_p.mean():.2f}', f'{sub_p.median():.2f}', f'{sub_p.max():.2f}'))
all_s = pairs['sali'].dropna()
t5_data.append(R('All pairs (sim≥0.6)', f'{len(all_s):,}',
                  f'{all_s.mean():.2f}', f'{all_s.median():.2f}', f'{all_s.max():.2f}'))
t5_data.append(R('SALI > 10', '335', '—', '—', '—'))
t5_data.append(R('SALI > 20', '19',  '—', '—', '65.88'))
story.append(make_table(t5_data, [130, 70, 70, 75, 65]))
story.append(Paragraph('SALI = |ΔpIC50| / (1 − Tanimoto). Pairs with Tanimoto = 1.0 excluded (NaN).', small_s))
story.append(Spacer(1, 6*mm))

# ── Table 6: MMP Summary ──────────────────────────────────────────────────────
story.append(Paragraph("<b>Table 6</b> — MMP analysis summary", head_s))
t6_data = [RB('Category','Metric','Value','Notes')]
t6_rows = [
    ('Overview', 'Total MMP pairs',               '707',             'All 99 cliff compounds'),
    ('Overview', 'MMP-validated severe cliffs',    '85 of 94 (90.4%)','Tanimoto ≥ 0.8, ΔpIC50 ≥ 2.0'),
    ('Overview', 'Non-MMP severe cliffs',          '9',               'Scaffold hops / fragment merges'),
    ('MMP types','Single R-group change',          '49 (57.6%)',      'Most common transformation'),
    ('MMP types','Small substituent change',       '28 (32.9%)',      ''),
    ('MMP types','Medium substituent change',      '8  (9.4%)',       ''),
]
for cat, met, val, note in t6_rows:
    t6_data.append(R(cat, met, val, note))
story.append(make_table(t6_data, [70, 145, 90, 105]))
story.append(Spacer(1, 6*mm))

# ── Table 7: Patent Summary ───────────────────────────────────────────────────
story.append(Paragraph("<b>Table 7</b> — Patent compound analysis", head_s))
from scipy.stats import mannwhitneyu
pat_p = df.loc[df['source_list']=='pubchem_confirmatory','pic50_consensus']
pub_p = df.loc[df['source_list']!='pubchem_confirmatory','pic50_consensus']
_, pv  = mannwhitneyu(pat_p, pub_p, alternative='two-sided')
pv_str = 'p < 0.001' if pv < 0.001 else f'p = {pv:.3f}'

t7_data = [RB('Metric','Value')]
for m, v in [
    ('Patent-exclusive compounds',     '233'),
    ('Published compounds',            '2,860'),
    ('Patent mean pIC50',              '6.082'),
    ('Published mean pIC50',           '6.588'),
    ('Difference (published − patent)','0.506 log units'),
    ('Mann-Whitney U p-value',         pv_str),
    ('Patent-exclusive scaffolds',     '103'),
    ('Patent cliff contribution',      '1 of 94 pairs (1.1%)'),
]:
    t7_data.append(R(m, v))
story.append(make_table(t7_data, [250, 160]))
story.append(Spacer(1, 6*mm))

# ── Table 8: Source Independence ──────────────────────────────────────────────
story.append(Paragraph("<b>Table 8</b> — Source independence by combination", head_s))
t8_data = [RB('Source combination','Score','N','is_multi (≥0.6)','Interpretation')]
combo_info = [
    ('BindingDB only',              1.0,  95,  '✓','Source-exclusive'),
    ('PubChem only',                1.0, 233,  '✓','Source-exclusive'),
    ('ChEMBL only',                 1.0,  10,  '✓','Source-exclusive'),
    ('ChEMBL + PubChem',            0.7,  23,  '✓','Likely independent'),
    ('BindingDB + ChEMBL',          0.6, 167,  '✓','Likely independent'),
    ('BindingDB + PubChem',         0.5,1199,  '✗','Partial redundancy'),
    ('BindingDB + ChEMBL + PubChem',0.3,1366,  '✗','High redundancy'),
]
for src, sc, n, flag, interp in combo_info:
    t8_data.append(R(src, str(sc), f'{n:,}', flag, interp))
story.append(make_table(t8_data, [145, 38, 45, 68, 110]))
story.append(Paragraph(
    'Threshold = 0.6: 528 true multi-source, 2,565 pipeline redundancy. '
    'Single-source score = 1.0 by definition (no cross-source comparison available).', small_s))
story.append(Spacer(1, 6*mm))

# ── Table 9: Top 20 SALI Pairs ────────────────────────────────────────────────
story.append(PageBreak())
story.append(Paragraph("<b>Table 9</b> — Top 20 SALI pairs (deduplicated)", head_s))
p_sali = pairs[pairs['sali'].notna()].copy()
p_sali['pair_key'] = p_sali.apply(
    lambda r: '|'.join(sorted([r['inchi_key_a'],r['inchi_key_b']])), axis=1)
top20 = p_sali.drop_duplicates('pair_key').nlargest(20,'sali').reset_index(drop=True)
t9_data = [RB('Rank','Compound A (14)','Compound B (14)','Tanimoto','ΔpIC50','SALI','Tier')]
for i, row in top20.iterrows():
    t9_data.append(R(
        i+1,
        row['inchi_key_a'][:14],
        row['inchi_key_b'][:14],
        f"{row['tanimoto']:.3f}",
        f"{row['delta_pic50']:.3f}",
        f"{row['sali']:.2f}",
        row['cliff_tier'].capitalize() if pd.notna(row['cliff_tier']) else 'Non-cliff',
    ))
tbl9 = make_table(t9_data, [28, 82, 82, 52, 45, 42, 50])
# Color tier cells
tier_bg = {'Severe': HexColor('#FFE5E0'), 'Moderate': HexColor('#FFF3E0'),
           'Broad': HexColor('#E3F2FD')}
for i, row in top20.iterrows():
    tier = row['cliff_tier'].capitalize() if pd.notna(row['cliff_tier']) else ''
    bg = tier_bg.get(tier)
    if bg:
        tbl9.setStyle(TableStyle([('BACKGROUND', (6,i+1),(6,i+1), bg)]))
story.append(tbl9)
story.append(Paragraph('Pairs are unique (A/B order normalized). SALI = |ΔpIC50| / (1 − Tanimoto).', small_s))

# ── Build PDF ─────────────────────────────────────────────────────────────────
pdf_path = 'outputs/figures/nature/PAD4_DB_v2_FULL_REVIEW.pdf'
doc = SimpleDocTemplate(
    pdf_path, pagesize=A4,
    leftMargin=MARGIN, rightMargin=MARGIN,
    topMargin=MARGIN, bottomMargin=MARGIN,
    title='PAD4-DB v2 — Full Review Package',
    author='PAD4-DB pipeline',
)
doc.build(story)
sz = os.path.getsize(pdf_path) / (1024*1024)
print(f'Saved: {pdf_path}')
print(f'Size:  {sz:.2f} MB')

# Count pages roughly
from reportlab.lib.pagesizes import A4
print('Contents: Cover + 11 figure pages + 9 table pages + section dividers')
print('DONE')
