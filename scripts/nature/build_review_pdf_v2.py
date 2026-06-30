"""
PAD4-DB v2 — Combined Review PDF (new figure package)
Embeds all 11 figures (6 main + 5 supp) and renders all 10 tables as pages.
Output: outputs/figures/nature_v2/PAD4_DB_v2_REVIEW_v2.pdf
"""
import os, textwrap
import pandas as pd
from PIL import Image
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image as RLImage,
    Table as RLTable, TableStyle, PageBreak, HRFlowable,
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT

ROOT   = "/home/nidhal/PAD4-db_V2"
FIGDIR = os.path.join(ROOT, "outputs/figures/nature_v2")
TABDIR = os.path.join(ROOT, "outputs/tables/nature_v2/csv")
OUT    = os.path.join(FIGDIR, "PAD4_DB_v2_REVIEW_v2.pdf")

PAGE_W, PAGE_H = letter          # 8.5 × 11 in
MARGIN = 0.55 * inch

# ── Styles ──────────────────────────────────────────────────────────────────
styles = getSampleStyleSheet()

def sty(name, parent='Normal', **kw):
    return ParagraphStyle(name, parent=styles[parent], **kw)

ST_COVER  = sty('cover',  fontSize=28, textColor=colors.HexColor('#004488'),
                leading=36, alignment=TA_CENTER, spaceAfter=12)
ST_SUB    = sty('sub',    fontSize=13, textColor=colors.HexColor('#555555'),
                leading=18, alignment=TA_CENTER, spaceAfter=8)
ST_BODY   = sty('body',   fontSize=8.5, textColor=colors.HexColor('#333333'),
                leading=12, spaceAfter=4)
ST_FIGCAP = sty('figcap', fontSize=8, textColor=colors.HexColor('#555555'),
                leading=11, alignment=TA_CENTER, spaceAfter=6)
ST_SECHDG = sty('sechd',  fontSize=12, textColor=colors.HexColor('#004488'),
                leading=16, spaceBefore=12, spaceAfter=6,
                fontName='Helvetica-Bold')
ST_TABCAP = sty('tabcap', fontSize=9, textColor=colors.HexColor('#333333'),
                leading=13, spaceBefore=8, spaceAfter=4,
                fontName='Helvetica-Bold')
ST_CELL   = sty('cell',   fontSize=7.5, textColor=colors.HexColor('#222222'),
                leading=10)
ST_HDR    = sty('hdr',    fontSize=8,   textColor=colors.white,
                leading=11, alignment=TA_CENTER, fontName='Helvetica-Bold')
ST_NOTE   = sty('note',   fontSize=7,   textColor=colors.HexColor('#888888'),
                leading=10, spaceAfter=4)

NAVY  = colors.HexColor('#004488')
BLUE  = colors.HexColor('#0077BB')
LGREY = colors.HexColor('#E8E8E8')
MGREY = colors.HexColor('#CCCCCC')
DGREY = colors.HexColor('#555555')
TEAL  = colors.HexColor('#009988')
RED   = colors.HexColor('#CC3311')
ORNG  = colors.HexColor('#EE7733')

# ── Helpers ──────────────────────────────────────────────────────────────────
def embed_figure(path, caption, max_w=None, max_h=None):
    """Return [RLImage, caption_para, PageBreak] for a figure PNG."""
    usable_w = PAGE_W - 2 * MARGIN
    usable_h = PAGE_H - 2 * MARGIN - 1.2 * inch   # leave room for caption

    if max_w: usable_w = min(usable_w, max_w)
    if max_h: usable_h = min(usable_h, max_h)

    with Image.open(path) as im:
        iw, ih = im.size
    ratio = min(usable_w / iw, usable_h / ih)
    draw_w = iw * ratio
    draw_h = ih * ratio
    img = RLImage(path, width=draw_w, height=draw_h)
    cap = Paragraph(caption, ST_FIGCAP)
    return [img, Spacer(1, 6), cap, PageBreak()]

def csv_table(path, title, note=None, max_rows=60, col_widths=None):
    """Return flowables for a table from CSV."""
    df = pd.read_csv(path)
    if len(df) > max_rows:
        df_show = df.head(max_rows)
        truncated = True
    else:
        df_show = df
        truncated = False

    out = []
    out.append(Paragraph(title, ST_TABCAP))

    # Build data list
    headers = [Paragraph(str(c).replace('_', ' ').title(), ST_HDR) for c in df_show.columns]
    rows = [headers]
    for _, row in df_show.iterrows():
        cells = []
        for v in row.values:
            txt = str(v) if not pd.isna(v) else '—'
            if len(txt) > 40:
                txt = txt[:38] + '…'
            cells.append(Paragraph(txt, ST_CELL))
        rows.append(cells)

    ncols = len(df_show.columns)
    usable_w = PAGE_W - 2 * MARGIN
    if col_widths is None:
        cw = [usable_w / ncols] * ncols
    else:
        cw = col_widths

    tbl = RLTable(rows, colWidths=cw, repeatRows=1)
    tbl.setStyle(TableStyle([
        ('BACKGROUND',  (0, 0), (-1, 0),  NAVY),
        ('TEXTCOLOR',   (0, 0), (-1, 0),  colors.white),
        ('FONTNAME',    (0, 0), (-1, 0),  'Helvetica-Bold'),
        ('FONTSIZE',    (0, 0), (-1, 0),  8),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, LGREY]),
        ('FONTNAME',    (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE',    (0, 1), (-1, -1), 7.5),
        ('GRID',        (0, 0), (-1, -1), 0.3, MGREY),
        ('TOPPADDING',  (0, 0), (-1, -1), 3),
        ('BOTTOMPADDING',(0,0), (-1, -1), 3),
        ('LEFTPADDING', (0, 0), (-1, -1), 4),
        ('RIGHTPADDING',(0, 0), (-1, -1), 4),
        ('VALIGN',      (0, 0), (-1, -1), 'MIDDLE'),
    ]))
    out.append(tbl)
    if truncated:
        out.append(Paragraph(
            f'(Showing first {max_rows} of {len(df)} rows — full data in CSV)',
            ST_NOTE))
    if note:
        out.append(Paragraph(note, ST_NOTE))
    out.append(Spacer(1, 10))
    return out

# ── Assemble document ────────────────────────────────────────────────────────
story = []

# ── Cover page ───────────────────────────────────────────────────────────────
story.append(Spacer(1, 1.8 * inch))
story.append(Paragraph('PAD4-DB v2', ST_COVER))
story.append(Paragraph('New Figure & Table Package — Review Copy', ST_SUB))
story.append(Spacer(1, 0.15 * inch))
story.append(HRFlowable(width="80%", thickness=1.5, color=NAVY, spaceAfter=16))
story.append(Spacer(1, 0.15 * inch))

cover_items = [
    ('6 Main-text figures', 'Fig 1–6 (overhaul per reviewer consensus)'),
    ('5 Supplementary figures', 'Supp S1–S5'),
    ('10 Data tables', '4 main text + 6 supplementary (CSV + LaTeX)'),
    ('Canonical numbers', '16/16 validated'),
    ('Color standard', 'Tol Vibrant palette, colorblind-safe'),
    ('Resolution', '300 DPI PNG + vector PDF per figure'),
]
tbl_data = [[Paragraph(f'<b>{k}</b>', ST_BODY), Paragraph(v, ST_BODY)]
            for k, v in cover_items]
cov_tbl = RLTable(tbl_data, colWidths=[2.2*inch, 4.2*inch])
cov_tbl.setStyle(TableStyle([
    ('ROWBACKGROUNDS', (0, 0), (-1, -1), [colors.white, LGREY]),
    ('TOPPADDING',    (0, 0), (-1, -1), 6),
    ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
    ('LEFTPADDING',   (0, 0), (-1, -1), 8),
    ('RIGHTPADDING',  (0, 0), (-1, -1), 8),
    ('GRID',          (0, 0), (-1, -1), 0.3, MGREY),
]))
story.append(cov_tbl)
story.append(Spacer(1, 0.5 * inch))
story.append(Paragraph(
    'Generated 2026-06-16 · scripts/nature/ · PAD4-DB v2 pipeline',
    sty('gen', fontSize=7.5, textColor=DGREY, alignment=TA_CENTER)))
story.append(PageBreak())

# ── Main figures ─────────────────────────────────────────────────────────────
main_figs = [
    ('fig01_headline.png',
     'Figure 1 — PAD4 Inhibitor Landscape. '
     '(a) t-SNE chemical space colored by scaffold family. '
     '(b) Top-5 scaffold families with RDKit 2D depictions. '
     '(c) pIC50 density by scaffold class. '
     '(d) PAD4-DB v2 key statistics at a glance.'),
    ('fig02_source_overlap.png',
     'Figure 2 — Source Database Overlap. '
     '(a) UpSet-style intersection plot (n=3,093 total; 7 combinations). '
     '(b) Source coverage totals; PubChem 91.2%, BindingDB 91.4%, ChEMBL 50.6%.'),
    ('fig03_potency.png',
     'Figure 3 — pIC50 Distribution. '
     '(a) Global histogram + KDE (median=6.84). '
     '(b) Multi-source vs PubChem-only split. '
     '(c) Violin by assay mechanism (enzymatic/confirmed/FP/covalent). '
     '(d) Cliff degree vs pIC50 (hub A=★ navy; hub B=◆ red).'),
    ('fig04_scaffold.png',
     'Figure 4 — Scaffold Landscape. '
     '(a) Top-30 series bar chart (largest n=174). '
     '(b) Series size distribution (375 series, log scale). '
     '(c) Lorenz concentration curve (Gini=0.532). '
     '(d) t-SNE colored by scaffold membership.'),
    ('fig05_cliff_network.png',
     'Figure 5 — Severe Cliff Network. '
     '(a) 99 nodes / 94 edges; hub A (★ navy, degree 15/12) and hub B (◆ red, degree 12/11). '
     '(b) Top-12 compounds by severe cliff degree.'),
    ('fig06_mmp.png',
     'Figure 6 — MMP Analysis. '
     '(a) Pairs by cliff tier and R-group change type (85/94 severe cliffs MMP-validated). '
     '(b) Discontinuity score vs pIC50. '
     '(c) ΔpIC50 for all severe vs MMP-validated. '
     '(d) Top-10 shared MMP cores by pair count.'),
]

story.append(Paragraph('Main-Text Figures', ST_SECHDG))
story.append(HRFlowable(width="100%", thickness=0.5, color=BLUE, spaceAfter=10))

for fname, cap in main_figs:
    fpath = os.path.join(FIGDIR, fname)
    if os.path.exists(fpath):
        story.extend(embed_figure(fpath, cap))
    else:
        story.append(Paragraph(f'[MISSING: {fname}]', ST_FIGCAP))
        story.append(PageBreak())

# ── Supplementary figures ─────────────────────────────────────────────────────
supp_figs = [
    ('supp_s01_pipeline.png',
     'Supplementary Figure S1 — Pipeline Workflow. '
     'Five-stage processing pipeline: source databases → standardize → deduplicate → analyse → PAD4-DB v2 (3,093 compounds).'),
    ('supp_s02_sali.png',
     'Supplementary Figure S2 — SALI Landscape. '
     '(a) SALI distribution (log y; SALI>10: n=335, SALI>20: n=19). '
     '(b) SALI vs ΔpIC50 colored by Tanimoto. '
     '(c) Top-20 SALI pairs by tier.'),
    ('supp_s03_patent.png',
     'Supplementary Figure S3 — Patent Scaffold Analysis. '
     't-SNE overlay of 233 patent-exclusive compounds (orange) on published compounds (blue). '
     'Patent chemistry occupies peripheral chemical space; contributes 1/94 severe cliffs.'),
    ('supp_s04_independence.png',
     'Supplementary Figure S4 — Source Independence Scoring. '
     '(a) Score distribution lollipop chart. '
     '(b) Threshold comparison: 528 compounds at ≥0.6 vs 361 at ≥0.7. '
     '(c) pIC50 KDE by independence tier.'),
    ('supp_s05_scaffold_structures.png',
     'Supplementary Figure S5 — Top-20 Scaffold 2D Structures (600 DPI). '
     'RDKit-generated Murcko scaffold depictions for the 20 largest series, '
     'ranked by compound count with pIC50 summary.'),
]

story.append(Paragraph('Supplementary Figures', ST_SECHDG))
story.append(HRFlowable(width="100%", thickness=0.5, color=TEAL, spaceAfter=10))

for fname, cap in supp_figs:
    fpath = os.path.join(FIGDIR, fname)
    if os.path.exists(fpath):
        story.extend(embed_figure(fpath, cap))
    else:
        story.append(Paragraph(f'[MISSING: {fname}]', ST_FIGCAP))
        story.append(PageBreak())

# ── Main tables ───────────────────────────────────────────────────────────────
story.append(Paragraph('Main-Text Tables', ST_SECHDG))
story.append(HRFlowable(width="100%", thickness=0.5, color=RED, spaceAfter=10))

main_tabs = [
    ('table1_hub_summary.csv',
     'Table 1 — Cliff Hub Compound Summary',
     'Hub A: high-degree severe cliffs in rank-1 azaindole-benzimidazole scaffold series. '
     'Hub B: structural singletons (cyclobutyl vs cyclopentyl sulfonamide, Tanimoto=0.975).'),
    ('table2_cliff_summary.csv',
     'Table 2 — Activity Cliff Summary by Tier',
     'Severe threshold: Tanimoto ≥ 0.8 AND |ΔpIC50| ≥ 2.0. '
     'Moderate: ≥ 0.8 AND |ΔpIC50| ≥ 1.5. Broad: ≥ 0.8 AND |ΔpIC50| ≥ 1.0.'),
    ('table3_mmp_summary.csv',
     'Table 3 — MMP Analysis Summary',
     '85/94 (90.4%) of severe cliffs validated by matched molecular pair analysis. '
     'Single-atom changes account for the majority of R-group transformations.'),
    ('table4_source_independence.csv',
     'Table 4 — Source Independence by Combination',
     'Score = 1 − (shared_source_measurements / total_measurements). '
     'Threshold 0.6 yields 528 truly independent compounds; 0.7 yields 361.'),
]

for fname, title, note in main_tabs:
    fpath = os.path.join(TABDIR, fname)
    if os.path.exists(fpath):
        story.extend(csv_table(fpath, title, note=note))
    else:
        story.append(Paragraph(f'[MISSING: {fname}]', ST_NOTE))
story.append(PageBreak())

# ── Supplementary tables ──────────────────────────────────────────────────────
story.append(Paragraph('Supplementary Tables', ST_SECHDG))
story.append(HRFlowable(width="100%", thickness=0.5, color=ORNG, spaceAfter=10))

supp_tabs = [
    ('tableS1_source_coverage.csv',
     'Table S1 — Source Database Coverage',
     None),
    ('tableS2_mechanism_pic50.csv',
     'Table S2 — pIC50 by Assay Mechanism',
     'Kruskal–Wallis p < 0.001. Post-hoc Dunn group letters where available.'),
    ('tableS3_sali_distribution.csv',
     'Table S3 — SALI Distribution by Cliff Tier',
     None),
    ('tableS4_patent_analysis.csv',
     'Table S4 — Patent-Exclusive Compound Analysis',
     None),
    ('tableS5_top20_sali_pairs.csv',
     'Table S5 — Top 20 SALI Pairs',
     'Ranked by SALI score. InChIKeys truncated to 14 characters.'),
    ('tableS6_full_compound_list.csv',
     'Table S6 — Full Compound List (first 50 rows shown; 3,093 total)',
     'Complete 3,093-row table available as CSV for data deposition. '
     'Columns: inchi_key, smiles_std, pIC50, source_list, mechanism_class, patent_flag, source_independence_score.',
     30),  # show only 30 rows for S6
]

for item in supp_tabs:
    fname, title = item[0], item[1]
    note = item[2] if len(item) > 2 else None
    max_r = item[3] if len(item) > 3 else 60
    fpath = os.path.join(TABDIR, fname)
    if os.path.exists(fpath):
        story.extend(csv_table(fpath, title, note=note, max_rows=max_r))
    else:
        story.append(Paragraph(f'[MISSING: {fname}]', ST_NOTE))
story.append(PageBreak())

# ── Back page ─────────────────────────────────────────────────────────────────
story.append(Spacer(1, 2 * inch))
story.append(Paragraph('End of Review Package', ST_COVER))
story.append(Spacer(1, 0.3 * inch))
story.append(Paragraph(
    'PAD4-DB v2 · 3,093 curated PAD4 inhibitors · '
    'scripts/nature/ · outputs/figures/nature_v2/ · outputs/tables/nature_v2/',
    sty('back', fontSize=7.5, textColor=DGREY, alignment=TA_CENTER)))

# ── Build PDF ─────────────────────────────────────────────────────────────────
doc = SimpleDocTemplate(
    OUT,
    pagesize=letter,
    leftMargin=MARGIN, rightMargin=MARGIN,
    topMargin=MARGIN, bottomMargin=MARGIN,
    title='PAD4-DB v2 — Review Package v2',
    author='PAD4-DB pipeline',
)
doc.build(story)
print(f"Written: {OUT}")
import os
sz = os.path.getsize(OUT) / 1e6
print(f"Size: {sz:.1f} MB")
print(f"Pages: approx {len([x for x in story if isinstance(x, PageBreak)]) + 1}")
