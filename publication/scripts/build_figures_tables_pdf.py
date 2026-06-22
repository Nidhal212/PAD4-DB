"""
build_figures_tables_pdf.py — generate publication/PAD4_DB_v2_figures_and_tables.pdf
Contains all main + supplementary figures with captions, and all main tables.
Run from project root: conda run -n pad4bench python3 publication/scripts/build_figures_tables_pdf.py
"""
import os, sys, textwrap
os.chdir('/home/nidhal/PAD4-db_V2')

from pathlib import Path
from PIL import Image as PILImage
import io

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
pt = 1.0  # reportlab native unit is already points
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, PageBreak, Image,
    Table, TableStyle, HRFlowable, KeepTogether
)
from reportlab.platypus.flowables import Flowable

# ─── Paths ───────────────────────────────────────────────────────────────────
ROOT   = Path('.')
FIG_M  = ROOT / 'publication/figures/main'
FIG_S  = ROOT / 'publication/figures/supplementary'
OUT    = ROOT / 'publication/PAD4_DB_v2_figures_and_tables.pdf'
OUT.parent.mkdir(parents=True, exist_ok=True)

# ─── Page geometry ───────────────────────────────────────────────────────────
W, H = A4  # 595.3 x 841.9 pt
MARGIN_L = MARGIN_R = 18*mm
MARGIN_T = MARGIN_B = 18*mm
CONTENT_W = W - MARGIN_L - MARGIN_R

# ─── Colours ─────────────────────────────────────────────────────────────────
NAVY   = colors.HexColor('#1A237E')
BLUE   = colors.HexColor('#0072B2')
ORANGE = colors.HexColor('#E69F00')
VERM   = colors.HexColor('#D55E00')
GREY   = colors.HexColor('#555555')
LGREY  = colors.HexColor('#DDDDDD')
BLACK  = colors.black
WHITE  = colors.white

# ─── Styles ──────────────────────────────────────────────────────────────────
styles = getSampleStyleSheet()

def _style(**kw):
    """Create a ParagraphStyle from keyword overrides."""
    s = ParagraphStyle('_', **kw)
    return s

S_TITLE = _style(fontName='Helvetica-Bold', fontSize=16, leading=22,
                 alignment=TA_CENTER, textColor=NAVY, spaceAfter=6)
S_SUBTITLE = _style(fontName='Helvetica', fontSize=10, leading=14,
                    alignment=TA_CENTER, textColor=GREY, spaceAfter=4)
S_SECTION = _style(fontName='Helvetica-Bold', fontSize=12, leading=16,
                   textColor=NAVY, spaceBefore=10, spaceAfter=4,
                   borderPad=0)
S_CAPTION_BOLD = _style(fontName='Helvetica-Bold', fontSize=8, leading=11,
                        textColor=BLACK, spaceAfter=1)
S_CAPTION = _style(fontName='Helvetica', fontSize=8, leading=11,
                   textColor=GREY, spaceAfter=8, alignment=TA_JUSTIFY)
S_FOOTNOTE = _style(fontName='Helvetica-Oblique', fontSize=7, leading=10,
                    textColor=GREY, spaceAfter=4)
S_NORMAL = _style(fontName='Helvetica', fontSize=8, leading=11,
                  textColor=BLACK, spaceAfter=4)
S_TABLE_TITLE = _style(fontName='Helvetica-Bold', fontSize=8.5, leading=12,
                       textColor=NAVY, spaceBefore=8, spaceAfter=3)
S_TABLE_FOOT = _style(fontName='Helvetica-Oblique', fontSize=7, leading=10,
                      textColor=GREY, spaceAfter=10)

# ─── Helper: load PNG, scale to max_width preserving aspect ──────────────────
def figure_flowable(png_path, max_w=CONTENT_W, max_h=None):
    """Return an Image flowable scaled to fit within (max_w, max_h)."""
    if max_h is None:
        max_h = H - MARGIN_T - MARGIN_B - 60  # leave room for caption
    img = PILImage.open(png_path)
    iw, ih = img.size
    scale = min(max_w / iw, max_h / ih)
    return Image(str(png_path), width=iw * scale, height=ih * scale)

# ─── Helper: section divider ─────────────────────────────────────────────────
def section_header(title):
    return [
        HRFlowable(width=CONTENT_W, thickness=1, color=NAVY, spaceAfter=4),
        Paragraph(title, S_SECTION),
    ]

# ─── Helper: figure block (image + caption label + caption text) ──────────────
def figure_block(png_path, label, caption_text, footnote=None, max_h=None):
    """Return a list of flowables: image, label, caption, optional footnote."""
    items = []
    if Path(png_path).exists():
        items.append(figure_flowable(png_path, max_h=max_h))
        items.append(Spacer(1, 3))
    else:
        items.append(Paragraph(f'[Figure not found: {png_path}]', S_FOOTNOTE))
    items.append(Paragraph(label, S_CAPTION_BOLD))
    items.append(Paragraph(caption_text, S_CAPTION))
    if footnote:
        items.append(Paragraph(footnote, S_FOOTNOTE))
    return items

# ─── Table style helper ───────────────────────────────────────────────────────
def make_table(data, col_widths=None, header_bg=NAVY, header_fg=WHITE,
               row_colors=True):
    """Build a styled reportlab Table from list-of-lists data."""
    n_cols = len(data[0])
    if col_widths is None:
        col_widths = [CONTENT_W / n_cols] * n_cols

    t = Table(data, colWidths=col_widths, repeatRows=1)
    cmd = [
        # Header
        ('BACKGROUND', (0, 0), (-1, 0), header_bg),
        ('TEXTCOLOR',  (0, 0), (-1, 0), header_fg),
        ('FONTNAME',   (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE',   (0, 0), (-1, 0), 7.5),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 4),
        ('TOPPADDING',    (0, 0), (-1, 0), 4),
        # Body
        ('FONTNAME',   (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE',   (0, 1), (-1, -1), 7),
        ('TOPPADDING',    (0, 1), (-1, -1), 3),
        ('BOTTOMPADDING', (0, 1), (-1, -1), 3),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        # Grid
        ('GRID', (0, 0), (-1, -1), 0.4, LGREY),
        ('LINEBELOW', (0, 0), (-1, 0), 1.0, NAVY),
        ('LINEBELOW', (0, -1), (-1, -1), 0.6, GREY),
    ]
    if row_colors:
        for i in range(1, len(data)):
            if i % 2 == 0:
                cmd.append(('BACKGROUND', (0, i), (-1, i),
                             colors.HexColor('#F5F5F5')))
    t.setStyle(TableStyle(cmd))
    return t

# ─────────────────────────────────────────────────────────────────────────────
# BUILD STORY
# ─────────────────────────────────────────────────────────────────────────────
story = []

# ══ COVER PAGE ════════════════════════════════════════════════════════════════
story += [
    Spacer(1, 40*mm),
    Paragraph('PAD4-DB v2', S_TITLE),
    Spacer(1, 4),
    Paragraph('A Provenance-First Database of PAD4 Inhibitors with<br/>'
              'Activity Cliff Characterization and Source Independence Scoring',
              _style(fontName='Helvetica-Bold', fontSize=13, leading=18,
                     alignment=TA_CENTER, textColor=NAVY, spaceAfter=8)),
    Spacer(1, 6),
    HRFlowable(width=80*mm, thickness=1.5, color=ORANGE, spaceAfter=6),
    Spacer(1, 6),
    Paragraph('Figures &amp; Tables', _style(fontName='Helvetica-Bold', fontSize=14,
              leading=18, alignment=TA_CENTER, textColor=ORANGE, spaceAfter=8)),
    Spacer(1, 12),
    Paragraph('3,093 curated PAD4 inhibitors · 95 PubChem AIDs · ChEMBL · BindingDB',
              S_SUBTITLE),
    Paragraph('1,244 Bemis-Murcko scaffolds · 94 severe activity cliff pairs',
              S_SUBTITLE),
    Paragraph('Four cliff-hub compounds · Two structural classes · Source independence scoring',
              S_SUBTITLE),
    Spacer(1, 20*mm),
    Paragraph('Generated 2026-06-18', S_SUBTITLE),
    Paragraph('RDKit 2025.09.5 · Python 3.10.19 · conda env: pad4bench',
              S_FOOTNOTE),
    PageBreak(),
]

# ══ CONTENTS ══════════════════════════════════════════════════════════════════
story += [
    Paragraph('Contents', S_SECTION),
    HRFlowable(width=CONTENT_W, thickness=0.5, color=LGREY, spaceAfter=6),
    Spacer(1, 4),
]
contents = [
    ('Main Figures', 'Fig. 1 — t-SNE landscape'),
    ('', 'Fig. 2 — Source overlap'),
    ('', 'Fig. 3 — Potency &amp; composition'),
    ('', 'Fig. 4 — Scaffold diversity'),
    ('', 'Fig. 5 — Activity cliff network'),
    ('', 'Fig. 6 — MMP analysis (panels a + b)'),
    ('Main Tables', 'Table 1 — Cliff-hub properties'),
    ('', 'Table 2 — Cliff tier summary'),
    ('', 'Table 3 — Source distribution'),
    ('', 'Table 4 — Source independence scores'),
    ('', 'Table 5 — Comparison with existing resources'),
    ('Supplementary Figures', 'Fig. S1 — Pipeline flowchart'),
    ('', 'Fig. S2 — SALI distribution'),
    ('', 'Fig. S3 — Patent-exclusive analysis'),
    ('', 'Fig. S4 — Reference compound recovery'),
    ('', 'Fig. S5 — Dominant scaffold structures'),
    ('', 'Fig. S6 — Permutation analysis'),
    ('', 'Fig. S7 — Physicochemical properties'),
    ('Supplementary Tables', 'Tables S1–S6 (descriptions + data)'),
]
for section, item in contents:
    row_txt = (f'<b>{section}</b>&nbsp;&nbsp;&nbsp;{item}' if section
               else f'&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;{item}')
    story.append(Paragraph(row_txt, S_NORMAL))
story.append(PageBreak())

# ══ MAIN FIGURES ══════════════════════════════════════════════════════════════
story += section_header('Main Figures')
story.append(Spacer(1, 4))

# Define figure order: (png filename, label, caption, footnote)
MAIN_FIGS = [
    ('fig01_headline.png',
     'Fig. 1. t-SNE landscape of the PAD4-DB v2 compound space.',
     'Two-dimensional t-SNE embedding of 3,093 compounds based on ECFP4 fingerprints '
     '(perplexity=30). (a) Chemical space coloured by scaffold type; clusters correspond '
     'to dominant scaffold series. (b) Headline statistics for PAD4-DB v2. '
     '(c) pIC50 density by scaffold class; the azaindole-benzimidazole series (blue) '
     'dominates the high-potency region. (d) Top-15 scaffold series by size.',
     None),

    ('fig02_source_overlap.png',
     'Fig. 2. Source overlap and independence scoring.',
     'UpSet plot showing the distribution of 3,093 compounds across PubChem confirmatory, '
     'BindingDB, and ChEMBL sources. (a) Intersection sizes; orange bar marks the '
     '233 patent-exclusive PubChem-only compounds. (b) Per-source compound counts '
     'with percentage of the full dataset.',
     None),

    ('fig03_potency.png',
     'Fig. 3. Potency and composition statistics.',
     '(a) Histogram of consensus pIC50 values (n=3,093; bin width=0.25). Bimodal '
     'distribution with primary peak at pIC50≈7.0 and patent-compound shoulder at '
     'pIC50≈5.0. Mean (6.55) and median (6.84) indicated. (b) pIC50 distribution by '
     'source type (published vs patent-exclusive; Mann-Whitney U p&lt;0.001). '
     '(c) Mechanism class distribution: violin plots of pIC50 by assay class.',
     None),

    ('fig04_scaffold.png',
     'Fig. 4. Scaffold diversity.',
     '(a) Top-30 scaffold series ranked by size; rank-1 azaindole-benzimidazole series '
     'contains 174 compounds. (b) Lorenz curve of scaffold size distribution; '
     'Gini coefficient=0.532 indicates moderate concentration. '
     '(c) Log-scale histogram of series sizes (n=375 series, range 2–174, median 3).',
     None),

    ('fig05_cliff_network.png',
     'Fig. 5. Severe activity cliff network.',
     '(a) Force-directed network of 99 nodes (severe cliff compounds) and 94 edges '
     '(Tanimoto ≥ 0.8, |ΔpIC50| ≥ 2.0). Hub A compounds (orange stars; pIC50≈5.4) '
     'are series-embedded floor compounds in the azaindole-benzimidazole chemotype. '
     'Hub B compounds (green diamonds; pIC50=4.30) are scaffold singletons whose '
     'structural promiscuity generates cross-chemotype cliff pairs. '
     '(b) Severe cliff partner count (degree) for the 12 highest-degree compounds; '
     'four hub compounds account for 50 of 94 severe cliff pairs (53.2%).',
     None),

    ('fig06_mmp.png',
     'Fig. 6a. MMP analysis of severe cliff pairs.',
     '(a) MMP change-type distribution: 80 MMP-validated severe pairs classified as '
     'single-atom change (49), small substituent (28), or medium substituent (8). '
     '(b) Compound-level discontinuity score vs. pIC50; Hub A (orange stars) and '
     'Hub B (green diamonds) identified. (c) ΔpIC50 histogram for all 94 severe pairs '
     'vs. the 80 MMP-validated subset. (d) Top-10 MMP cores by number of cliff pairs '
     'sharing the same molecular core.',
     None),

    ('fig06b_cliff_pairs.png',
     'Fig. 6b. Representative MMP-confirmed severe cliff pairs.',
     'Four severe cliff pairs selected for change-type balance (2 single-atom + '
     '2 small-substituent), warhead-free, non-ecfp4-only, no compound repeated. '
     'Pairs (i)–(ii): single-atom-change pairs involving Hub Class B compounds '
     '(cyclobutyl vs cyclopentyl sulfonamide variants; Tanimoto=0.975). '
     'Pairs (iii)–(iv): small-substituent non-hub pairs. Differing atoms/bonds '
     'highlighted in vermillion; scaffold aligned across each pair. '
     'Labels: ΔpIC50, Tanimoto, MMP type, hub class.',
     None),
]

for fname, label, caption, fn in MAIN_FIGS:
    path = FIG_M / fname
    items = figure_block(str(path), label, caption, fn, max_h=H - MARGIN_T - MARGIN_B - 80)
    story += items
    story.append(PageBreak())

# ══ MAIN TABLES ═══════════════════════════════════════════════════════════════
story += section_header('Main Tables')
story.append(Spacer(1, 6))

# ── Table 1 — Hub compounds ───────────────────────────────────────────────────
story.append(Paragraph('Table 1. Cliff-hub compound properties.', S_TABLE_TITLE))
t1_data = [
    ['Compound ID\n(InChIKey)', 'Class', 'pIC50', 'MW\n(Da)', 'Severe\npairs', '% of 94',
     'In scaffold\nseries', 'Archetype'],
    ['SMADULGDNOCLOP-GISFHXKWSA-N', 'A', '5.390', '611', '15', '16.0%', 'Yes (n=174)',
     'Series-embedded floor'],
    ['RAVBZQAQTVGKIV-XBPDSQQVSA-N', 'A', '5.341', '591', '12', '12.8%', 'Yes (n=174)',
     'Series-embedded floor'],
    ['UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'B', '4.301', '606', '12', '12.8%', 'Singleton',
     'Scaffold-singleton attractor'],
    ['DVCKJOQIVOGXEI-XMMPIXPASA-N', 'B', '4.301', '620', '11', '11.7%', 'Singleton',
     'Scaffold-singleton attractor'],
]
col_w = [90, 25, 32, 30, 32, 30, 48, CONTENT_W - 287]
story.append(make_table(t1_data, col_widths=col_w))
story.append(Paragraph(
    'Class A inter-hub Tanimoto = 0.761; Class B inter-hub Tanimoto = 0.975 '
    '(cyclobutyl vs cyclopentyl; 1 CH₂ difference). Cross-class A–B Tanimoto ≈ 0.49. '
    'Hub compounds collectively account for 50/94 severe cliff pairs (53.2%).',
    S_TABLE_FOOT))
story.append(Spacer(1, 8))

# ── Table 2 — Cliff tier summary ─────────────────────────────────────────────
story.append(Paragraph('Table 2. Activity cliff tier summary.', S_TABLE_TITLE))
t2_data = [
    ['Tier', 'Tanimoto\nthreshold', '|ΔpIC50|\nthreshold', 'Pairs', 'Compounds',
     'MMP-confirmed', 'Max |ΔpIC50|'],
    ['Severe',   '≥ 0.8', '≥ 2.0', '94',  '99',  '80 (85.1%)', '3.045'],
    ['Moderate', '≥ 0.8', '1.5–<2.0', '193', '209', '—',        '1.987'],
    ['Broad',    '≥ 0.8', '1.0–<1.5', '580', '654', '—',        '1.499'],
]
col_w2 = [55, 55, 65, 35, 55, 80, CONTENT_W - 345]
story.append(make_table(t2_data, col_widths=col_w2))
story.append(Paragraph(
    'MMP confirmation: shared-core matched molecular pair analysis using FMCS (RDKit). '
    'Severe cliffs: ECFP4 (radius=2, 2048 bits); 13 pairs (13.8%) flagged as ECFP4-only '
    '(not confirmed by ECFP6 ≥ 0.8 or MMP). Fingerprint-sensitivity analysis: '
    '80% of ECFP4-only borderline pairs (51/64) are MMP-confirmed.',
    S_TABLE_FOOT))
story.append(Spacer(1, 8))

# ── Table 3 — Source distribution ────────────────────────────────────────────
story.append(Paragraph('Table 3. Source distribution summary.', S_TABLE_TITLE))
t3_data = [
    ['Source combination', 'n compounds', 'Independence\nscore', 'Interpretation'],
    ['BindingDB + ChEMBL + PubChem', '1,366 (44.2%)', '0.3', 'Pipeline re-curation'],
    ['BindingDB + PubChem',          '1,199 (38.8%)', '0.5', 'Pipeline re-curation'],
    ['BindingDB + ChEMBL',           '167  (5.4%)',   '0.6', 'Genuinely multi-source'],
    ['ChEMBL + PubChem',             '23   (0.7%)',   '0.7', 'Genuinely multi-source'],
    ['PubChem only (patent-exclusive)','233  (7.5%)', '1.0', 'Single-source'],
    ['BindingDB only',               '95   (3.1%)',   '1.0', 'Single-source'],
    ['ChEMBL only',                  '10   (0.3%)',   '1.0', 'Single-source'],
]
col_w3 = [140, 80, 60, CONTENT_W - 280]
story.append(make_table(t3_data, col_widths=col_w3))
story.append(Paragraph(
    'Source independence score reflects the degree to which multi-database presence '
    'represents genuinely independent experimental replication rather than shared '
    'PubChem–ChEMBL–BindingDB curation pipelines. Scores ≥ 0.6 are free of '
    'pipeline re-curation redundancy (528 compounds, 17.1%).',
    S_TABLE_FOOT))
story.append(Spacer(1, 8))

# ── Table 4 — Independence score summary ─────────────────────────────────────
story.append(Paragraph('Table 4. Source independence score summary.', S_TABLE_TITLE))
t4_data = [
    ['Score', 'n compounds', 'Percentage', 'Interpretation'],
    ['0.3', '1,366', '44.2%', 'BindingDB + ChEMBL + PubChem redundancy'],
    ['0.5', '1,199', '38.8%', 'BindingDB + PubChem redundancy'],
    ['0.6', '167',   '5.4%',  'Genuinely multi-source (BindingDB + ChEMBL)'],
    ['0.7', '23',    '0.7%',  'Genuinely multi-source (ChEMBL + PubChem)'],
    ['1.0', '338',   '10.9%', 'Single-source (no cross-source redundancy)'],
    ['≥ 0.6 (free of re-curation)', '528', '17.1%', 'Free of pipeline redundancy'],
    ['< 0.6 (pipeline redundancy)', '2,565', '82.9%', 'Pipeline-redundant'],
]
col_w4 = [75, 60, 60, CONTENT_W - 195]
# Build table with highlighted summary rows
t4_base = [
    ['Score', 'n compounds', 'Percentage', 'Interpretation'],
    ['0.3', '1,366', '44.2%', 'BindingDB + ChEMBL + PubChem redundancy'],
    ['0.5', '1,199', '38.8%', 'BindingDB + PubChem redundancy'],
    ['0.6', '167',   '5.4%',  'Genuinely multi-source (BindingDB + ChEMBL)'],
    ['0.7', '23',    '0.7%',  'Genuinely multi-source (ChEMBL + PubChem)'],
    ['1.0', '338',   '10.9%', 'Single-source (no cross-source redundancy)'],
    ['≥ 0.6', '528', '17.1%', 'Free of pipeline re-curation redundancy'],
    ['< 0.6', '2,565', '82.9%', 'Pipeline-redundant'],
]
t = Table(t4_base, colWidths=col_w4, repeatRows=1)
t.setStyle(TableStyle([
    ('BACKGROUND',    (0, 0), (-1, 0), NAVY),
    ('TEXTCOLOR',     (0, 0), (-1, 0), WHITE),
    ('FONTNAME',      (0, 0), (-1, 0), 'Helvetica-Bold'),
    ('FONTSIZE',      (0, 0), (-1, -1), 7),
    ('TOPPADDING',    (0, 0), (-1, -1), 3),
    ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
    ('GRID',          (0, 0), (-1, -1), 0.4, LGREY),
    ('LINEBELOW',     (0, 0), (-1, 0), 1.0, NAVY),
    ('BACKGROUND',    (0, 2), (-1, 2), colors.HexColor('#F5F5F5')),
    ('BACKGROUND',    (0, 4), (-1, 4), colors.HexColor('#F5F5F5')),
    ('BACKGROUND',    (0, 6), (-1, 7), colors.HexColor('#FFF8E1')),
    ('FONTNAME',      (0, 6), (-1, 7), 'Helvetica-Bold'),
    ('LINEABOVE',     (0, 6), (-1, 6), 0.8, GREY),
]))
story.append(t)
story.append(Paragraph(
    'Of 528 compounds with score ≥ 0.6: 190 are genuinely multi-source '
    '(BindingDB+ChEMBL n=167, ChEMBL+PubChem n=23); 338 are single-source '
    '(score=1.0, no cross-source redundancy — absence of replication, '
    'not confirmation of independence).',
    S_TABLE_FOOT))
story.append(Spacer(1, 8))

# ── Table 5 — Comparison with existing resources ─────────────────────────────
story.append(Paragraph(
    'Table 5. Comparison with existing PAD4 and activity-cliff resources.',
    S_TABLE_TITLE))
t5_data = [
    ['Feature', 'ChEMBL\n(CHEMBL6111)', 'BindingDB\n(Q9UM07)', 'PAD4-DB v2'],
    ['Multi-source integration\n(PubChem+ChEMBL+BindingDB)', '✗', '✗', '✓'],
    ['Provenance / six-layer architecture', '✗', '✗', '✓'],
    ['Source-independence scoring', '✗', '✗', '✓'],
    ['HTS structural reference (327,336 compounds)', '✗', '✗', '✓'],
    ['Activity-cliff characterisation\n(MMP + permutation-validated)', '✗', '✗', '✓'],
    ['Cliff-hub annotations', '✗', '✗', '✓'],
    ['PAD4 Golden Set (≥2 AIDs, spread ≤0.5)', '✗', '✗', '✓ (n=47)'],
    ['Benchmark-ready splits (via PAD4-Bench)', '✗', '✗', '✓'],
]
col_w5 = [CONTENT_W - 180, 55, 55, 70]
story.append(make_table(t5_data, col_widths=col_w5))
story.append(Paragraph(
    '✓ = feature present; ✗ = not available. '
    'ChEMBL and BindingDB rows reflect the raw target exports without additional curation. '
    'PAD4-DB v2 integrates and curates all three sources.',
    S_TABLE_FOOT))
story.append(PageBreak())

# ══ SUPPLEMENTARY FIGURES ═════════════════════════════════════════════════════
story += section_header('Supplementary Figures')
story.append(Spacer(1, 4))

SUPP_FIGS = [
    ('fig_s01_pipeline.png',
     'Supplementary Fig. S1. PAD4-DB v2 curation pipeline.',
     'Flowchart of the six-step curation pipeline from raw source data to '
     'PAD4-DB v2 outputs. Step 00: raw inventory QC (95 unique AIDs). '
     'Step 01: SMILES standardisation (341,282 rows, 328,976 unique InChIKeys). '
     'Step 02: activity normalisation (99.0% OK). '
     'Steps 03–04: replicate aggregation and deduplication yielding 3,093 unique SAR compounds. '
     'Step 05: SAR analysis (scaffolds, cliffs, MMP, permutation test).',
     None),

    ('fig_s02_sali.png',
     'Supplementary Fig. S2. Structure-Activity Landscape Index (SALI) analysis.',
     '(a) Distribution of SALI values for all pairs at Tanimoto ≥ 0.6 '
     '(n=358,416 pairs); vertical dashed lines mark SALI>10 (n=335) and SALI>20 (n=19). '
     '(b) SALI vs |ΔpIC50| scatter coloured by Tanimoto; high-SALI points (orange) '
     'cluster at high similarity and large potency gap. '
     '(c) Top-20 compound pairs by SALI value with cliff tier annotation.',
     None),

    ('fig_s03_patent.png',
     'Supplementary Fig. S3. Patent-exclusive compound chemical space.',
     't-SNE projection of all 3,093 compounds with patent-exclusive compounds '
     'highlighted (orange; n=233) against published compounds (blue; n=2,860). '
     'Patent-exclusive compounds span multiple chemical space clusters, '
     'contributing 103 unique scaffolds not found in published ChEMBL/BindingDB data. '
     'Mean pIC50 of patent-exclusive compounds = 6.13 vs 6.53 for published.',
     None),

    ('fig_s04_reference_recovery.png',
     'Supplementary Fig. S4. Reference compound recovery.',
     'Concordance scatter (x = published pIC50, y = PAD4-DB v2 consensus pIC50) '
     'for seven recovered PAD4 reference inhibitors. '
     'Five compounds (filled circles) have |ΔpIC50| ≤ 0.15 log units; '
     'mean |ΔpIC50| = 0.061 for these five. '
     'GSK484 (open circle) has |ΔpIC50| = 0.25 (inter-assay variability; '
     'salt-form standardisation to free base). '
     'JBI-589 (orange diamond) has |ΔpIC50| = 0.91, attributed to different '
     'Ca²⁺ concentrations between assay formats [Knuckley2010]. '
     'Grey band = ±0.3 log unit concordance window.',
     None),

    ('fig_s05_scaffold_structures.png',
     'Supplementary Fig. S5. Top-20 scaffold series — 2D structures.',
     'Bemis-Murcko scaffold structures for the 20 largest scaffold series, '
     'rendered using RDKit. Rank-1: azaindole-benzimidazole core (n=174 compounds). '
     'Labels: rank, series size (n), and mean consensus pIC50.',
     None),

    ('fig_s06_permutation.png',
     'Supplementary Fig. S6. Permutation analysis of the activity-cliff landscape.',
     '10,000 permutations of consensus pIC50 values with the Tanimoto similarity '
     'structure held fixed (seed=42). '
     '(a) Null distribution of severe cliff pair counts (grey); '
     'observed count = 94 (vermillion line, far left). '
     'Null mean = 1,923 ± 125; depletion ratio = 0.049; p < 0.0001. '
     'The observed landscape is ~20-fold depleted relative to a random potency assignment, '
     'confirming the cliffs are genuine SAR discontinuities. '
     '(b) Null distribution of hub concentration (fraction of count-matched cliff '
     'pairs incident to the four highest-degree nodes); observed = 53.2% (vermillion, far right). '
     'Null mean = 15.2% ± 2.6%; p < 0.0001. Hub structure is 3.5-fold above chance.',
     None),

    ('fig_s07_physicochemical.png',
     'Supplementary Fig. S7. Physicochemical property landscape.',
     'Eight-panel histogram grid (2×4, constrained layout) of RDKit descriptors for '
     'all 3,093 PAD4-DB v2 compounds (0 parse failures). '
     '(a) Molecular weight (median 590.7 Da; above 500 Da Lipinski cutoff reflecting '
     'the azaindole-benzimidazole scaffold). '
     '(b) Crippen cLogP (median 4.60). '
     '(c) TPSA (median 111.4 Å²). '
     '(d) H-bond acceptors (median 6). '
     '(e) H-bond donors (median 2). '
     '(f) Rotatable bonds (median 6). '
     '(g) Fraction Csp3 (median 0.38). '
     '(h) Aromatic rings (median 5). '
     'Dashed vermillion line = median. Lipinski Ro5: 19.8% strict (0 violations), '
     '68.9% classic (≤1 violation). Veber: 90.7% compliant (RotB≤10, TPSA≤140).',
     None),
]

for fname, label, caption, fn in SUPP_FIGS:
    path = FIG_S / fname
    items = figure_block(str(path), label, caption, fn, max_h=H - MARGIN_T - MARGIN_B - 90)
    story += items
    story.append(PageBreak())

# ══ SUPPLEMENTARY TABLES ══════════════════════════════════════════════════════
story += section_header('Supplementary Tables')
story.append(Spacer(1, 6))
story.append(Paragraph(
    'Full supplementary table data are provided as separate CSV files in '
    '<b>outputs/tables/nature_v2/csv/</b> and as the complete compound list '
    '<b>publication/data/pad4_compounds.parquet</b>. '
    'Descriptions below summarise the content of each table.',
    S_NORMAL))
story.append(Spacer(1, 8))

supp_tables = [
    ('Supplementary Table S1. Source record counts by AID and layer.',
     'Complete AID inventory: 95 unique PubChem AIDs, 1 ChEMBL target (CHEMBL6111), '
     '1 BindingDB target (Q9UM07). For each AID: layer assignment (A–F), raw row count, '
     'potency-space rows (use_in_potency_model=True), and norm_status distribution. '
     'File: outputs/tables/nature_v2/csv/tableS1_source_coverage.csv'),

    ('Supplementary Table S2. Mechanism class pIC50 statistics.',
     '4 mechanism classes: enzymatic (n=2,079), enzymatic_confirmed (n=878), '
     'fp_ic50 (n=115), covalent (n=21). Columns: mean pIC50, median pIC50, SD, min, max. '
     'File: outputs/tables/nature_v2/csv/tableS2_mechanism_pic50.csv'),

    ('Supplementary Table S3. SALI distribution statistics.',
     'SALI statistics stratified by cliff tier (severe/moderate/broad/non-cliff). '
     'Columns: n pairs, SALI>10, SALI>20, max SALI, mean |ΔpIC50|. '
     'File: outputs/tables/nature_v2/csv/tableS3_sali_distribution.csv'),

    ('Supplementary Table S4. Patent-exclusive compound analysis.',
     'Summary metrics for the 233 patent-exclusive (PubChem-only) compounds: '
     'scaffold count, mean/median pIC50, severe cliff contribution (1 pair, 1.06%), '
     'and 103 unique scaffolds not present in published ChEMBL/BindingDB space. '
     'File: outputs/tables/nature_v2/csv/tableS4_patent_analysis.csv'),

    ('Supplementary Table S5. Top-20 SALI pairs.',
     'Top 20 compound pairs ranked by SALI value. '
     'Columns: rank, InChIKey A, InChIKey B, SALI, |ΔpIC50|, Tanimoto, cliff tier. '
     'Maximum SALI = 65.88 (near-identical pair, non-cliff). '
     'File: outputs/tables/nature_v2/csv/tableS5_top20_sali_pairs.csv'),

    ('Supplementary Table S6. Fingerprint sensitivity analysis — 94 severe cliff pairs.',
     'All 94 severe cliff pairs with ECFP4 and ECFP6 Tanimoto values, '
     'MMP validation status (mmp_confirmed), hub involvement (hub_a/hub_b), '
     'and ecfp4_only_cliff flag (True for 13 pairs, 13.8%). '
     '80% of borderline pairs (ECFP4 0.80–0.85, ECFP6 <0.80) are MMP-confirmed (51/64). '
     'File: outputs/tables/nature_v2/table_s6_fingerprint_sensitivity.csv'),
]

for title, desc in supp_tables:
    story.append(Paragraph(title, S_TABLE_TITLE))
    story.append(Paragraph(desc, S_TABLE_FOOT))
    story.append(Spacer(1, 6))

# ── Render S1 source coverage table (first few rows as sample) ───────────────
story.append(Spacer(1, 4))
story.append(Paragraph('Supplementary Table S1 — sample rows:', S_TABLE_TITLE))
import pandas as pd
try:
    s1 = pd.read_csv('outputs/tables/nature_v2/csv/tableS1_source_coverage.csv')
    rows = [list(s1.columns)] + [list(str(v) for v in r) for _, r in s1.iterrows()]
    n_cols = len(rows[0])
    cw = [CONTENT_W / n_cols] * n_cols
    story.append(make_table(rows[:8], col_widths=cw))
    story.append(Paragraph(f'Showing first {min(7, len(s1))} of {len(s1)} rows.', S_TABLE_FOOT))
except Exception as e:
    story.append(Paragraph(f'[Table S1 not available: {e}]', S_FOOTNOTE))

story.append(Spacer(1, 8))

# ── PAD4 Golden Set ───────────────────────────────────────────────────────────
story.append(Paragraph(
    'PAD4 Golden Set (publication/data/PAD4_Golden_Set.csv)',
    S_TABLE_TITLE))
story.append(Paragraph(
    '47 compounds with ≥2 distinct PubChem AIDs and cross-assay pIC50 spread ≤0.5 log units. '
    'Provides a high-confidence internal reproducibility subset. '
    'SHA256: 568075ecb03c59325be0a0dfb10ef3b355ca61209e6b1498286b58989c7d8015.',
    S_TABLE_FOOT))
try:
    gs = pd.read_csv('publication/data/PAD4_Golden_Set.csv')
    gs_show = gs[['inchi_key', 'consensus_pic50', 'n_assays',
                  'max_cross_assay_delta', 'mechanism_class']].head(10)
    gs_rows = [list(gs_show.columns)] + [list(str(v) for v in r)
                                          for _, r in gs_show.iterrows()]
    cw_gs = [120, 60, 45, 70, CONTENT_W - 295]
    story.append(make_table(gs_rows, col_widths=cw_gs))
    story.append(Paragraph(
        f'Showing first 10 of 47 rows. All 47 compounds have n_assays ≥ 2 '
        f'and max_cross_assay_delta ≤ 0.5.',
        S_TABLE_FOOT))
except Exception as e:
    story.append(Paragraph(f'[Golden Set not available: {e}]', S_FOOTNOTE))

# ─── BUILD ───────────────────────────────────────────────────────────────────
def _on_page(canvas, doc):
    """Footer on every page."""
    canvas.saveState()
    canvas.setFont('Helvetica', 6)
    canvas.setFillColor(GREY)
    canvas.drawString(MARGIN_L, 10*mm, 'PAD4-DB v2 — Figures & Tables')
    canvas.drawRightString(W - MARGIN_R, 10*mm, f'Page {doc.page}')
    canvas.restoreState()

doc = SimpleDocTemplate(
    str(OUT),
    pagesize=A4,
    leftMargin=MARGIN_L, rightMargin=MARGIN_R,
    topMargin=MARGIN_T, bottomMargin=MARGIN_B,
)
doc.build(story, onFirstPage=_on_page, onLaterPages=_on_page)
size_kb = OUT.stat().st_size // 1024
print(f"\nDone. Output: {OUT}")
print(f"Size: {size_kb} KB ({size_kb/1024:.1f} MB)")
