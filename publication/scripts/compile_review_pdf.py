"""
compile_review_pdf.py  —  Compiles all figures + tables into a single review PDF.

Output: outputs/review/PAD4DB_v2_figures_and_tables_review.pdf
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages
from textwrap import wrap

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT  = ROOT / 'outputs/review'
OUT.mkdir(parents=True, exist_ok=True)
PDF  = OUT / 'PAD4DB_v2_figures_and_tables_review.pdf'

# ─────────────────────────────────────────────────────────────────────────────
# Figure registry — (filename, caption)
# ─────────────────────────────────────────────────────────────────────────────
MAIN_FIGURES = [
    ('fig01_headline.png',
     'Figure 1. Chemical space landscape of 3,093 PAD4 inhibitors (t-SNE). '
     'Background: all compounds (grey, n=3,093). Coloured circles: severe cliff '
     'compounds (n=99, viridis scale by pIC50). Stars: Hub A (navy, n=2, series floor). '
     'Diamonds: Hub B (red, n=2, singleton attractor).'),
    ('fig02_source_overlap.png',
     'Figure 2. Source overlap and independence. (a) UpSet plot by source combination. '
     '(b) Source coverage totals. (c) Stacked bar: 528 non-redundant (17.1%) vs '
     '2,565 pipeline-redundant (82.9%) compounds.'),
    ('fig03_potency.png',
     'Figure 3. Potency distribution. (a) Global histogram + KDE; mean=6.55, median=6.84. '
     '(b) Violin plots by source (Mann-Whitney p < 0.001, PubChem vs Patent). '
     '(c) Compound counts by mechanism class with mean pIC50.'),
    ('fig04_scaffold.png',
     'Figure 4. Scaffold landscape. (a) Top-20 scaffold series by size. '
     '(b) Cliff density scatter (severe cliff pairs / possible pairs within scaffold; '
     'series with n≥4 only). (c) Lorenz curve of scaffold size distribution (Gini=0.532).'),
    ('fig05_cliff_network.png',
     'Figure 5. Severe activity cliff network. (a) 99-node network; node size ∝ degree; '
     'node colour = pIC50 (viridis); edge colour = ΔpIC50 (grey→dark magenta). '
     'Hub A (navy ★), Hub B (red ◆). (b) Top-12 degree bar chart.'),
    ('fig06b_cliff_pairs.png',
     'Figure 6. MMP analysis. (a) MMP-confirmed severe cliff pairs by change type '
     '(80 total: 45 single-atom, 27 small substituent, 8 medium). '
     '(b) Four representative cliff pairs with colourblind-safe diff-atom highlighting: '
     'blue = gain-of-potency atoms, orange = loss-of-potency atoms.'),
]

SUPP_FIGURES = [
    ('fig_s01_scaffold_cliff_density.png',
     'Figure S1. Scaffold cliff-density ranking — which PAD4 chemotypes produce rugged SAR. '
     '(a) The 11/155 scaffold series (n>=4) with any within-scaffold cliff. '
     '(b) Ruggedness is scaffold-intrinsic: S1 (n=174, rugged) vs S2 (n=102, smooth) at matched sampling.'),
    ('fig_s02_assay_enrichment.png',
     'Figure S2. Assay-class cliff enrichment (Fisher exact; no class significantly enriched — '
     'negative result arguing against an assay-format origin for cliffs).'),
    ('fig_s03_sas_map.png',
     'Figure S3. Structure-Activity Similarity (SAS) map. (a) Pairwise Tanimoto '
     'similarity vs |ΔpIC50| for all 358,416 related pairs (hexbin, log density); '
     '94 severe cliffs highlighted in the upper-right activity-cliff quadrant. '
     '(b) |ΔpIC50| distribution by similarity bin; cliff rate annotated per bin '
     '(diagonal absence: only 0.61% of near-identical pairs exceed the cliff threshold).'),
]

# ─────────────────────────────────────────────────────────────────────────────
# Table registry — (title, csv_path, notes)
# ─────────────────────────────────────────────────────────────────────────────
TABLES = [
    ('Table 1. Hub compound summary',
     ROOT / 'outputs/tables/nature_v2/csv/table1_hub_summary.csv',
     '4 hub compounds; A=series floor, B=singleton attractor.'),

    ('Table 2. Activity cliff tier summary',
     ROOT / 'outputs/tables/nature_v2/csv/table2_cliff_summary.csv',
     'Severe: Tan≥0.8, |ΔpIC50|≥2.0. Moderate: ≥1.5. Broad: ≥1.0.'),

    ('Table 3. MMP change-type breakdown',
     ROOT / 'outputs/tables/nature_v2/csv/table3_mmp_summary.csv',
     'Note: "In severe cliffs" from table generation script (85); '
     'bar chart uses canonical pair-key join (80). Discrepancy logged for review.'),

    ('Table 4. Source independence by combination',
     ROOT / 'outputs/tables/nature_v2/csv/table4_source_independence.csv',
     'Score ≥0.6 = non-redundant. 528 / 3,093 (17.1%) meet threshold.'),

    ('Supplementary Table S1. Source coverage',
     ROOT / 'outputs/tables/nature_v2/csv/tableS1_source_coverage.csv', ''),

    ('Supplementary Table S2. Mechanism class pIC50',
     ROOT / 'outputs/tables/nature_v2/csv/tableS2_mechanism_pic50.csv', ''),

    ('Supplementary Table S3. SALI distribution',
     ROOT / 'outputs/tables/nature_v2/csv/tableS3_sali_distribution.csv', ''),

    ('Supplementary Table S4. Patent analysis',
     ROOT / 'outputs/tables/nature_v2/csv/tableS4_patent_analysis.csv', ''),

    ('Supplementary Table S5. Top SALI pairs (deduplicated)',
     ROOT / 'outputs/tables/nature_v2/csv/tableS5_top20_sali_pairs.csv',
     'Deduplicated in this revision: 20 rows → 17 rows (3 duplicate pairs removed).'),

    ('Supplementary Table: Hub physicochemical properties',
     ROOT / 'outputs/tables/supp_hub_properties.csv',
     'Hub (n=4) vs non-hub cliff (n=95). Only pIC50 is significantly different (p=0.007).'),

    ('Supplementary Table: Assay class cliff enrichment',
     ROOT / 'outputs/tables/supp_assay_enrichment.csv',
     "Fisher's exact test; no mechanism class significantly enriched in cliffs."),

    ('Supplementary Table: Patent cliff odds ratio',
     ROOT / 'outputs/tables/supp_patent_cliff_odds.csv',
     'OR=0.121, p=0.006 (patent depleted in cliffs); phi=0.045 (negligible effect). Exploratory.'),

    ('Supplementary Table: SAS quadrant distribution (NEW)',
     ROOT / 'outputs/tables/supp_sas_quadrants.csv',
     'Companion to Fig S10. 96.1% of related pairs are non-descript; '
     'activity cliffs are 0.026% of all 358,416 related pairs.'),
]


def section_page(pdf, title, subtitle='', color='#1A237E'):
    """Full-page section divider."""
    fig, ax = plt.subplots(figsize=(11.69, 8.27))
    ax.set_facecolor(color)
    fig.patch.set_facecolor(color)
    ax.text(0.5, 0.58, title, transform=ax.transAxes,
            fontsize=26, fontweight='bold', color='white',
            ha='center', va='center', wrap=True)
    if subtitle:
        ax.text(0.5, 0.42, subtitle, transform=ax.transAxes,
                fontsize=13, color='#DDDDDD',
                ha='center', va='center')
    ax.axis('off')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def figure_page(pdf, img_path, caption):
    """One figure per page (A4 landscape)."""
    fig, axes = plt.subplots(2, 1, figsize=(11.69, 8.27),
                              gridspec_kw={'height_ratios': [0.88, 0.12]})
    ax_img, ax_cap = axes

    img = mpimg.imread(img_path)
    ax_img.imshow(img)
    ax_img.axis('off')

    wrapped = '\n'.join(wrap(caption, width=140))
    ax_cap.text(0.02, 0.85, wrapped, transform=ax_cap.transAxes,
                fontsize=7.5, va='top', ha='left',
                color='#222222', wrap=False)
    ax_cap.axis('off')
    ax_cap.set_facecolor('#F8F8F8')
    fig.patch.set_facecolor('white')
    plt.tight_layout(pad=0.3)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def table_page(pdf, title, csv_path, notes=''):
    """Render a CSV table as a formatted matplotlib table page."""
    df = pd.read_csv(csv_path)

    # Single wide row (e.g. odds-ratio summary) reads far better transposed
    if len(df) == 1 and len(df.columns) > 5:
        df = df.T.reset_index()
        df.columns = ['Metric', 'Value']

    # Truncate very wide dataframes to fit page
    MAX_COLS = 8
    if len(df.columns) > MAX_COLS:
        df = df.iloc[:, :MAX_COLS]

    MAX_ROWS = 35
    truncated = len(df) > MAX_ROWS
    if truncated:
        shown = df.head(MAX_ROWS)
    else:
        shown = df

    n_rows = len(shown)
    n_cols = len(shown.columns)
    fig_h  = max(4.5, min(8.27, 1.2 + n_rows * 0.22 + 1.0))
    fig, ax = plt.subplots(figsize=(11.69, fig_h))
    ax.axis('off')
    fig.patch.set_facecolor('white')

    # Title
    ax.text(0.0, 1.00, title,
            transform=ax.transAxes, fontsize=11, fontweight='bold',
            va='top', ha='left', color='#1A237E')

    # Build cell text
    col_labels = ['\n'.join(wrap(str(c), 16)) for c in shown.columns]
    cell_text  = []
    for _, row in shown.iterrows():
        cell_text.append([str(v) for v in row.values])

    tbl = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        loc='center',
        cellLoc='left',
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.auto_set_column_width(col=list(range(n_cols)))

    # Style header row
    for j in range(n_cols):
        tbl[0, j].set_facecolor('#1A237E')
        tbl[0, j].set_text_props(color='white', fontweight='bold')

    # Alternate row shading
    for i in range(1, n_rows + 1):
        bg = '#EEF2FF' if i % 2 == 0 else 'white'
        for j in range(n_cols):
            tbl[i, j].set_facecolor(bg)

    # Notes + truncation warning
    foot_lines = []
    if truncated:
        foot_lines.append(f'[Showing first {MAX_ROWS} of {len(df)} rows]')
    if notes:
        foot_lines.extend(wrap(f'Note: {notes}', width=160))
    if foot_lines:
        ax.text(0.0, -0.02, '\n'.join(foot_lines),
                transform=ax.transAxes, fontsize=7, va='top', ha='left',
                color='#555555', style='italic')

    plt.tight_layout(pad=0.5)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Build PDF
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 60)
print("COMPILING REVIEW PDF")
print("=" * 60)

with PdfPages(str(PDF)) as pdf:

    # ── Cover page ────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(11.69, 8.27))
    ax.set_facecolor('#0D1B2A')
    fig.patch.set_facecolor('#0D1B2A')
    ax.text(0.5, 0.72, 'PAD4-DB v2', transform=ax.transAxes,
            fontsize=34, fontweight='bold', color='white', ha='center')
    ax.text(0.5, 0.60, 'Figure & Table Review Package', transform=ax.transAxes,
            fontsize=18, color='#90CAF9', ha='center')
    ax.text(0.5, 0.50, 'Nature Scientific Data — Submission Draft', transform=ax.transAxes,
            fontsize=13, color='#AAAAAA', ha='center')
    ax.text(0.5, 0.38, '3,093 PAD4 inhibitors  ·  94 severe cliffs  ·  4 hub compounds  ·  1,244 scaffolds',
            transform=ax.transAxes, fontsize=11, color='#80DEEA', ha='center')
    ax.text(0.5, 0.28, 'Generated 2026-06-19',
            transform=ax.transAxes, fontsize=10, color='#777777', ha='center')
    ax.axis('off')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ── Main figures ──────────────────────────────────────────────────────────
    section_page(pdf, 'MAIN FIGURES',
                 'Figures 1–6  ·  Nature Methods / Scientific Data format (600 dpi)',
                 color='#1A237E')

    for fname, caption in MAIN_FIGURES:
        img_path = ROOT / 'publication/figures/main' / fname
        if not img_path.exists():
            print(f"  SKIP (not found): {fname}")
            continue
        print(f"  Adding main figure: {fname}")
        figure_page(pdf, str(img_path), caption)

    # ── Supplementary figures ─────────────────────────────────────────────────
    section_page(pdf, 'SUPPLEMENTARY FIGURES',
                 'Figures S1–S9', color='#004D40')

    for fname, caption in SUPP_FIGURES:
        img_path = ROOT / 'publication/figures/supplementary' / fname
        if not img_path.exists():
            print(f"  SKIP (not found): {fname}")
            continue
        print(f"  Adding supp figure: {fname}")
        figure_page(pdf, str(img_path), caption)

    # ── Tables ────────────────────────────────────────────────────────────────
    section_page(pdf, 'TABLES',
                 'Main tables 1–4  ·  Supplementary tables S1–S6  ·  New supplementary analyses',
                 color='#4A148C')

    for title, csv_path, notes in TABLES:
        if not Path(csv_path).exists():
            print(f"  SKIP (not found): {csv_path}")
            continue
        print(f"  Adding table: {title[:50]}")
        table_page(pdf, title, csv_path, notes)

    # ── Change log ────────────────────────────────────────────────────────────
    section_page(pdf, 'CHANGE LOG',
                 'Summary of all modifications made in this revision session',
                 color='#B71C1C')

    changes = [
        ('FIG 4 v7',     'Complete layout rework: was a squished 12.8x2.9in 1x3 strip → now a '
                         'balanced 2x2 (DOUBLE x 5.4in). New panels: (a) top-15 series coloured by '
                         'mean pIC50 + patent-exclusive outlines, (b) Lorenz curve (Gini=0.532), '
                         '(c) NEW SAR ruggedness panel (series size vs intra-scaffold pIC50 spread) '
                         'directly visualising the "scaffold-dependent SAR ruggedness" claim, '
                         '(d) cliff density (series n>=4).'),
        ('FIG 1 v7',     'Replaced noisy multi-ring KDE line contours with a smooth filled density '
                         'underlay (Blues, alpha 0.22). Tightened axis limits to data envelope to '
                         'remove dead whitespace. Reads as a clean landscape cloud.'),
        ('FIG 6b v7',    'Removed redundant "(n=45)" from y-labels; bar tips now show count + % of 80 '
                         '(e.g. "45 (56%)"). Taller bar panel; tighter molecule grid spacing.'),
        ('FIG 5 v7',     'Panel b non-hub bars now labelled with real InChIKey skeletons (traceable) '
                         'instead of generic "Compound 5-12". X-axis clarified as network degree.'),
        ('FIG S10 NEW',  'Added the canonical Structure-Activity Similarity (SAS) map: (a) hexbin '
                         'density of all 358,416 related pairs with the 94 severe cliffs highlighted; '
                         '(b) |ΔpIC50| distribution by similarity bin showing diagonal absence. '
                         'Companion table supp_sas_quadrants.csv (cliffs = 0.026% of related pairs).'),
        ('VALIDATION',   '00_validate_canonical_numbers.py created — pre-flight guard asserting '
                         '3,093 compounds, 94 severe cliffs, 4 hub compounds, 1,244 scaffolds, '
                         '80 MMP-confirmed pairs, 13 ecfp4-only pairs.'),
        ('STYLE',        'figure_style.py: SEM["cliff"] changed from #CC3311 (red) to #AA0044 '
                         '(dark magenta) to separate cliff edges from Hub B red. Added SEM["background"], '
                         'SEM["gain"], SEM["loss"] keys.'),
        ('DATA',         'pad4_compounds.parquet + shared_assets.parquet: hub_class column added '
                         '("A"/"B"/"none") for the 4 hub compounds.'),
        ('FIG 1',        't-SNE landscape redesigned: layered background scatter + gaussian_kde '
                         'contours + viridis-coloured cliff compounds + hub markers (★ / ◆). '
                         'Hub A star s=100, Hub B diamond s=40 (proportional). '
                         'Legend cliff-proxy now uses viridis midpoint color.'),
        ('FIG 2',        'Source overlap: 3-panel (UpSet / source bars / independence stacked bar). '
                         'Panel b x-axis widened (1.70×) to stop label clipping. '
                         'Redundant legend removed from panel b (y-axis labels sufficient). '
                         'Panel c legend repositioned.'),
        ('FIG 3',        'Potency distribution: 3-panel. Enzymatic BAEE bar changed from '
                         '#BBBBBB (invisible grey) to #0077BB (blue). '
                         'All mechanism bars now distinctively coloured.'),
        ('FIG 4',        'Scaffold landscape: cliff density scatter (panel b) now filters to '
                         'series with n≥4 compounds. Previously n=2 scaffolds produced trivially '
                         'high cliff density = 1.0, creating a misleading cluster at y=1.0. '
                         'Max cliff density: 1.0 → 0.167. Annotation repositioned inside axis.'),
        ('FIG 5',        'Cliff network: node size scale reduced (18–300 → 10–140) to prevent '
                         'hub nodes dominating the plot. Hub A fixed s=180 (star), '
                         'Hub B fixed s=52 (diamond), lw=1.5. Edge alpha raised to 0.60 '
                         'baseline and lw=0.75 so edges visible on white background.'),
        ('FIG 6b',       'CRITICAL BUG FIXED: MMP bar chart previously showed counts from all '
                         '707 MMP pairs in mmp_pairs_cliff99.csv (345/259/101), but labels said '
                         'n=45/27/8. Fixed by joining severe cliff pairs × MMP file on canonical '
                         'pair key. Now correctly shows 45/27/8 = 80 MMP-confirmed severe pairs. '
                         'Labels dynamically generated from actual join counts.'),
        ('TABLE S5',     'tableS5_top20_sali_pairs.csv: 3 duplicate rows removed. '
                         'Dedup used canonical pair key (sorted A/B) not rank number. '
                         'Result: 20 rows → 17 rows.'),
        ('MANUSCRIPT',   'Tables 3+4 merged into single source distribution table with '
                         'Mean pIC50 column. Table 5 standardised (Yes/No/N/A, no checkmarks). '
                         'S6 condensed to "Top 10 shown; full 94-pair dataset in CSV". '
                         'Conclusion updated with locked scaffold-dependent SAR ruggedness statement.'),
        ('NEW ANALYSES', 'Three new supplementary scripts: '
                         '(1) supp_hub_properties.py — hub vs non-hub cliff physicochemical comparison; '
                         'only pIC50 significant (p=0.007). '
                         '(2) supp_assay_cliff_enrichment.py — no mechanism class enriched in cliffs. '
                         '(3) supp_patent_analysis.py — patent depleted in cliffs (OR=0.121, p=0.006, '
                         'phi=0.045 negligible effect); fig_s09 generated.'),
    ]

    fig, ax = plt.subplots(figsize=(11.69, 8.27))
    ax.axis('off')
    fig.patch.set_facecolor('white')

    y = 0.97
    ax.text(0.0, y, 'Change Log — PAD4-DB v2 Manuscript Revision', transform=ax.transAxes,
            fontsize=13, fontweight='bold', color='#B71C1C', va='top')
    y -= 0.06

    for tag, desc in changes:
        ax.text(0.0, y, f'[{tag}]', transform=ax.transAxes,
                fontsize=8.5, fontweight='bold', color='#1A237E', va='top')
        wrapped = wrap(desc, width=130)
        for line in wrapped:
            y -= 0.038
            ax.text(0.01, y, line, transform=ax.transAxes,
                    fontsize=7.5, color='#222222', va='top')
        y -= 0.025
        if y < 0.05:
            break

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # PDF metadata
    d = pdf.infodict()
    d['Title']   = 'PAD4-DB v2 — Figure & Table Review'
    d['Author']  = 'PAD4-DB v2 Pipeline'
    d['Subject'] = 'Nature Scientific Data submission review'
    d['Keywords']= 'PAD4, SAR, activity cliffs, scaffold, MMP'

print(f"\nPDF saved: {PDF}")
print(f"File size: {PDF.stat().st_size / 1024:.0f} KB")
