"""
build_manuscript_from_md.py — Render PAD4_DB_manuscript_DRAFT_v7.md to .docx.

The Markdown file is the single source of truth. This renderer parses it directly
(headings, paragraphs with **bold**/*italic*, bullet & numbered lists, pipe tables,
and `> **[Figure N. ...]**` caption blockquotes that trigger image insertion),
so the docx never drifts from the manuscript text again.

Run: conda run -n pad4bench python publication/scripts/build_manuscript_from_md.py
"""
import re
from pathlib import Path
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from PIL import Image

ROOT = Path('/home/nidhal/PAD4-db_V2')
MD   = ROOT / 'publication/manuscript/PAD4_DB_manuscript_DRAFT_v7.md'
OUT  = ROOT / 'publication/manuscript/PAD4_DB_manuscript_DRAFT_v7.docx'
FIGS = {'main': ROOT / 'publication/figures/main', 'supp': ROOT / 'publication/figures/supplementary'}

FIG_MAP = {
    'Figure 1': 'main/fig01_headline', 'Figure 2': 'main/fig02_source_overlap',
    'Figure 3': 'main/fig03_potency', 'Figure 4': 'main/fig04_scaffold',
    'Figure 5': 'main/fig05_cliff_network', 'Figure 6': 'main/fig06b_cliff_pairs',
    'Supplementary Figure S1': 'supplementary/fig_s01_scaffold_cliff_density',
    'Supplementary Figure S2': 'supplementary/fig_s02_assay_enrichment',
    'Supplementary Figure S3': 'supplementary/fig_s03_sas_map',
    'Supplementary Figure S4': 'supplementary/fig_s04_ruggedness_panels',
    'Supplementary Figure S5': 'supplementary/fig_s05_null_models',
}

doc = Document()
s = doc.sections[0]
s.page_width, s.page_height = Cm(21.0), Cm(29.7)
s.left_margin = s.right_margin = s.top_margin = s.bottom_margin = Cm(2.4)
doc.styles['Normal'].font.name = 'Arial'
doc.styles['Normal'].font.size = Pt(11)


def font(run, size=11, bold=False, italic=False, color=None):
    run.font.name = 'Arial'; run.font.size = Pt(size); run.bold = bold; run.italic = italic
    if color:
        run.font.color.rgb = RGBColor(*color)


_TOKEN = re.compile(r'(\*\*.+?\*\*|\*.+?\*|`.+?`)')
def add_inline(p, text, size=11, base_italic=False):
    for seg in _TOKEN.split(text):
        if not seg:
            continue
        if seg.startswith('**') and seg.endswith('**'):
            font(p.add_run(seg[2:-2]), size=size, bold=True, italic=base_italic)
        elif seg.startswith('*') and seg.endswith('*'):
            font(p.add_run(seg[1:-1]), size=size, italic=True)
        elif seg.startswith('`') and seg.endswith('`'):
            r = p.add_run(seg[1:-1]); font(r, size=size - 0.5, italic=base_italic); r.font.name = 'Consolas'
        else:
            font(p.add_run(seg), size=size, italic=base_italic)


def para(text, size=11, align=WD_ALIGN_PARAGRAPH.JUSTIFY, after=6, indent=None, hanging=False):
    p = doc.add_paragraph(); p.alignment = align
    p.paragraph_format.space_after = Pt(after)
    if indent is not None:
        p.paragraph_format.left_indent = Cm(indent)
        if hanging:
            p.paragraph_format.first_line_indent = Cm(-indent)
    add_inline(p, text, size=size)
    return p


def heading(text, level):
    p = doc.add_paragraph(); sizes = {1: 13, 2: 13, 3: 11.5}
    font(p.add_run(text), size=sizes.get(level, 11), bold=True, italic=(level == 3))
    p.paragraph_format.space_before = Pt(12 if level <= 2 else 9)
    p.paragraph_format.space_after = Pt(4)


def _borders(tbl):
    tblPr = tbl._tbl.find(qn('w:tblPr'))
    if tblPr is None:
        tblPr = OxmlElement('w:tblPr'); tbl._tbl.insert(0, tblPr)
    b = OxmlElement('w:tblBorders')
    for side in ('top', 'bottom'):
        e = OxmlElement(f'w:{side}'); e.set(qn('w:val'), 'single'); e.set(qn('w:sz'), '8'); e.set(qn('w:color'), '000000'); b.append(e)
    ih = OxmlElement('w:insideH'); ih.set(qn('w:val'), 'single'); ih.set(qn('w:sz'), '2'); ih.set(qn('w:color'), 'CCCCCC'); b.append(ih)
    tblPr.append(b)


def shade(cell, color):
    tcPr = cell._tc.get_or_add_tcPr(); sh = OxmlElement('w:shd')
    sh.set(qn('w:fill'), color); sh.set(qn('w:val'), 'clear'); tcPr.append(sh)


def add_table(rows):
    ncol = len(rows[0])
    tbl = doc.add_table(rows=len(rows), cols=ncol); tbl.alignment = WD_TABLE_ALIGNMENT.CENTER
    _borders(tbl)
    for i, row in enumerate(rows):
        for j in range(ncol):
            cell = tbl.rows[i].cells[j]; cell.text = ''
            val = row[j] if j < len(row) else ''
            add_inline(cell.paragraphs[0], val, size=8.5)
            for r in cell.paragraphs[0].runs:
                r.font.size = Pt(8.5)
            if i == 0:
                shade(cell, '1A237E')
                for r in cell.paragraphs[0].runs:
                    r.font.bold = True; r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
            elif i % 2 == 0:
                shade(cell, 'EEF2FF')
    doc.add_paragraph().paragraph_format.space_after = Pt(2)


def add_figure(fig_key):
    rel = FIG_MAP.get(fig_key)
    if not rel:
        return
    png = (ROOT / 'publication/figures' / rel).with_suffix('.png')
    if not png.exists():
        para(f'[FIGURE NOT FOUND: {fig_key}]'); return
    w, h = Image.open(png).size
    width_cm = min(16.0, 21.0 * w / h)
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(8)
    p.add_run().add_picture(str(png), width=Cm(width_cm))


# ── Parse ─────────────────────────────────────────────────────────────────────
lines = MD.read_text().split('\n')
i, first_h1 = 0, True
fig_label_re = re.compile(r'\[((?:Supplementary )?Figure S?\d+)[a-z]?\.')
while i < len(lines):
    ln = lines[i].rstrip('\n')
    st = ln.strip()
    if st == '' or st == '---':
        i += 1; continue
    # tables: consecutive pipe rows
    if st.startswith('|'):
        block = []
        while i < len(lines) and lines[i].strip().startswith('|'):
            block.append(lines[i].strip()); i += 1
        rows = []
        for r in block:
            if re.match(r'^\|[\s:|-]+\|$', r):  # separator
                continue
            cells = [c.strip() for c in r.strip('|').split('|')]
            rows.append(cells)
        if rows:
            add_table(rows)
        continue
    # figure caption blockquote
    if st.startswith('>'):
        body = st.lstrip('> ').strip()
        m = fig_label_re.search(body)
        if m:
            add_figure(m.group(1))
        cap = body
        # strip surrounding [ ] of the bold label for display
        cap = cap.replace('**[', '**').replace(']**', '**', 1)
        p = doc.add_paragraph(); p.paragraph_format.space_after = Pt(10); p.paragraph_format.space_before = Pt(3)
        add_inline(p, cap, size=9)
        i += 1; continue
    # headings
    if st.startswith('### '):
        heading(st[4:], 3); i += 1; continue
    if st.startswith('## '):
        heading(st[3:], 2); i += 1; continue
    if st.startswith('# '):
        if first_h1:
            p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            font(p.add_run(st[2:]), size=16, bold=True); first_h1 = False
        else:
            heading(st[2:], 1)
        i += 1; continue
    # numbered list (references)
    if re.match(r'^\d+\.\s', st):
        para(re.sub(r'^\d+\.\s', '', st), size=9, align=WD_ALIGN_PARAGRAPH.LEFT,
             after=3, indent=0.6, hanging=True)
        # re-add the number
        last = doc.paragraphs[-1]
        num = st.split('.', 1)[0] + '. '
        run = last.runs[0]; run.text = num + run.text
        i += 1; continue
    # bullets
    if st.startswith('- ') or st.startswith('* '):
        para('•  ' + st[2:], size=10, align=WD_ALIGN_PARAGRAPH.LEFT, after=3, indent=0.5, hanging=True)
        i += 1; continue
    # plain paragraph
    para(st)
    i += 1

doc.save(str(OUT))
import zipfile
n_img = len([n for n in zipfile.ZipFile(OUT).namelist() if n.startswith('word/media/')])
print(f"Saved {OUT.name}: {OUT.stat().st_size/1024:.0f} KB · {len(doc.tables)} tables · {n_img} images")
