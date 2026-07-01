"""
build_docx_with_figures.py — Render PAD4_DB_manuscript_FINAL.md to DOCX with figures.

Mirrors build_pdf_with_figures.py's markdown preprocessing (Figure N | caption
lines -> pandoc image + italic caption) but targets pandoc's native docx
writer instead of xelatex, so no LaTeX engine is involved. Word will show
each figure as an inline image followed by its italic caption paragraph,
matching the numbering/captions already used in the PDF build.

Run: conda run -n pad4bench python publication/scripts/build_docx_with_figures.py
"""
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path('/home/nidhal/PAD4-db_V2')
MD   = ROOT / 'publication/manuscript/PAD4_DB_manuscript_FINAL.md'
OUT  = ROOT / 'publication/manuscript/PAD4_DB_manuscript_FINAL.docx'
FIGS = ROOT / 'publication/figures'

FIG_MAP = {
    'Figure 1':                  'main/fig01_headline.png',
    'Figure 2':                  'main/fig02_source_overlap.png',
    'Figure 3':                  'main/fig03_potency.png',
    'Figure 4':                  'main/fig04_scaffold.png',
    'Figure 5':                  'main/fig05_cliff_network.png',
    'Figure 6':                  'main/fig06b_cliff_pairs.png',
    'Supplementary Figure S1':   'supplementary/fig_s01_scaffold_cliff_density.png',
    'Supplementary Figure S2':   'supplementary/fig_s02_assay_enrichment.png',
    'Supplementary Figure S3':   'supplementary/fig_s03_sas_map.png',
    'Supplementary Figure S4':   'supplementary/fig_s04_ruggedness_panels.png',
    'Supplementary Figure S5':   'supplementary/fig_s05_null_models.png',
}

FIG_LINE_RE = re.compile(r'^\*\*((?:Supplementary )?Figure S?\d+)\s*\|(.+)')


def preprocess(src: str) -> str:
    """Convert figure caption lines to pandoc image + caption syntax."""
    out_lines = []
    lines = src.split('\n')

    for ln in lines:
        st = ln.strip()

        m = FIG_LINE_RE.match(st)
        if m:
            fig_key  = m.group(1)          # e.g. "Figure 1"
            caption  = st                  # full caption including the bold label

            rel = FIG_MAP.get(fig_key)
            if rel:
                img_path = FIGS / rel
                if img_path.exists():
                    # docx: pandoc scales inline images to fit the page automatically
                    # when no explicit width is given and the image exceeds text width;
                    # cap at 6.5in (standard US-letter text width at 1in margins) to be
                    # safe across reference templates.
                    out_lines.append(f'![]({img_path}){{width=6.5in}}')
                    out_lines.append('')
                    cap_clean = re.sub(r'\*\*', '', caption)
                    out_lines.append(f'*{cap_clean}*')
                    out_lines.append('')
                    continue
                else:
                    out_lines.append(f'**[IMAGE NOT FOUND: {fig_key}]**')
                    out_lines.append(caption)
                    out_lines.append('')
                    continue
            else:
                out_lines.append(f'**[FIGURE KEY NOT MAPPED: {fig_key}]**')
                out_lines.append(caption)
                out_lines.append('')
                continue

        out_lines.append(ln)

    return '\n'.join(out_lines)


def main():
    src = MD.read_text(encoding='utf-8')
    processed = preprocess(src)

    with tempfile.NamedTemporaryFile(
        suffix='.md', mode='w', encoding='utf-8', delete=False
    ) as tmp:
        tmp.write(processed)
        tmp_path = tmp.name

    pandoc_cmd = [
        'pandoc', tmp_path,
        '-o', str(OUT),
        '--from=markdown+smart',
    ]

    print(f'Running pandoc → {OUT.name} ...')
    result = subprocess.run(pandoc_cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print('ERROR:', result.stderr[-3000:])
        raise RuntimeError('pandoc failed')

    if result.stderr:
        for line in result.stderr.split('\n'):
            if line.strip():
                print('WARN:', line)

    Path(tmp_path).unlink(missing_ok=True)

    size_kb = OUT.stat().st_size / 1024
    print(f'Saved {OUT.name}: {size_kb:.0f} KB')


if __name__ == '__main__':
    main()
