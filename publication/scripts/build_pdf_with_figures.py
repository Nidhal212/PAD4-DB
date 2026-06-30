"""
build_pdf_with_figures.py — Render PAD4_DB_manuscript_FINAL.md to PDF with figures.

Preprocesses the markdown:
  - Converts **Figure N | caption** lines to pandoc ![](image) + italic caption
  - Strips bare --- separators (replaced by figure context)
  - Passes result to pandoc/xelatex

Run: conda run -n pad4bench python publication/scripts/build_pdf_with_figures.py
"""
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path('/home/nidhal/PAD4-db_V2')
MD   = ROOT / 'publication/manuscript/PAD4_DB_manuscript_FINAL.md'
OUT  = ROOT / 'publication/manuscript/PAD4_DB_manuscript_FINAL.pdf'
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
            # full caption including the bold label
            caption  = st

            rel = FIG_MAP.get(fig_key)
            if rel:
                img_path = FIGS / rel
                if img_path.exists():
                    # pandoc image: centre, 95% text width
                    out_lines.append(
                        f'![]({img_path}){{width=95%}}'
                    )
                    out_lines.append('')
                    # caption as italic paragraph (strip ** for readability)
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
        '--pdf-engine=xelatex',
        '-V', 'geometry:margin=2.4cm',
        '-V', 'fontsize=11pt',
        '-V', 'mainfont=DejaVu Serif',
        '-V', 'monofont=DejaVu Sans Mono',
        '--highlight-style=tango',
    ]

    print(f'Running pandoc → {OUT.name} ...')
    result = subprocess.run(pandoc_cmd, capture_output=True, text=True)

    if result.returncode != 0:
        # fallback: try without font specification (uses LaTeX defaults)
        print('Font fallback: retrying without mainfont/monofont...')
        pandoc_cmd2 = [
            'pandoc', tmp_path,
            '-o', str(OUT),
            '--pdf-engine=xelatex',
            '-V', 'geometry:margin=2.4cm',
            '-V', 'fontsize=11pt',
            '--highlight-style=tango',
        ]
        result = subprocess.run(pandoc_cmd2, capture_output=True, text=True)

    if result.returncode != 0:
        print('ERROR:', result.stderr[-2000:])
        raise RuntimeError('pandoc failed')

    if result.stderr:
        # print only non-trivial warnings
        for line in result.stderr.split('\n'):
            if line.strip() and 'Missing character' not in line:
                print('WARN:', line)

    Path(tmp_path).unlink(missing_ok=True)

    size_kb = OUT.stat().st_size / 1024
    print(f'Saved {OUT.name}: {size_kb:.0f} KB')


if __name__ == '__main__':
    main()
