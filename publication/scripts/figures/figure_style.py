"""
figure_style.py — shared style for all PAD4-DB v2 figures.
Import and call set_style() at the top of every figure script.

Nature Methods / Scientific Data submission standards.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

# Nature colorblind-safe palette (CLAUDE.md locked)
C = {
    'black':      '#000000',
    'blue':       '#0077BB',   # primary data / main series
    'orange':     '#EE7733',   # patent-exclusive / highlight
    'red':        '#CC3311',   # severe cliffs / Hub Class B
    'teal':       '#009988',   # enzymatic_confirmed
    'cyan':       '#33BBEE',   # fp_ic50
    'navy':       '#1A237E',   # Hub Class A
    'grey':       '#BBBBBB',   # background / non-highlighted
    'gray_dark':  '#555555',   # secondary data
    'lightgrey':  '#DDDDDD',   # borders / boxes
    # Legacy aliases kept for backward compatibility
    'sky':        '#33BBEE',
    'green':      '#009988',
    'vermillion': '#CC3311',
    'purple':     '#555555',
    'yellow':     '#EE7733',
}

# Semantic colour map — same meaning in every figure (CLAUDE.md locked)
SEM = {
    'sar':                  C['blue'],
    'hts':                  C['grey'],
    'active':               C['blue'],
    'inactive':             C['red'],
    'classA':               C['navy'],    # Hub Class A — navy (locked 2026-06-16)
    'classB':               C['red'],     # Hub Class B — red  (locked 2026-06-16)
    'cliff':                '#AA0044',    # Cliff edges / cliff highlight — dark magenta (v6)
    'noncliff':             C['cyan'],
    'single_atom':          C['blue'],
    'small_subst':          C['orange'],
    'medium':               C['teal'],
    'patent':               C['orange'],
    'published':            C['blue'],
    'enzymatic':            C['grey'],
    'enzymatic_confirmed':  C['teal'],
    'fp_ic50':              C['cyan'],
    'covalent':             C['red'],
    'multi_source':         C['red'],
    'redundant':            C['blue'],
    'background':           '#AAAAAA',    # Background scatter (alpha=0.4)
    'gain':                 '#00CC33',    # Gain-of-function atoms in MMP panels
    'loss':                 '#CC0000',    # Loss-of-function atoms in MMP panels
}

# Nature column widths in inches (89 mm / 183 mm)
SINGLE   = 3.504
ONEHALF  = 5.354   # ~136 mm
DOUBLE   = 7.205

def set_style():
    """Apply Nature rcParams. Call once at top of every figure script."""
    mpl.rcParams.update({
        'font.family':          'sans-serif',
        'font.sans-serif':      ['Arial', 'Helvetica', 'Liberation Sans',
                                 'Arimo', 'DejaVu Sans'],
        'font.size':            7,
        'axes.titlesize':       7,
        'axes.labelsize':       7,
        'xtick.labelsize':      6,
        'ytick.labelsize':      6,
        'legend.fontsize':      6,
        'legend.handlelength':  1.5,
        'legend.handletextpad': 0.5,
        'figure.titlesize':     8,
        'axes.linewidth':       0.75,
        'xtick.major.width':    0.5,
        'ytick.major.width':    0.5,
        'xtick.minor.width':    0.5,
        'ytick.minor.width':    0.5,
        'xtick.major.size':     3,
        'ytick.major.size':     3,
        'xtick.minor.size':     2,
        'ytick.minor.size':     2,
        'xtick.direction':      'in',
        'ytick.direction':      'in',
        'lines.linewidth':      0.75,
        'lines.markersize':     4,
        'patch.linewidth':      0.5,
        'axes.spines.top':      False,
        'axes.spines.right':    False,
        'axes.grid':            False,
        'legend.frameon':       False,
        'figure.dpi':           150,
        'savefig.dpi':          600,
        'savefig.bbox':         'tight',
        'savefig.pad_inches':   0.02,
        'savefig.facecolor':    'white',
        'figure.facecolor':     'white',
        'pdf.fonttype':         42,
        'ps.fonttype':          42,
        'svg.fonttype':         'none',
    })
    import matplotlib.font_manager as fm
    resolved = fm.findfont(mpl.rcParams['font.sans-serif'][0], fallback_to_default=False)
    family = fm.FontProperties(fname=resolved).get_name() if resolved else 'FALLBACK'
    print(f"  [style] resolved font: {family} ({resolved})")
    return family

def panel_label(ax, s, x=-0.12, y=1.04):
    """Bold panel letter at consistent position (axes coordinates)."""
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8,
            fontweight='bold', va='bottom', ha='right')

def save_fig(fig, stem):
    """Save PDF (vector, editable fonts) and PNG @ 600 dpi."""
    fig.savefig(f'{stem}.pdf')
    fig.savefig(f'{stem}.png', dpi=600)
    print(f"  Saved: {stem}.pdf + .png")
