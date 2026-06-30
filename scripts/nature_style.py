"""
PAD4-DB v2 — Nature Figure Design System
Apply with: from nature_style import apply_nature_style, COLORS, panel_label

All figures must use these settings for journal submission.
"""
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np

# ── Dimensions (in inches) ────────────────────────────────────────────────────
SINGLE_COL_W = 3.504   # 89 mm
DOUBLE_COL_W = 7.205   # 183 mm
MAX_HEIGHT   = 9.724   # 247 mm

# ── Color palette (colorblind-safe) ──────────────────────────────────────────
COLORS = {
    'blue':       '#0077BB',   # primary data / main series
    'orange':     '#EE7733',   # patent-exclusive / highlight
    'red':        '#CC3311',   # severe cliffs / Hub Class B
    'teal':       '#009988',   # enzymatic_confirmed / secondary
    'cyan':       '#33BBEE',   # fp_ic50
    'navy':       '#1A237E',   # Hub Class A
    'gray_light': '#BBBBBB',   # background / non-highlighted
    'gray_dark':  '#555555',   # secondary data
    'black':      '#000000',   # hub stars / emphasis
    # Named aliases for readability
    'severe':              '#CC3311',
    'moderate':            '#EE7733',
    'broad':               '#0077BB',
    'enzymatic':           '#BBBBBB',
    'enzymatic_confirmed': '#009988',
    'fp_ic50':             '#33BBEE',
    'covalent':            '#CC3311',
    'hub_a':               '#1A237E',
    'hub_b':               '#CC3311',
    'patent':              '#EE7733',
    'published':           '#0077BB',
    'multi_source':        '#CC3311',
    'redundant':           '#0077BB',
}

# ── RC param block ─────────────────────────────────────────────────────────────
RC_PARAMS = {
    # Font
    'font.family':          'sans-serif',
    'font.sans-serif':      ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':            7,
    # Axes
    'axes.labelsize':       7,
    'axes.titlesize':       7,
    'axes.linewidth':       0.75,
    'axes.spines.top':      False,
    'axes.spines.right':    False,
    'axes.grid':            False,
    # Ticks
    'xtick.labelsize':      6,
    'ytick.labelsize':      6,
    'xtick.direction':      'in',
    'ytick.direction':      'in',
    'xtick.major.size':     3,
    'ytick.major.size':     3,
    'xtick.minor.size':     2,
    'ytick.minor.size':     2,
    'xtick.major.width':    0.5,
    'ytick.major.width':    0.5,
    # Lines and markers
    'lines.linewidth':      0.75,
    'lines.markersize':     4,
    'patch.linewidth':      0.5,
    # Legend
    'legend.fontsize':      6,
    'legend.frameon':       False,
    'legend.handlelength':  1.5,
    'legend.handletextpad': 0.5,
    # Saving
    'savefig.dpi':          600,
    'savefig.bbox':         'tight',
    'savefig.facecolor':    'white',
    'figure.facecolor':     'white',
    # Layout
    'figure.constrained_layout.use': True,
}


def apply_nature_style():
    """Apply all Nature rcParams. Call once at script top."""
    matplotlib.rcParams.update(RC_PARAMS)


def panel_label(ax, letter, x=-0.12, y=1.04):
    """Add bold panel label (A, B, C...) outside upper-left of axes."""
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=8, fontweight='bold',
            va='bottom', ha='right',
            fontfamily='Arial')


def clean_axes(ax):
    """Remove top/right spines (redundant if rcParams set, but safe to call)."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def save_figure(fig, base_path, dpi=600):
    """Save PNG (600 dpi) + SVG with Nature white background."""
    png_path = base_path if base_path.endswith('.png') else base_path + '.png'
    svg_path = png_path.replace('.png', '.svg')
    fig.savefig(png_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    fig.savefig(svg_path,          bbox_inches='tight', facecolor='white')
    import os
    print(f"Saved: {png_path}  ({os.path.getsize(png_path)/1024:.1f} KB)")
    print(f"Saved: {svg_path}  ({os.path.getsize(svg_path)/1024:.1f} KB)")
    return png_path, svg_path


def mannwhitney_label(a, b, alternative='two-sided'):
    """Return formatted p-value string from Mann-Whitney U test."""
    from scipy.stats import mannwhitneyu
    _, p = mannwhitneyu(a, b, alternative=alternative)
    if p < 0.001:
        return 'p < 0.001'
    elif p < 0.01:
        return f'p = {p:.3f}'
    else:
        return f'p = {p:.2f}'


def mean_sd_band(ax, x, values_list, color, label=None, alpha=0.15):
    """Plot mean ± SD as line + shaded band for a list of value arrays."""
    arr = np.array(values_list)
    mu = arr.mean(axis=0)
    sd = arr.std(axis=0)
    ax.plot(x, mu, color=color, linewidth=0.75, label=label)
    ax.fill_between(x, mu - sd, mu + sd, color=color, alpha=alpha)


# ── Colormap shortcuts ────────────────────────────────────────────────────────
CMAP_ACTIVITY = 'viridis'    # pIC50, activity
CMAP_DENSITY  = 'Blues'      # pair density
CMAP_DIVERGE  = 'RdBu_r'     # diverging (if needed)

# ── Quick legend helpers ──────────────────────────────────────────────────────
def color_patch(color, label):
    return mpatches.Patch(color=color, label=label)

def color_line(color, label, lw=0.75, ls='-'):
    return mlines.Line2D([], [], color=color, linewidth=lw, linestyle=ls, label=label)
