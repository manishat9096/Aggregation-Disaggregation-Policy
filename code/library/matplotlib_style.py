"""
Plot styling extracted from the microgrid model's plot_functions.py.

Palettes, the contrast machinery behind the shade ramps, and the axis/legend
cosmetics that were repeated inline in every plotting function.

    import matplotlib_style as ms

    fig, ax = plt.subplots(figsize = ms.FIG_SIZE)
    ...
    ms.style_axes(ax, ylabel = 'Loss (MWh)', title = 'Losses by stage')

Colour conventions this encodes:
  - one hue per entity, shade for a sub-stage within it (stage_colours)
  - linestyle or hatch carries the second dimension, never a second hue
"""

import colorsys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

FIG_SIZE = (10, 4)
GRID_KW  = dict(axis = 'y', lw = 0.4, alpha = 0.4)

PALETTE = ['#4c72b0', '#dd8452', '#55a868', '#c44e52', '#8172b3',
           '#937860', '#da8bc3', '#8c8c8c']

SURFACE   = "#020202"                             #chart surface the shades are checked against
COMP_HUES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4',
             '#008300', '#4a3aa7', '#e34948']          #hue per component
HATCHES   = ['', '///', '...', 'xxx', '\\']       #texture per model
SHADES    = ['#555', '#9a9a9a', '#2a78d6', '#eb6834']    #secondary axis traces
STYLES    = ['-', '--', '-.', ':']                #linestyle per model

TITLE_KW  = dict(loc = 'left', fontsize = 11)
LEGEND_KW = dict(frameon = False, fontsize = 9)

# ---------------------------------------------------------------- paper mode

FIG_COL = (3.5, 2.0)                              #IEEE single column
INK, MUTED = "#1a1a1a", "#8a8a8a"

ROLE = {'frechet': "#D55E00", 'comon': "#0072B2",  #semantic, colourblind-safe (Okabe-Ito)
        'indep'  : "#009E73", 'rho'  : "#7B3294"}

PAPER_RC = {
    "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "dejavuserif",
    "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.linewidth": 0.6, "lines.linewidth": 1.2,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": "#e4e4e4", "grid.linewidth": 0.5,
    "legend.frameon": False, "figure.dpi": 150, "savefig.bbox": "tight",
}


def use_paper(**over):
    """Serif, single-column rcParams for figures headed into a manuscript.

    Global, so call it once near the imports rather than per figure.
    """
    plt.rcParams.update({**PAPER_RC, **over})


# ------------------------------------------------------------------ colours

def colours(keys, palette = PALETTE):
    """{key: colour}, cycling the palette in first-seen order."""
    return {k: palette[i % len(palette)] for i, k in enumerate(keys)}


def luminance(rgb):
    c = np.asarray(rgb, dtype = float)
    c = np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    return float(0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2])


def contrast(rgb, surface = SURFACE):
    hi, lo = sorted((luminance(rgb), luminance(mcolors.to_rgb(surface))), reverse = True)
    return (hi + 0.05) / (lo + 0.05)


def top_lightness(hue, floor = 2.05, surface = SURFACE):
    """Lightest step of a hue that still clears `floor`:1 against the surface.

    Bisected rather than tabulated: a saturated yellow goes invisible far earlier
    than a blue, so one shared ceiling would either wash out half the hues or
    waste the range of the other half.
    """
    h, _, sat = colorsys.rgb_to_hls(*mcolors.to_rgb(hue))
    lo, hi = 0.20, 0.95
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        if contrast(colorsys.hls_to_rgb(h, mid, sat), surface) >= floor: lo = mid
        else: hi = mid
    return lo


def stage_colours(cols, hues = COMP_HUES, sep = '.', spread = 0.5):
    """Hue identifies the entity, shade the stage within it.

    Keys are 'entity.stage', so all stages of one entity come out as shades of a
    single hue rather than unrelated colours. Ten stages never need ten hues.
    """
    ents = list(dict.fromkeys(c.split(sep)[0] for c in cols))
    out = {}
    for i, ent in enumerate(ents):
        mine = [c for c in cols if c.split(sep)[0] == ent]
        rgb  = np.array(mcolors.to_rgb(hues[i % len(hues)]))
        for j, c in enumerate(mine):
            f = 0.0 if len(mine) == 1 else spread * j / (len(mine) - 1)
            out[c] = tuple(rgb + (1.0 - rgb) * f)      #lighten toward white
    return out


def shade_ramp(hue, n, floor = 2.05, surface = SURFACE):
    """`n` steps of one hue, from full saturation up to the readable ceiling."""
    h, l0, sat = colorsys.rgb_to_hls(*mcolors.to_rgb(hue))
    top = top_lightness(hue, floor, surface)
    return [colorsys.hls_to_rgb(h, l0 + (top - l0) * (0 if n == 1 else j / (n - 1)), sat)
            for j in range(n)]


# ------------------------------------------------------------------- axes

def style_axes(ax, xlabel = '', ylabel = None, title = None, rotate = None,
               grid = 'y', box = False, fontsize = None):
    """Grid behind the data, top and right spines off, left-aligned title.

    grid = 'y' (house default), 'both' or False. box = True keeps all four
    spines instead, for panels read as a frame rather than an L.
    """
    ax.set_xlabel(xlabel, **({'fontsize': fontsize} if fontsize else {}))
    if ylabel is not None: ax.set_ylabel(ylabel, **({'fontsize': fontsize} if fontsize else {}))
    if title  is not None: ax.set_title(title, **{**TITLE_KW, **({'fontsize': fontsize + 1} if fontsize else {})})
    if grid:
        ax.grid(**{**GRID_KW, 'axis': 'both' if grid == 'both' else 'y'})
        ax.set_axisbelow(True)
    if box:
        ax.spines[:].set_visible(True); ax.spines[:].set_linewidth(0.6)
    else:
        ax.spines[['top', 'right']].set_visible(False)
    if fontsize: ax.tick_params(labelsize = fontsize, width = 0.6)
    if rotate is not None:
        plt.setp(ax.get_xticklabels(), rotation = rotate, ha = 'right')
    return ax


def style_all(axes, **kw):
    """style_axes over a tuple or array of axes."""
    for a in np.ravel(np.asarray(axes, dtype = object)): style_axes(a, **kw)
    return axes


# ------------------------------------------------------------------ legends

def colour_handles(labels, col = None, lw = 2):
    """Line handles keyed by colour, for the primary dimension."""
    col = col or colours(labels)
    return [mlines.Line2D([0], [0], color = col[l], lw = lw, label = l) for l in labels]


def style_handles(labels, color = '#555', lw = 1.6, styles = STYLES):
    """Line handles keyed by linestyle, for the secondary dimension."""
    return [mlines.Line2D([0], [0], color = color, ls = styles[i % len(styles)], lw = lw, label = l)
            for i, l in enumerate(labels)]


def hatch_handles(labels, facecolor = '#dddddd', hatches = HATCHES):
    """Patch handles keyed by hatch, for the secondary dimension on bars."""
    return [mpatches.Patch(facecolor = facecolor, hatch = hatches[i % len(hatches)], label = l)
            for i, l in enumerate(labels)]


def patch_handles(labels, hues = COMP_HUES):
    """Patch handles keyed by hue, for the primary dimension on bars."""
    return [mpatches.Patch(facecolor = hues[i % len(hues)], label = l)
            for i, l in enumerate(labels)]


def legend_above(ax, handles, ncol = None, **kw):
    """Single row of keys above the axes, out of the data."""
    return ax.legend(handles = handles, ncol = ncol or len(handles),
                     bbox_to_anchor = (0.0, 1.02), loc = 'lower left',
                     handlelength = 1.4, columnspacing = 1.4, **{**LEGEND_KW, **kw})


def legend_right(ax, handles, **kw):
    """Keys stacked to the right; pair with fig.subplots_adjust(right = 0.80)."""
    return ax.legend(handles = handles, ncol = 1,
                     bbox_to_anchor = (1.01, 1.0), loc = 'upper left',
                     **{**LEGEND_KW, 'fontsize': 8, **kw})
