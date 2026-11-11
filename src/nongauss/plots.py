"""Static matplotlib figures: light theme, hillshaded, at most 2 panels each.

Every `plot_*` function takes a SurfaceResult (or a pair) and returns a figure. Pass
`interactive=True` to get the Plotly version instead of matplotlib (see interactive.py);
use `save(fig, path)` to write PNG (static) or HTML (interactive).

Credits: figure design by Sebastian Korsak, with Claude; surface model after Yang et al. (2014).
"""
from functools import wraps

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LightSource
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
from .johnson import johnson_boundaries
from .surface import hhcf_1d, abbott_curve, radial_psd, report_text, pdf_and_qq, log_power_2d

CMAP = 'inferno'
C1, C2, C3, INK, GRID = '#1971c2', '#e8590c', '#868e96', '#212529', '#e9ecef'   # blue, orange, grey
STYLE = {'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
         'axes.edgecolor': '#adb5bd', 'axes.labelcolor': INK, 'text.color': INK,
         'xtick.color': INK, 'ytick.color': INK, 'axes.spines.top': False, 'axes.spines.right': False,
         'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.8, 'axes.titleweight': 'bold',
         'axes.titlesize': 11.5, 'font.size': 10, 'legend.frameon': False}


def styled_static(fn):
    """Light style only (static-only plots)."""
    @wraps(fn)
    def wrapper(*a, **k):
        with plt.rc_context(STYLE):
            return fn(*a, **k)
    return wrapper


def switchable(fn):
    """Light style for matplotlib; `interactive=True` dispatches to the same-named Plotly function."""
    @wraps(fn)
    def wrapper(*a, interactive=False, **k):
        if interactive:
            from . import interactive as ia      # lazy: plotly is optional
            return getattr(ia, fn.__name__)(*a, **k)
        with plt.rc_context(STYLE):
            return fn(*a, **k)
    return wrapper


def _row(n=2, size=(11, 4.6), **kw):
    return plt.subplots(1, n, figsize=size, constrained_layout=True, **kw)


def _clean(ax):
    """Image-style axis: no ticks, spines or grid."""
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)


def _shade(z, cmap=CMAP, vmin=None, vmax=None, exag=2.0):
    """RGB hillshaded rendering of a heightmap (soft light from the north-west)."""
    ls = LightSource(azdeg=315, altdeg=40)
    return ls.shade(z, cmap=plt.get_cmap(cmap), blend_mode='soft', vert_exag=exag,
                    vmin=z.min() if vmin is None else vmin, vmax=z.max() if vmax is None else vmax)


def _surf3d(ax, z, title, zmax, cmap=CMAP, step=2):
    """Hillshaded 3D surface with a contour 'shadow' on the floor."""
    zs = z[::step, ::step]
    X, Y = np.meshgrid(np.arange(zs.shape[1]) * step, np.arange(zs.shape[0]) * step)
    ax.plot_surface(X, Y, zs, facecolors=_shade(zs, cmap, 0, zmax, 1.5), shade=False,
                    linewidth=0, antialiased=False, rstride=1, cstride=1)
    ax.contourf(X, Y, zs, zdir='z', offset=-0.25 * zmax, levels=14, cmap=cmap, alpha=0.9, vmin=0, vmax=zmax)
    ax.set_zlim(-0.25 * zmax, zmax); ax.view_init(32, -58); ax.set_box_aspect((1, 1, 0.45))
    ax.set_title(title, pad=2)
    ax.set_axis_off()


@switchable
def plot_surfaces(res, cmap=CMAP):
    """Gaussian scaffold vs. final surface in 3D (shared height scale)."""
    fig = plt.figure(figsize=(11.5, 5), constrained_layout=True)
    zmax = max(res.z.max(), res.z_gauss.max())
    for i, (z, t) in enumerate([(res.z_gauss, 'Gaussian scaffold'),
                                (res.z, f'Non-Gaussian (Johnson {res.j_type})')]):
        _surf3d(fig.add_subplot(1, 2, i + 1, projection='3d'), z, t, zmax, cmap)
    return fig


@switchable
def plot_heightmap_mask(res, threshold=None, cmap=CMAP):
    """Hillshaded heightmap with colorbar, and the binarized mask with its area fraction."""
    thr = res.z.mean() if threshold is None else threshold
    fig, (a, b) = _row()
    a.imshow(_shade(res.z, cmap), origin='lower')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(res.z.min(), res.z.max())),
                 ax=a, fraction=0.046, pad=0.02, label='height')
    b.imshow(res.z > thr, cmap='Greys', origin='lower')
    a.set_title('Heightmap (hillshaded)')
    b.set_title(f'Above z = {thr:.2f}  ({100 * (res.z > thr).mean():.0f}% of area)')
    _clean(a); _clean(b)
    return fig


@switchable
def plot_distribution(res):
    """Height histogram with fitted Johnson PDF and mean marker, and Q-Q plot against the fit."""
    v = res.z.ravel()
    fit = pdf_and_qq(res)
    fig, (a, b) = _row()
    a.hist(v, bins=100, density=True, color=C1, alpha=0.35, edgecolor='white', lw=0.3, label='surface')
    a.axvline(v.mean(), color=C3, ls='--', lw=1.2, label=f'mean {v.mean():.2f}')
    if fit:
        a.plot(fit['x'], fit['pdf'], color=C2, lw=2.4, label=f'Johnson {res.j_type} fit')
        b.scatter(fit['osm'], fit['osr'], s=5, color=C1, alpha=0.4, edgecolor='none')
        b.plot(fit['osm'], fit['slope'] * fit['osm'] + fit['intercept'], color=C2, lw=2)
        b.set_title(f"Q-Q vs. fit (R = {fit['r']:.4f})")
        b.set_xlabel('theoretical quantiles'); b.set_ylabel('sample quantiles')
    else:
        b.text(0.5, 0.5, 'Q-Q n/a for the\ntwo-point (ST) case', ha='center', va='center'); b.set_axis_off()
    a.set(title='Height distribution', xlabel='height', ylabel='density'); a.legend()
    return fig


@switchable
def plot_profile_bearing(res):
    """Line profile through the middle of the surface, and the Abbott-Firestone bearing curve."""
    fig, (a, b) = _row(size=(11.5, 4.4))
    row = res.z.shape[0] // 2
    a.plot(res.z_gauss[row], color=C3, lw=1.1, label='Gaussian scaffold')
    a.plot(res.z[row], color=C2, lw=1.4, label='non-Gaussian')
    a.set(title=f'Height profile along x (row {row})', xlabel='x (pixels)', ylabel='height'); a.legend()
    for z, c, lab in ((res.z_gauss, C3, 'Gaussian scaffold'), (res.z, C2, 'non-Gaussian')):
        pct, h = abbott_curve(z)
        b.plot(pct, h, color=c, lw=2, label=lab)
    b.set(title='Bearing-area (Abbott-Firestone) curve', xlabel='% of area above height', ylabel='height')
    b.set_xlim(0, 100); b.legend()
    return fig


@switchable
def plot_hhcf_report(res):
    """Height-height correlation along x and y, next to the requested-vs-achieved report."""
    fig, (a, b) = _row(size=(11.5, 4.4))
    for axis, name, c in ((1, 'x', C1), (0, 'y', C2)):
        g = hhcf_1d(res.z, axis)
        a.loglog(np.arange(1, len(g)), g[1:], color=c, lw=2, label=f'along {name}')
    a.legend(); a.grid(True, which='both', alpha=0.4)
    a.set(title='Height-height correlation', xlabel='r (pixels)', ylabel='G(r)')
    b.axis('off')
    b.text(0, 1, 'Statistical report', fontsize=12, fontweight='bold', va='top')
    b.text(0, 0.9, report_text(res), family='monospace', fontsize=9.5, va='top')
    return fig


@switchable
def plot_spectrum(res):
    """2D power spectrum of the surface, and the radial spectrum with the expected high-k power law."""
    fig, (a, b) = _row(size=(11.5, 4.6))
    p = log_power_2d(res.z)
    im = a.imshow(p, origin='lower', cmap='magma', extent=[-0.5, 0.5, -0.5, 0.5], vmin=np.percentile(p, 5))
    fig.colorbar(im, ax=a, fraction=0.046, pad=0.02, label='log10 power')
    a.set(title='2D power spectrum', xlabel='kx (cycles/pixel)', ylabel='ky (cycles/pixel)'); a.grid(False)
    for z, c, lab in ((res.z_gauss, C3, 'Gaussian scaffold'), (res.z, C2, 'non-Gaussian')):
        k, psd = radial_psd(z)
        b.loglog(k, psd, color=c, lw=1.8, label=lab)
    alpha = res.requested.get('alpha')
    if alpha:                                      # exp(-r^(2a)) has PSD ~ k^-(2+2a) at high k
        k, psd = radial_psd(res.z_gauss)
        i0 = len(k) // 2
        ks = k[len(k) // 3:]
        b.loglog(ks, psd[i0] * (ks / k[i0]) ** (-(2 + 2 * alpha)), 'k:', lw=1.2,
                 label=f'slope -(2+2α) = {-(2 + 2 * alpha):.1f}')
    b.set(title='Radial power spectrum', xlabel='k (cycles/pixel)', ylabel='power'); b.legend()
    return fig


@switchable
def plot_feasibility(points, k_max=None):
    """Johnson (skewness, kurtosis) plane: feasible regions with your request(s) marked.

    `points` is a SurfaceResult or a list of (skewness, kurtosis[, label]).
    """
    if hasattr(points, 'requested'):
        q = points.requested
        points = [(q['skewness'], q['kurtosis'], f'request ({points.j_type})')]
    s = np.linspace(0, 3, 400)
    k_lo, k_ln = johnson_boundaries(s)                       # impossible below k_lo; SL line k_ln
    top = k_max or max(14, 1.2 * max(p[1] for p in points))
    fig, ax = plt.subplots(figsize=(7.5, 5.2), constrained_layout=True)
    ax.fill_between(s, 0, k_lo, color='#f1f3f5', hatch='//', edgecolor='#ced4da', label='impossible')
    ax.fill_between(s, k_lo, k_ln, color='#a5d8ff', label='SB (bounded)')
    ax.fill_between(s, k_ln, top, color='#ffd8a8', label='SU (unbounded)')
    ax.plot(s, k_lo, color=C3, lw=1.4); ax.plot(s, k_ln, color=C2, lw=1.8, label='SL (lognormal) line')
    ax.scatter([0], [3], marker='s', s=60, color=INK, zorder=5, label='Gaussian (SN)')
    for p in points:
        ax.scatter(abs(p[0]), p[1], s=90, color=C1, edgecolor='white', zorder=6)
        ax.annotate(p[2] if len(p) > 2 else f'({p[0]:g}, {p[1]:g})', (abs(p[0]), p[1]),
                    xytext=(8, 6), textcoords='offset points', fontsize=9)
    ax.set(xlim=(0, 3), ylim=(1, top), xlabel='|skewness|', ylabel='kurtosis',
           title='Which Johnson family does your request fall in?')
    ax.legend(loc='upper left')
    return fig


@switchable
def plot_compare(res_a, res_b, labels, cmap=CMAP):
    """Two hillshaded heightmaps side by side on a shared color scale."""
    vmax = max(res_a.z.max(), res_b.z.max())
    fig, axs = _row(size=(11.5, 4.8))
    for ax, r, lab in zip(axs, (res_a, res_b), labels):
        ax.imshow(_shade(r.z, cmap, 0, vmax), origin='lower')
        ax.set_title(lab); _clean(ax)
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, vmax)), ax=axs,
                 fraction=0.025, pad=0.02, label='height')
    return fig


@styled_static
def plot_logo(res, size=3):
    """Small square logo: hillshaded surface, no axes."""
    fig = plt.figure(figsize=(size, size), facecolor='none')
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(_shade(res.z, 'inferno', exag=3), origin='lower'); ax.axis('off')
    return fig


@styled_static
def plot_banner(icon, texture, text='NonGaussian-Surfaces', height=2.2):
    """Logo banner: small surface icon + monospace wordmark filled with a non-Gaussian heightmap.

    `icon` and `texture` are SurfaceResults (texture should be large, e.g. N=512).
    """
    fp = FontProperties(family='DejaVu Sans Mono', weight='bold')
    tp = TextPath((0, 0), text, size=1, prop=fp)
    x0, y0, x1, y1 = tp.get_extents().extents
    w = (x1 - x0) / (y1 - y0) * height * 0.62            # text width for a 0.62*height tall word
    fig = plt.figure(figsize=(height + w + 0.5, height), facecolor='none')
    ax_i = fig.add_axes([0.02 * height / fig.get_figwidth(), 0.06, height * 0.88 / fig.get_figwidth(), 0.88])
    ax_i.imshow(_shade(icon.z, 'inferno', exag=3), origin='lower'); ax_i.axis('off')

    ax = fig.add_axes([(height + 0.12) / fig.get_figwidth(), 0.2, w / fig.get_figwidth(), 0.6])
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.axis('off')
    patch = PathPatch(tp, facecolor='none', edgecolor='none'); ax.add_patch(patch)
    z = texture.z                                         # crop the texture to the word's aspect ratio
    nx = z.shape[1]; ny = max(8, int(nx * (y1 - y0) / (x1 - x0)))
    img = ax.imshow(_shade(z[:ny], 'inferno', z.min(), np.percentile(z, 99.5), 3) ** 0.8,
                    extent=[x0, x1, y0, y1], origin='lower', aspect='auto', interpolation='bicubic', zorder=2)
    img.set_clip_path(patch)
    ax.add_patch(PathPatch(tp, facecolor='none', edgecolor='#2b0a3d', lw=0.5, zorder=3))
    return fig


KINDS = {'surfaces': plot_surfaces, 'heightmap': plot_heightmap_mask, 'distribution': plot_distribution,
         'profile': plot_profile_bearing, 'hhcf': plot_hhcf_report, 'spectrum': plot_spectrum,
         'feasibility': plot_feasibility}


def plot(res, kind='surfaces', interactive=False, **kw):
    """Plot by name: kind in KINDS ('surfaces', 'heightmap', 'distribution', 'profile', 'hhcf', 'spectrum', 'feasibility')."""
    if kind not in KINDS:
        raise ValueError(f'kind must be one of {sorted(KINDS)}, got {kind!r}')
    return KINDS[kind](res, interactive=interactive, **kw)


def save(fig, path, dpi=140, transparent=False):
    """Save a matplotlib figure as an image (closes it; `transparent=True` for no background) or a Plotly figure as HTML."""
    if hasattr(fig, 'write_html'):
        fig.write_html(path, include_plotlyjs='cdn')
        return
    fig.savefig(path, dpi=dpi, facecolor='none' if transparent else 'white', transparent=transparent, bbox_inches='tight')
    plt.close(fig)
