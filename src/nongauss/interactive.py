"""Interactive Plotly counterparts of every plot in `plots.py` (same names, white theme).

You normally call them through `plots.plot_xxx(res, interactive=True)`.
Save any figure as standalone HTML with `fig.write_html(path)`.

Credits: Sebastian Korsak, with Claude; surface model after Yang et al. (2014).
"""
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .johnson import johnson_boundaries
from .surface import (generate_surface, hhcf_1d, abbott_curve, radial_psd, report_rows,
                      pdf_and_qq, log_power_2d)

BLUE, ORANGE, GREY = '#1971c2', '#e8590c', '#868e96'
CS = 'Inferno'
SCENE = dict(aspectmode='manual', aspectratio=dict(x=1, y=1, z=0.4),
             xaxis_visible=False, yaxis_visible=False)


def _layout(fig, **kw):
    fig.update_layout(template='plotly_white', margin=dict(l=50, r=30, t=70, b=50), **kw)
    return fig


def _grid(cols=2, specs=None, titles=(), widths=None):
    return make_subplots(rows=1, cols=cols, specs=[specs or [{}] * cols], subplot_titles=titles,
                         column_widths=widths, horizontal_spacing=0.08)


def _image_axes(fig, n=2):
    """Square pixels and hidden ticks for heatmap panels."""
    for i in range(1, n + 1):
        fig.update_xaxes(visible=False, row=1, col=i)
        fig.update_yaxes(visible=False, scaleanchor=f'x{i if i > 1 else ""}', row=1, col=i)


def plot_surfaces(res, step=2):
    """Gaussian scaffold and final surface as rotatable 3D surfaces with floor contours."""
    zmax = max(res.z.max(), res.z_gauss.max())
    fig = _grid(specs=[{'type': 'scene'}] * 2, titles=('Gaussian scaffold', f'Non-Gaussian (Johnson {res.j_type})'))
    for col, z in enumerate((res.z_gauss, res.z), 1):
        fig.add_trace(go.Surface(z=z[::step, ::step], colorscale=CS, cmin=0, cmax=zmax, showscale=col == 2,
                                 colorbar=dict(title='height', len=0.7),
                                 contours_z=dict(show=True, usecolormap=True, project_z=True)), 1, col)
    fig.update_layout(scene=dict(**SCENE, zaxis=dict(range=[0, zmax])), scene2=dict(**SCENE, zaxis=dict(range=[0, zmax])))
    return _layout(fig)


def plot_heightmap_mask(res, threshold=None):
    """Zoomable heightmap and binarized mask (hover shows exact values)."""
    thr = res.z.mean() if threshold is None else threshold
    fig = _grid(titles=('Heightmap', f'Above z = {thr:.2f} ({100 * (res.z > thr).mean():.0f}% of area)'))
    fig.add_trace(go.Heatmap(z=res.z, colorscale=CS, colorbar=dict(title='height', x=0.46, len=0.9)), 1, 1)
    fig.add_trace(go.Heatmap(z=(res.z > thr).astype(int), colorscale='Greys', showscale=False), 1, 2)
    _image_axes(fig)
    return _layout(fig)


def plot_distribution(res):
    """Height histogram with fitted Johnson PDF, and Q-Q plot against the fit."""
    fit = pdf_and_qq(res)
    fig = _grid(titles=('Height distribution', 'Q-Q vs. fit' + (f' (R = {fit["r"]:.4f})' if fit else '')))
    fig.add_trace(go.Histogram(x=res.z.ravel(), nbinsx=100, histnorm='probability density',
                               marker_color=BLUE, opacity=0.45, name='surface'), 1, 1)
    if fit:
        fig.add_trace(go.Scatter(x=fit['x'], y=fit['pdf'], line=dict(color=ORANGE, width=3),
                                 name=f'Johnson {res.j_type} fit'), 1, 1)
        fig.add_trace(go.Scatter(x=fit['osm'], y=fit['osr'], mode='markers', name='quantiles',
                                 marker=dict(color=BLUE, size=4, opacity=0.5)), 1, 2)
        fig.add_trace(go.Scatter(x=fit['osm'], y=fit['slope'] * fit['osm'] + fit['intercept'],
                                 line=dict(color=ORANGE, width=2), showlegend=False), 1, 2)
    fig.update_xaxes(title_text='height', row=1, col=1); fig.update_yaxes(title_text='density', row=1, col=1)
    fig.update_xaxes(title_text='theoretical quantiles', row=1, col=2)
    fig.update_yaxes(title_text='sample quantiles', row=1, col=2)
    return _layout(fig, barmode='overlay')


def plot_profile_bearing(res):
    """Height profile through the middle row, and the Abbott-Firestone bearing curve."""
    row = res.z.shape[0] // 2
    fig = _grid(titles=(f'Height profile along x (row {row})', 'Bearing-area (Abbott-Firestone) curve'))
    for z, c, lab in ((res.z_gauss, GREY, 'Gaussian scaffold'), (res.z, ORANGE, 'non-Gaussian')):
        fig.add_trace(go.Scatter(y=z[row], line=dict(color=c, width=1.5), name=lab, legendgroup=lab), 1, 1)
        pct, h = abbott_curve(z)
        fig.add_trace(go.Scatter(x=pct, y=h, line=dict(color=c, width=2.5), name=lab, legendgroup=lab,
                                 showlegend=False), 1, 2)
    fig.update_xaxes(title_text='x (pixels)', row=1, col=1); fig.update_yaxes(title_text='height', row=1, col=1)
    fig.update_xaxes(title_text='% of area above height', range=[0, 100], row=1, col=2)
    return _layout(fig)


def plot_hhcf_report(res):
    """Height-height correlation along x and y, next to the requested-vs-achieved table."""
    fig = _grid(specs=[{'type': 'xy'}, {'type': 'table'}], titles=('Height-height correlation', ''), widths=[0.58, 0.42])
    for axis, name, c in ((1, 'x', BLUE), (0, 'y', ORANGE)):
        g = hhcf_1d(res.z, axis)
        fig.add_trace(go.Scatter(x=np.arange(1, len(g)), y=g[1:], line=dict(color=c, width=2.5),
                                 name=f'along {name}'), 1, 1)
    rows = report_rows(res)
    fmt = lambda e: 'n/a' if e is None else f'{e:+.1f}%'
    fig.add_trace(go.Table(header=dict(values=['metric', 'requested', 'achieved', 'error'], fill_color='#f1f3f5'),
                           cells=dict(values=[[r[0] for r in rows], [f'{r[1]:.3f}' for r in rows],
                                              [f'{r[2]:.3f}' for r in rows], [fmt(r[3]) for r in rows]],
                                      fill_color='white')), 1, 2)
    fig.update_xaxes(type='log', title_text='r (pixels)', row=1, col=1)
    fig.update_yaxes(type='log', title_text='G(r)', row=1, col=1)
    return _layout(fig)


def plot_spectrum(res):
    """2D power spectrum heatmap and radial power spectra (log-log)."""
    fig = _grid(titles=('2D power spectrum', 'Radial power spectrum'))
    p = log_power_2d(res.z)
    fig.add_trace(go.Heatmap(z=p, x=np.linspace(-0.5, 0.5, p.shape[1]), y=np.linspace(-0.5, 0.5, p.shape[0]),
                             colorscale='Magma', zmin=np.percentile(p, 5), colorbar=dict(title='log10 power', x=0.46, len=0.9)), 1, 1)
    for z, c, lab in ((res.z_gauss, GREY, 'Gaussian scaffold'), (res.z, ORANGE, 'non-Gaussian')):
        k, psd = radial_psd(z)
        fig.add_trace(go.Scatter(x=k, y=psd, line=dict(color=c, width=2.5), name=lab), 1, 2)
    fig.update_xaxes(title_text='kx (cycles/pixel)', row=1, col=1); fig.update_yaxes(title_text='ky (cycles/pixel)', row=1, col=1)
    fig.update_xaxes(type='log', title_text='k (cycles/pixel)', row=1, col=2)
    fig.update_yaxes(type='log', title_text='power', row=1, col=2)
    return _layout(fig)


def plot_feasibility(points, k_max=None):
    """Johnson (skewness, kurtosis) plane with feasible regions; hover on points for values."""
    if hasattr(points, 'requested'):
        q = points.requested
        points = [(q['skewness'], q['kurtosis'], f'request ({points.j_type})')]
    s = np.linspace(0, 3, 300)
    lo, ln = johnson_boundaries(s)
    top = k_max or max(14, 1.2 * max(p[1] for p in points))
    line = lambda y, name, fill, color, dash=None: go.Scatter(
        x=s, y=y, name=name, fill=fill, fillcolor=color, line=dict(color=GREY if fill else ORANGE, width=1.5, dash=dash))
    fig = go.Figure([line(np.ones_like(s), '', None, None).update(showlegend=False, line_width=0),
                     line(lo, 'impossible', 'tonexty', '#e9ecef'),
                     line(ln, 'SB (bounded)', 'tonexty', '#a5d8ff').update(line=dict(color=ORANGE, width=2)),
                     line(np.full_like(s, top), 'SU (unbounded)', 'tonexty', '#ffd8a8').update(line_width=0)])
    fig.add_trace(go.Scatter(x=[0], y=[3], mode='markers', name='Gaussian (SN)', marker=dict(symbol='square', size=10, color='#212529')))
    fig.add_trace(go.Scatter(x=[abs(p[0]) for p in points], y=[p[1] for p in points], mode='markers+text',
                             text=[p[2] if len(p) > 2 else '' for p in points], textposition='top right', name='requests',
                             marker=dict(size=12, color=BLUE, line=dict(color='white', width=1.5))))
    fig.update_xaxes(title_text='|skewness|', range=[0, 3]); fig.update_yaxes(title_text='kurtosis', range=[1, top])
    return _layout(fig, title='Which Johnson family does your request fall in?')


def plot_compare(res_a, res_b, labels):
    """Two heightmaps side by side on a shared color scale."""
    vmax = max(res_a.z.max(), res_b.z.max())
    fig = _grid(titles=labels)
    for col, r in enumerate((res_a, res_b), 1):
        fig.add_trace(go.Heatmap(z=r.z, colorscale=CS, zmin=0, zmax=vmax, showscale=col == 2,
                                 colorbar=dict(title='height')), 1, col)
    _image_axes(fig)
    return _layout(fig)


def explorer(param='skewness', values=(-1.5, -1, 0, 1, 1.5), N=96, **fixed):
    """Slider over one parameter (skewness, kurtosis, rms or alpha): 3D surface + histogram.

    All surfaces share the seed, so only the chosen parameter changes.
    """
    fixed = {'seed': 0, 'corlength_x': 8, 'corlength_y': 8, **fixed}
    results = [generate_surface(N=N, **{**fixed, param: v}) for v in values]
    zmax = max(r.z.max() for r in results)
    traces = lambda r: [go.Surface(z=r.z, colorscale=CS, cmin=0, cmax=zmax, showscale=False),
                        go.Histogram(x=r.z.ravel(), nbinsx=80, marker_color=BLUE)]

    fig = _grid(specs=[{'type': 'scene'}, {'type': 'xy'}], titles=('Surface', 'Height distribution'), widths=[0.62, 0.38])
    for col, t in enumerate(traces(results[0]), 1):
        fig.add_trace(t, 1, col)
    fig.frames = [go.Frame(data=traces(r), traces=[0, 1], name=str(v)) for r, v in zip(results, values)]
    steps = [dict(method='animate', label=str(v), args=[[str(v)], dict(mode='immediate', frame=dict(duration=0, redraw=True))])
             for v in values]
    fig.update_layout(scene=dict(**SCENE, zaxis=dict(range=[0, zmax])), showlegend=False,
                      sliders=[dict(active=0, currentvalue=dict(prefix=f'{param}: '), steps=steps)])
    return _layout(fig)
