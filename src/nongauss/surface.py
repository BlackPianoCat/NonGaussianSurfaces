"""Surface generation (Gaussian scaffold + Johnson marginals via rank mapping) and statistics.

Credits
-------
- Method: Yang et al., CMES 103(4), 251-279 (2014).
- Original MATLAB surface generator: Dr. Vasilis Constantoudis.
- Python port, extensions and maintenance: Sebastian Korsak (rewritten with Claude).
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.stats import skew, kurtosis as _kurtosis

from scipy.stats import probplot

from .johnson import f_johnson_M, sample_johnson, johnson_dist


@dataclass
class SurfaceResult:
    """Generated surface and the metadata needed for plots/reports."""
    z: np.ndarray                 # final surface (min height = 0)
    z_gauss: np.ndarray           # Gaussian scaffold (same spatial structure)
    coef: tuple = None            # Johnson (gamma, delta, xi, lam)
    j_type: str = 'Gaussian'      # SN / SL / SB / SU / ST / Gaussian
    msg: str = ''                 # fit warning, if any
    requested: dict = field(default_factory=dict)


def gaussian_surface(N=500, rms=3.0, corlength_x=10, corlength_y=10, alpha=0.9, rng=None):
    """Gaussian surface with autocorrelation rms^2 * exp(-r^(2 alpha)), via FFT filtering."""
    rng = np.random.default_rng(rng)
    t = np.arange(-N / 2, N / 2)                       # unit grid spacing
    TX, TY = np.meshgrid(t, t)                      # x runs along columns (axis 1), y along rows
    r = np.hypot(TX / corlength_x, TY / corlength_y)
    R = np.zeros((N + 1, N + 1))
    R[:N, :N] = rms ** 2 * np.exp(-r ** (2 * alpha))
    amp = np.sqrt(2) * np.abs(np.fft.fft2(R, s=[N, N]))  # sqrt(dtx^2 + dty^2) = sqrt(2)

    noise = rng.random((N, N))
    noise = (noise - noise.mean()) / noise.std()
    z = np.fft.ifft2(np.fft.fft2(noise) * np.sqrt(amp)).real
    return (z - z.mean()) * rms / z.std()


def rank_map(z_spatial, values):
    """Give each pixel of z_spatial the value of equal rank from `values` (keeps spatial layout)."""
    out = np.empty(z_spatial.size)
    out[np.argsort(z_spatial, axis=None)] = np.sort(values, axis=None)
    return out.reshape(z_spatial.shape)


def generate_surface(N=500, rms=3.0, skewness=1.0, kurtosis=3.0, corlength_x=10, corlength_y=10,
                     alpha=0.9, non_Gauss=True, invert=False, clip_percentiles=None, seed=None):
    """Generate a (non-)Gaussian self-affine surface and return a SurfaceResult.

    non_Gauss=False      -> Gaussian surface only.
    invert               -> flip sign before shifting (also flips skewness).
    clip_percentiles     -> optional (low, high) clip of the noise (biases kurtosis down).
    """
    rng = np.random.default_rng(seed)
    z_gs = gaussian_surface(N, rms, corlength_x, corlength_y, alpha, rng)
    req = dict(rms=rms, skewness=skewness, kurtosis=kurtosis,
               corlength_x=corlength_x, corlength_y=corlength_y, alpha=alpha)

    if not non_Gauss:
        return SurfaceResult(z_gs - z_gs.min(), z_gs - z_gs.min(), requested=dict(req, skewness=0, kurtosis=3))

    coef, j_type, msg = f_johnson_M(0, rms, skewness, kurtosis)
    noise = sample_johnson(coef, j_type, (N, N), rng)
    if clip_percentiles is not None:
        noise = np.clip(noise, *np.percentile(noise, clip_percentiles))
    noise = (noise - noise.mean()) / noise.std() * rms

    z = rank_map(z_gs, noise)
    z = (-z if invert else z)
    return SurfaceResult(z - z.min(), z_gs - z_gs.min(), coef, j_type, msg, req)


def SAimage_fft_2(N=500, rms=3, skewness=1, kurtosis=3, corlength_x=10, corlength_y=10, alpha=0.9,
                  threshold=None, non_Gauss=True, corr=True, invert=False, seed=None,
                  show=True, clip_percentiles=None, plot_style='rich'):
    """Backward-compatible wrapper: returns the surface array, optionally plotting it.

    plot_style: 'rich' (summary figures), 'classic' (3D + mask), 'none'.
    `corr` is kept for compatibility (the HHCF figure is part of 'rich').
    """
    from . import plots  # lazy: keeps matplotlib optional for pure generation
    res = generate_surface(N, rms, skewness, kurtosis, corlength_x, corlength_y, alpha,
                           non_Gauss, invert, clip_percentiles, seed)
    if plot_style in ('rich', 'classic'):
        plots.plot_surfaces(res)
        plots.plot_heightmap_mask(res, threshold)
        if plot_style == 'rich' and non_Gauss:
            plots.plot_distribution(res)
            if corr:
                plots.plot_hhcf_report(res)
        if show:
            import matplotlib.pyplot as plt
            plt.show()
    return res.z


# ----------------------------------------------------------------------
# Statistics and reports
# ----------------------------------------------------------------------


def moments(z):
    """(mean, std, skewness, Pearson kurtosis) of an array."""
    v = np.ravel(z)
    return v.mean(), v.std(), skew(v), _kurtosis(v, fisher=False)


def hhcf_1d(z, axis=1):
    """1-D height-height correlation G(r) = sqrt(<(z(x+r)-z(x))^2>) along an axis."""
    z = np.moveaxis(z, axis, 1)
    n = z.shape[1]
    return np.array([np.sqrt(np.mean((z[:, :n - r] - z[:, r:]) ** 2)) for r in range(n // 2)])


def report_rows(res):
    """Rows (name, requested, achieved, % error or None) for a SurfaceResult."""
    req = res.requested
    ach = moments(res.z)
    rows = []
    for name, r, a in zip(('mean', 'rms', 'skewness', 'kurtosis'),
                          (0, req['rms'], req['skewness'], req['kurtosis']), ach):
        err = None if name == 'mean' or abs(r) < 1e-9 else 100 * (a - r) / abs(r)
        rows.append((name, r, a, err))
    return rows


def report_text(res):
    """Plain-text report card (requested vs. achieved)."""
    head = f"{'metric':<10}{'requested':>11}{'achieved':>11}{'error':>9}"
    lines = [head, '-' * len(head)]
    for name, r, a, e in report_rows(res):
        lines.append(f"{name:<10}{r:>11.3f}{a:>11.3f}{(f'{e:+.1f}%' if e is not None else 'n/a'):>9}")
    q = res.requested
    lines += ['', f"Johnson type : {res.j_type}" + (f" ({res.msg})" if res.msg else ''),
              f"Corr. lengths: x={q['corlength_x']}, y={q['corlength_y']}, alpha={q['alpha']}",
              f"Grid         : {res.z.shape[0]} x {res.z.shape[1]}",
              '(mean n/a: surface is shifted so min height = 0)']
    return '\n'.join(lines)


def print_report(res):
    print(report_text(res))


def abbott_curve(z, n=300):
    """Bearing-area (Abbott-Firestone) curve: (% of surface above height h, h)."""
    h = np.sort(z.ravel())
    idx = np.linspace(0, h.size - 1, n).astype(int)
    return 100 * (1 - idx / h.size), h[idx]


def radial_psd(z):
    """Isotropic power spectrum: (wavenumber k in cycles/pixel, radially averaged |FFT|^2)."""
    n = z.shape[0]
    p = np.abs(np.fft.fftshift(np.fft.fft2(z - z.mean()))) ** 2 / z.size
    ky, kx = np.indices(p.shape) - n // 2
    r = np.hypot(kx, ky).astype(int)
    m = r < n // 2
    psd = np.bincount(r[m], p[m]) / np.bincount(r[m])
    return np.arange(1, n // 2) / n, psd[1:]


def pdf_and_qq(res, n=400):
    """Fitted Johnson PDF on a grid and Q-Q data vs. the fit, in the surface's shifted frame.

    Returns dict(x, pdf, osm, osr, slope, intercept, r) or None for the two-point 'ST' case.
    """
    dist = johnson_dist(res.coef, res.j_type) if res.coef is not None else None
    if dist is None:
        return None
    v = res.z.ravel()
    shift = v.mean() - dist.mean()            # surface is shifted so min height = 0
    x = np.linspace(v.min(), v.max(), n)
    (osm, osr), (slope, intercept, r) = probplot(v - shift, dist=dist, fit=True)
    return dict(x=x, pdf=dist.pdf(x - shift), osm=osm, osr=osr, slope=slope, intercept=intercept, r=r)


def log_power_2d(z):
    """log10 of the centered 2D power spectrum (axes in cycles/pixel, range -0.5..0.5)."""
    return np.log10(np.abs(np.fft.fftshift(np.fft.fft2(z - z.mean()))) ** 2 + 1e-12)
