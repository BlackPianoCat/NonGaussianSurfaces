"""Regenerate every figure used in the README and mkdocs site (mkdocs/docs/assets/).

Run from the repo root:  python src/examples/make_readme_figures.py
Credits: Sebastian Korsak, with Claude. Surface model: Yang et al. (2014).
"""
from pathlib import Path

from nongauss import generate_surface, plots, interactive

ROOT = Path(__file__).resolve().parents[2]
ASSETS = ROOT / 'mkdocs' / 'docs' / 'assets'
IMG, HTML = ASSETS / 'img', ASSETS / 'interactive'
BASE = dict(N=256, rms=3, corlength_x=12, corlength_y=12, alpha=0.8, seed=7)


def main():
    IMG.mkdir(parents=True, exist_ok=True); HTML.mkdir(parents=True, exist_ok=True)

    res = generate_surface(skewness=1, kurtosis=4, **BASE)
    plots.save(plots.plot_surfaces(res), IMG / 'surfaces_3d.png')
    plots.save(plots.plot_heightmap_mask(res), IMG / 'heightmap_mask.png')
    plots.save(plots.plot_distribution(res), IMG / 'distribution.png')
    plots.save(plots.plot_hhcf_report(res), IMG / 'hhcf_report.png')
    plots.save(plots.plot_profile_bearing(res), IMG / 'profile_bearing.png')
    plots.save(plots.plot_spectrum(res), IMG / 'spectrum.png')
    fig = plots.plot_feasibility([(1, 4, 'this README (SB)'), (0.5, 4, 'SU'), (-1, 8, 'SU, kurt 8'), (0, 2.2, 'flat (SB)')])
    plots.save(fig, IMG / 'feasibility.png')
    plots.save(plots.plot_logo(generate_surface(skewness=1.2, kurtosis=5, **{**BASE, 'N': 160, 'seed': 3})),
               IMG / 'logo.png', dpi=100)

    # skewness: negative (pits) vs positive (peaks); kurtosis: flat vs spiky; anisotropy
    neg, pos = (generate_surface(skewness=s, kurtosis=4, **BASE) for s in (-1, 1))
    plots.save(plots.plot_compare(neg, pos, ['skewness = -1 (deep pits)', 'skewness = +1 (sharp peaks)']),
               IMG / 'compare_skewness.png')
    flat, spiky = (generate_surface(skewness=0, kurtosis=k, **BASE) for k in (2.2, 8))
    plots.save(plots.plot_compare(flat, spiky, ['kurtosis = 2.2 (flat, bounded)', 'kurtosis = 8 (spiky, heavy tails)']),
               IMG / 'compare_kurtosis.png')
    iso, aniso = (generate_surface(skewness=0.5, kurtosis=3.5, **{**BASE, 'corlength_x': cx}) for cx in (12, 40))
    plots.save(plots.plot_compare(iso, aniso, ['isotropic (xi_x = xi_y = 12)', 'anisotropic (xi_x = 40, xi_y = 12)']),
               IMG / 'compare_anisotropy.png')

    rough, smooth = (generate_surface(skewness=0.5, kurtosis=3.5, **{**BASE, 'alpha': a}) for a in (0.3, 1.0))
    plots.save(plots.plot_compare(rough, smooth, ['alpha = 0.3 (jagged)', 'alpha = 1.0 (smooth)']),
               IMG / 'compare_alpha.png')

    # interactive versions of every plot, plus parameter explorers with a slider
    for kind in plots.KINDS:
        plots.save(plots.plot(res, kind, interactive=True), HTML / f'{kind}.html')
    for p, vals, fixed in (('skewness', (-1.5, -1, 0, 1, 1.5), {'kurtosis': 4}), ('kurtosis', (2.5, 3, 4, 6, 9), {'skewness': 0.5})):
        plots.save(interactive.explorer(p, vals, **fixed), HTML / f'explorer_{p}.html')
    print('figures written to', IMG, 'and', HTML)


if __name__ == '__main__':
    main()
