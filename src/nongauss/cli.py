"""Command-line interface: generate a surface, show/save figures and an interactive HTML.

Example:
    nongauss --skew 1 --kurt 4 --alpha 0.7 --out out/ [--interactive]
Credits: Sebastian Korsak, with Claude. Surface model: Yang et al. (2014).
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from . import plots
from .surface import generate_surface, print_report


def parse():
    p = argparse.ArgumentParser(prog='nongauss', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name, default in dict(N=256, rms=3.0, skew=1.0, kurt=3.0, clx=10.0, cly=10.0, alpha=0.9).items():
        p.add_argument(f'--{name}', type=type(default), default=default)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--out', type=Path, default=None, help='save figures here instead of showing them')
    p.add_argument('--interactive', action='store_true', help='interactive Plotly figures (HTML) instead of PNG; needs plotly')
    return p.parse_args()


def main():
    a = parse()
    res = generate_surface(N=a.N, rms=a.rms, skewness=a.skew, kurtosis=a.kurt,
                           corlength_x=a.clx, corlength_y=a.cly, alpha=a.alpha, seed=a.seed)
    print_report(res)
    ia = a.interactive
    figs = {k: plots.plot(res, k, interactive=ia) for k in plots.KINDS if k != 'feasibility'}

    if a.out is None:
        if ia:
            for f in figs.values():
                f.show()
        else:
            plt.show()
        return
    a.out.mkdir(parents=True, exist_ok=True)
    for name, fig in figs.items():
        plots.save(fig, a.out / f"{name}.{'html' if ia else 'png'}")
    print(f'saved to {a.out}')


if __name__ == '__main__':
    main()
