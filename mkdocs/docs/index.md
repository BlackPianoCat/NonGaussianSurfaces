![NonGaussian-Surfaces](assets/img/logo_banner.png){ style="width:100%" }

Python library that generates rough 2D surfaces with a prescribed **height distribution** (rms, skewness, kurtosis) and **spatial structure** (correlation lengths, roughness exponent), plus tools to inspect them with static (matplotlib) or interactive (Plotly) plots.

## Requirements

- Python 3.10 or newer
- Required: `numpy`, `scipy`, `matplotlib`
- Optional: `plotly` (interactive plots)
- Docs and tests: `mkdocs`, `mkdocs-material`, `mkdocstrings[python]`, `pytest`

## Installation

### From PyPI

```bash
pip install nongaussian-surfaces                  # core + static plots
pip install "nongaussian-surfaces[interactive]"   # + Plotly
```

Install name: `nongaussian-surfaces`; import name: `nongauss`.

### From source

```bash
git clone https://github.com/BlackPianoCat/NonGaussianSurfaces
cd NonGaussianSurfaces
pip install -e ".[dev]"        # editable install: package, plotly, pytest, mkdocs
```

Extras defined in `pyproject.toml`:

| Extra | Adds |
|---|---|
| `interactive` | `plotly` |
| `docs` | `mkdocs`, `mkdocs-material`, `mkdocstrings[python]` |
| `dev` | everything above plus `pytest`, `build`, `twine` |

### Virtual environment (recommended)

```bash
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install "nongaussian-surfaces[interactive]"
```

### Verify

```bash
python -c "import nongauss; print(nongauss.__all__[:3])"
nongauss --help                          # command-line tool installed
python src/run_tests.py                  # from a source checkout: runs the test suite
```

In headless environments (servers, CI) use a non-interactive matplotlib backend: `export MPLBACKEND=Agg`.

## First run

```python
from nongauss import generate_surface, print_report, plots

res = generate_surface(N=256, rms=3, skewness=1, kurtosis=4,
                       corlength_x=12, corlength_y=12, alpha=0.8, seed=7)
print_report(res)                                    # requested vs. achieved statistics
fig = plots.plot_heightmap_mask(res)                 # matplotlib figure
plots.save(fig, "surface.png")
```

or from the shell:

```bash
nongauss --skew 1 --kurt 4 --alpha 0.8 --out out/
```

Next: [Guide](usage.md) for every function and option, [API reference](api.md) for signatures.

## Credits

Algorithm by Dr. Vasilis Constantoudis after Yang et al., *CMES* 103(4), 251-279 (2014). Johnson fitting ported from the MATLAB *Johnson Curve Toolbox* by Dave (2021), Python translation by Max Pierini ([EpiData.it](https://epidata.it)). Package, extensions and maintenance by Sebastian Korsak, with Claude. Apache 2.0 license.
