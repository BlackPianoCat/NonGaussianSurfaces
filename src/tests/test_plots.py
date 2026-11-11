import matplotlib.figure
import plotly.graph_objects as go
import pytest

from nongauss import generate_surface, plots, interactive

KINDS = list(plots.KINDS)


@pytest.mark.parametrize('kind', KINDS)
def test_static_plots_have_at_most_two_panels(res, kind):
    fig = plots.plot(res, kind)
    assert isinstance(fig, matplotlib.figure.Figure)
    panels = [a for a in fig.axes if a.get_label() != '<colorbar>']
    assert 1 <= len(panels) <= 2


@pytest.mark.parametrize('kind', KINDS)
def test_interactive_plots(res, kind):
    fig = plots.plot(res, kind, interactive=True)
    assert isinstance(fig, go.Figure) and len(fig.data) >= 1


def test_compare_logo_and_explorer(res):
    other = generate_surface(N=64, skewness=-1, kurtosis=4, seed=0)
    assert plots.plot_compare(res, other, ['a', 'b']) is not None
    assert isinstance(plots.plot_compare(res, other, ['a', 'b'], interactive=True), go.Figure)
    assert plots.plot_logo(res) is not None
    ex = interactive.explorer('skewness', (-1, 0, 1), N=32)
    assert len(ex.frames) == 3


def test_banner(res):
    big = generate_surface(N=128, seed=1)
    assert plots.plot_banner(res, big) is not None


def test_feasibility_accepts_points():
    assert plots.plot_feasibility([(1, 4, 'x'), (0, 2.2)]) is not None


def test_save_png_and_html(res, tmp_path):
    plots.save(plots.plot_distribution(res), tmp_path / 'a.png')
    plots.save(plots.plot_distribution(res, interactive=True), tmp_path / 'a.html')
    assert (tmp_path / 'a.png').stat().st_size > 0 and (tmp_path / 'a.html').stat().st_size > 0


def test_unknown_kind(res):
    with pytest.raises(ValueError):
        plots.plot(res, 'nope')
