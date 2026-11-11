import matplotlib
import pytest

matplotlib.use('Agg')   # headless plotting


@pytest.fixture(scope='module')
def res():
    """Small non-Gaussian surface shared by the plotting tests."""
    from nongauss import generate_surface
    return generate_surface(N=64, rms=2, skewness=1, kurtosis=4, seed=0)
