import numpy as np
import pytest

from nongauss import f_johnson_M, johnson_dist, sample_johnson
from nongauss.johnson import johnson_boundaries


@pytest.mark.parametrize('skew,kurt,expected', [
    (0, 3, 'SN'), (1, 3, 'SB'), (0, 2.2, 'SB'), (0.5, 4, 'SU'), (0, 8, 'SU'), (2, 4, 'ST')])
def test_family_selection(skew, kurt, expected):
    assert f_johnson_M(0, 1, skew, kurt)[1] == expected


@pytest.mark.parametrize('skew,kurt', [(1, 3), (-1, 5), (0.5, 4), (0, 2.5)])
def test_sampling_matches_moments(skew, kurt):
    from scipy.stats import skew as sk, kurtosis as ku
    coef, t, _ = f_johnson_M(0, 2, skew, kurt)
    x = sample_johnson(coef, t, 400_000, np.random.default_rng(0))
    assert x.std() == pytest.approx(2, rel=0.03)
    assert sk(x) == pytest.approx(skew, abs=0.1)
    assert ku(x, fisher=False) == pytest.approx(kurt, rel=0.1)


def test_dist_and_errors():
    coef, t, _ = f_johnson_M(0, 1, 0.5, 4)
    assert johnson_dist(coef, t).pdf(0.0) > 0
    assert johnson_dist(coef, 'ST') is None
    with pytest.raises(ValueError):
        johnson_dist(coef, 'XX')
    with pytest.raises(Exception):
        f_johnson_M(0, -1, 0, 3)


def test_boundaries_order():
    lo, ln = johnson_boundaries(np.array([0, 1, 2]))
    assert np.all(ln > lo) and lo[0] == 1 and ln[0] == pytest.approx(3)
