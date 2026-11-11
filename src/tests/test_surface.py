import numpy as np
import pytest

from nongauss import (generate_surface, gaussian_surface, rank_map, moments, hhcf_1d, abbott_curve,
                      radial_psd, report_text, SAimage_fft_2)


@pytest.mark.parametrize('skew,kurt', [(0, 3), (1, 3), (-1, 5), (1.5, 6), (0.5, 2.2)])
def test_moments_match_request(skew, kurt):
    res = generate_surface(N=256, rms=2, skewness=skew, kurtosis=kurt, seed=1)
    _, std, s, k = moments(res.z)
    assert std == pytest.approx(2, rel=0.02)
    assert s == pytest.approx(skew, abs=0.1)
    assert k == pytest.approx(kurt, rel=0.1)


def test_seed_and_shift():
    a, b = generate_surface(N=64, seed=3).z, generate_surface(N=64, seed=3).z
    assert a.min() == 0 and np.array_equal(a, b)
    assert not np.array_equal(a, generate_surface(N=64, seed=4).z)


def test_gaussian_only_and_invert():
    g = generate_surface(N=64, non_Gauss=False, seed=0)
    assert g.j_type == 'Gaussian' and g.z.shape == (64, 64)
    pos = generate_surface(N=128, skewness=1, kurtosis=4, seed=0)
    neg = generate_surface(N=128, skewness=1, kurtosis=4, seed=0, invert=True)
    assert moments(pos.z)[2] > 0 > moments(neg.z)[2]


def test_rank_map_keeps_order_and_values():
    z, v = np.random.default_rng(0).normal(size=(8, 8)), np.random.default_rng(1).exponential(size=(8, 8))
    out = rank_map(z, v)
    assert np.allclose(np.sort(out, axis=None), np.sort(v, axis=None))
    assert np.array_equal(np.argsort(out, axis=None), np.argsort(z, axis=None))


def test_correlation_length_anisotropy():
    z = gaussian_surface(128, 1, corlength_x=30, corlength_y=4, alpha=0.9, rng=0)
    assert hhcf_1d(z, axis=1)[5] < hhcf_1d(z, axis=0)[5]   # smoother along the long direction


def test_stats_helpers(res):
    pct, h = abbott_curve(res.z)
    assert pct[0] > pct[-1] and h[0] <= h[-1]
    k, psd = radial_psd(res.z)
    assert len(k) == len(psd) and np.all(psd > 0)
    txt = report_text(res)
    assert 'skewness' in txt and res.j_type in txt


def test_legacy_wrapper_returns_array():
    assert SAimage_fft_2(N=32, seed=0, show=False, plot_style='none').shape == (32, 32)
