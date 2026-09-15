"""The discrete Fréchet distance, checked against an independent brute force."""
import numpy as np
import pytest

from specleash.frechet import FrechetResult, discrete_frechet, frechet_coupling


def all_couplings(n, m):
    """Every monotone path from (0, 0) to (n-1, m-1) with steps (1,0), (0,1), (1,1)."""
    def walk(i, j):
        if (i, j) == (n - 1, m - 1):
            yield [(i, j)]
            return
        for di, dj in ((1, 0), (0, 1), (1, 1)):
            a, b = i + di, j + dj
            if a < n and b < m:
                for rest in walk(a, b):
                    yield [(i, j)] + rest
    yield from walk(0, 0)


def brute_force(p, q):
    """Minimum over every coupling of its longest link. Exponential; tiny curves only."""
    d = np.hypot(p[:, None, 0] - q[None, :, 0], p[:, None, 1] - q[None, :, 1])
    return min(max(d[i, j] for i, j in path) for path in all_couplings(len(p), len(q)))


@pytest.mark.parametrize("n,m", [(1, 1), (1, 4), (4, 1), (2, 3), (3, 3), (4, 5), (5, 5)])
def test_matches_brute_force_on_small_random_curves(n, m):
    rng = np.random.default_rng(n * 10 + m)
    for _ in range(20):
        p, q = rng.normal(size=(n, 2)), rng.normal(size=(m, 2))
        assert discrete_frechet(p, q) == pytest.approx(brute_force(p, q), abs=1e-12)
        assert frechet_coupling(p, q).distance == pytest.approx(brute_force(p, q), abs=1e-12)


def test_symmetric():
    rng = np.random.default_rng(1)
    p, q = rng.normal(size=(40, 2)), rng.normal(size=(55, 2))
    assert discrete_frechet(p, q) == pytest.approx(discrete_frechet(q, p), abs=1e-12)


def test_identical_curves_are_zero():
    p = np.c_[np.linspace(0, 1, 30), np.sin(np.linspace(0, 3, 30))]
    assert discrete_frechet(p, p) == 0.0


def test_a_constant_vertical_offset_is_the_offset():
    x = np.linspace(0, 1, 50)
    assert discrete_frechet(np.c_[x, x], np.c_[x, x + 0.3]) == pytest.approx(0.3, abs=1e-12)


def test_one_displaced_point_sets_the_distance():
    x = np.linspace(0, 1, 50)
    y = np.zeros(50)
    y[25] = 0.9
    assert discrete_frechet(np.c_[x, np.zeros(50)], np.c_[x, y]) == pytest.approx(0.9, abs=1e-12)


def test_the_separation_of_the_first_points_always_counts():
    x = np.linspace(0, 1, 10)
    a = np.where(np.arange(10) == 0, 1.0, 0.0)
    assert discrete_frechet(np.c_[x, a], np.c_[x, np.zeros(10)]) == pytest.approx(1.0, abs=1e-12)


def test_either_curve_may_wait():
    """Steps offset by one sample: the matching lets one curve hold still for a sample."""
    x = np.linspace(0, 1, 11)
    a = (np.arange(11) >= 6).astype(float)
    b = (np.arange(11) >= 5).astype(float)
    assert discrete_frechet(np.c_[x, a], np.c_[x, b]) == pytest.approx(0.1, abs=1e-12)
    assert discrete_frechet(np.c_[x, b], np.c_[x, a]) == pytest.approx(0.1, abs=1e-12)


def test_the_distance_is_at_least_the_separation_of_either_end():
    rng = np.random.default_rng(5)
    for _ in range(50):
        p, q = rng.normal(size=(12, 2)), rng.normal(size=(9, 2))
        ends = max(np.hypot(*(p[0] - q[0])), np.hypot(*(p[-1] - q[-1])))
        assert discrete_frechet(p, q) >= ends - 1e-12


def test_the_path_is_a_valid_coupling_whose_longest_link_is_the_distance():
    rng = np.random.default_rng(9)
    p, q = rng.normal(size=(17, 2)), rng.normal(size=(23, 2))
    r = frechet_coupling(p, q)
    assert isinstance(r, FrechetResult)
    assert tuple(r.path[0]) == (0, 0)
    assert tuple(r.path[-1]) == (16, 22)
    steps = {tuple(s) for s in np.diff(r.path, axis=0)}
    assert steps <= {(1, 0), (0, 1), (1, 1)}
    expected = np.hypot(p[r.path[:, 0], 0] - q[r.path[:, 1], 0], p[r.path[:, 0], 1] - q[r.path[:, 1], 1])
    assert np.allclose(r.link_lengths, expected)
    assert r.link_lengths.max() == pytest.approx(r.distance, abs=1e-12)
    assert r.link_lengths[r.bottleneck] == pytest.approx(r.distance, abs=1e-12)


@pytest.mark.parametrize("bad", [np.zeros((0, 2)), np.zeros((5, 3)), np.zeros(5), np.array([[0.0, np.nan]])])
def test_bad_curves_are_refused(bad):
    with pytest.raises(ValueError):
        discrete_frechet(bad, np.zeros((3, 2)))
