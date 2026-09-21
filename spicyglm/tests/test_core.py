"""Check the C++ core against dense implementations built from the definitions."""

import numpy as np
import pytest
from scipy.spatial import ConvexHull

from spicyglm import _core


def _inv_sqrt(M):
    vals, vecs = np.linalg.eigh(M)
    return vecs @ np.diag(vals ** -0.5) @ vecs.T


def _dense_cr2(cluster, group, var, resid):
    """CR2 v_hat and Satterthwaite df with N x N matrices (Sections 4 and 9)."""
    N = len(cluster)
    X = np.zeros((N, 2))
    X[np.arange(N), group] = 1.0
    V = np.diag(var)
    B = np.linalg.inv(X.T @ V @ X)
    IH = np.eye(N) - V @ X @ B @ X.T
    L = np.array([-1.0, 1.0])
    loadings, e = [], []
    for c in dict.fromkeys(cluster):
        idx = np.flatnonzero(cluster == c)
        Q = IH[idx] @ V @ IH[idx].T
        D = np.diag(np.sqrt(var[idx]))
        A = D @ _inv_sqrt(D @ Q @ D) @ D
        row = L @ B @ X[idx].T @ A
        e.append(row @ resid[idx])
        loadings.append(row @ IH[idx])
    G = np.array(loadings).T
    P = G.T @ V @ G
    return float(np.sum(np.square(e))), np.trace(P) ** 2 / np.sum(P ** 2)


def _random_design(rng, one_image_per_cluster=False):
    cluster, image, group, var = [], [], [], []
    img = 0
    for c in range(rng.integers(4, 9)):
        g = c % 2
        for _ in range(1 if one_image_per_cluster else rng.integers(1, 4)):
            n_cells = rng.integers(2, 10)
            v = rng.uniform(0.2, 4.0)
            cluster += [c] * n_cells
            image += [img] * n_cells
            group += [g] * n_cells
            var += [v] * n_cells
            img += 1
    var = np.array(var)
    resid = rng.poisson(var) - var
    return np.array(cluster), np.array(image), np.array(group), var, resid


@pytest.mark.parametrize("seed", range(6))
@pytest.mark.parametrize("one_image", [False, True])
def test_cr2_matches_dense(seed, one_image):
    rng = np.random.default_rng(seed)
    cluster, image, group, var, resid = _random_design(rng, one_image)
    out = _core.cr2_wald(cluster, image, group, var, resid)
    v_hat, df = _dense_cr2(cluster, group, var, resid)
    assert out["v_hat"] == pytest.approx(v_hat, rel=1e-10)
    assert out["df"] == pytest.approx(df, rel=1e-10)
    assert sum(np.square(out["e"])) == pytest.approx(out["v_hat"], rel=1e-12)


def test_cr2_rejects_single_cluster_group():
    cluster = np.array([0, 0, 1, 1])
    group = np.array([0, 0, 1, 1])
    var = np.ones(4)
    with pytest.raises(RuntimeError):
        _core.cr2_wald(cluster, cluster, group, var, np.zeros(4))


def test_fit_poisson_closed_forms():
    n = np.array([0, 0, 3, 5])
    density = np.array([1.0, 2.0, 1.5, 0.5])
    group = np.array([0, 0, 1, 1])
    firth = _core.fit_poisson(n, density, group, "firth")
    assert firth["beta"] == pytest.approx((np.log(0.5 / 3.0), np.log(8.5 / 2.0)))
    assert firth["mu"] == pytest.approx(np.exp(np.array(firth["beta"])[group]) * density)


def _images(rng):
    """Images with varied geometry: square, a thin strip, collinear, duplicates, one cell."""
    square = rng.uniform(0, 100, size=(400, 2))
    strip = rng.uniform(0, 1, size=(300, 2)) * [5000, 3]
    line = np.column_stack([np.linspace(0, 900, 150), np.full(150, 7.0)])
    dupes = np.repeat(rng.uniform(0, 30, size=(40, 2)), 3, axis=0)
    single = np.array([[1.0, 1.0]])
    return [square, strip, line, dupes, single]


@pytest.mark.parametrize("r", [0.01, 3.0, 25.0, 1e4])
def test_dataset_radius_counts_match_brute_force(r):
    rng = np.random.default_rng(int(r * 100))
    images = _images(rng)
    xy = np.vstack(images)
    offsets = np.concatenate([[0], np.cumsum([len(i) for i in images])]).astype(np.int32)
    types = rng.integers(0, 3, size=len(xy)).astype(np.int32)
    data = _core.Dataset(xy[:, 0], xy[:, 1], types, offsets, 3)
    data.build_radius_index(r)
    areas = data.image_areas("rectangle")
    for f, t in [(0, 1), (1, 1), (2, 0)]:
        md = data.poisson_model_data(areas, f, t)
        expected_rows, expected_n = [], []
        for img in range(len(images)):
            rows = np.arange(offsets[img], offsets[img + 1])
            ref, tgt = rows[types[rows] == f], rows[types[rows] == t]
            if len(ref) == 0 or len(tgt) == 0:
                continue
            d = np.linalg.norm(xy[ref, None, :] - xy[None, tgt, :], axis=2)
            expected_rows += list(ref)
            expected_n += list((d <= r).sum(axis=1))
        np.testing.assert_array_equal(md["row"], expected_rows)
        np.testing.assert_array_equal(md["n"], expected_n)


@pytest.mark.parametrize("k", [1, 5, 12])
def test_knn_matches_brute_force(k):
    rng = np.random.default_rng(k)
    sizes = [40, k, 300, 7]  # images with at most k cells get -1 rows
    pts = [rng.uniform(0, rng.uniform(10, 500), size=(s, 2)) * [1, rng.uniform(0.05, 1)] for s in sizes]
    xy = np.vstack(pts)
    offsets = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int32)
    got = _core.knn_indices(xy[:, 0], xy[:, 1], offsets, k)
    for img, s in enumerate(sizes):
        rows = slice(offsets[img], offsets[img + 1])
        if s <= k:
            assert (got[rows] == -1).all()
            continue
        d = np.linalg.norm(xy[rows, None, :] - xy[None, rows, :], axis=2)
        np.fill_diagonal(d, np.inf)
        expected = np.sort(d, axis=1)[:, :k]
        found = np.sort(np.take_along_axis(d, got[rows] - offsets[img], axis=1), axis=1)
        np.testing.assert_allclose(found, expected)
    np.testing.assert_array_equal(_core.knn_indices(xy[:, 0], xy[:, 1], offsets, k, n_threads=3), got)


def test_knn_breaks_ties_like_spatstat():
    from pathlib import Path

    import pandas as pd

    path = Path(__file__).parent / "r_reference" / "knn_ties.csv"
    if not path.exists():
        pytest.skip("run tests/r_reference/make_fixtures.R first")
    ref = pd.read_csv(path)
    expected = ref.filter(like="nn").to_numpy()
    got = _core.knn_indices(ref["x"].to_numpy(float), ref["y"].to_numpy(float),
                            np.array([0, len(ref)], np.int32), expected.shape[1])
    np.testing.assert_array_equal(got, expected)


def _binomial_objective(beta, n, k, offset, firth):
    """Negative (optionally Jeffreys-penalised) log-likelihood for one group."""
    p = 1 / (1 + np.exp(-(beta + offset)))
    ll = np.sum(n * np.log(p) + (k - n) * np.log1p(-p))
    if firth:
        ll += 0.5 * np.log(np.sum(k * p * (1 - p)))
    return -ll


@pytest.mark.parametrize("estimator", ["mle", "firth"])
@pytest.mark.parametrize("seed", range(3))
def test_fit_binomial_maximises_likelihood(estimator, seed):
    from scipy.optimize import minimize_scalar

    rng = np.random.default_rng(seed)
    k = 8
    image = rng.integers(0, 6, size=120)
    p0 = rng.uniform(0.05, 0.6, size=6)[image]
    group = (image >= 3).astype(int)
    n = rng.binomial(k, p0)
    if estimator == "firth":
        n[group == 0] = 0  # structural zero: finite only under Firth
    fit = _core.fit_binomial(n, k, p0, group, estimator)
    for g in range(2):
        sel = group == g
        offset = np.log(p0[sel] / (1 - p0[sel]))
        best = minimize_scalar(_binomial_objective, bounds=(-30, 30), method="bounded",
                               args=(n[sel], k, offset, estimator == "firth"),
                               options={"xatol": 1e-10})
        assert fit["beta"][g] == pytest.approx(best.x, abs=1e-6)


def test_window_area():
    rng = np.random.default_rng(1)
    pts = rng.normal(size=(500, 2))
    assert _core.window_area(pts[:, 0], pts[:, 1], "convex") == pytest.approx(ConvexHull(pts).volume)
    span = pts.max(axis=0) - pts.min(axis=0)
    assert _core.window_area(pts[:, 0], pts[:, 1], "rectangle") == pytest.approx(span[0] * span[1])
