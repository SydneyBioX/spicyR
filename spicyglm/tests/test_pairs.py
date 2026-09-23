"""Pair enumeration: poisson is direction-invariant, binomial is directional."""

import numpy as np
import pandas as pd
import pytest

from spicyglm import spicy_glm
from spicyglm.api import _pairs

TYPES = ["A", "B", "C"]


def _cells(seed=1, n_img=8, per_type=60):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_img):
        n = per_type * len(TYPES)
        rows.append(pd.DataFrame({"x": rng.uniform(0, 1000, n), "y": rng.uniform(0, 1000, n),
                                  "cellType": np.repeat(TYPES, per_type), "imageID": f"i{i}",
                                  "response": "beta" if i % 2 else "alpha"}))
    return pd.concat(rows, ignore_index=True)


def test_poisson_pairs_are_unordered():
    assert _pairs(None, None, TYPES, "poisson") == [("A", "B"), ("A", "C"), ("B", "C"),
                                                     ("A", "A"), ("B", "B"), ("C", "C")]
    assert _pairs("A", ["B", "C"], TYPES, "poisson") == [("A", "B"), ("A", "C"), ("B", "C"),
                                                          ("A", "A"), ("B", "B"), ("C", "C")]


def test_binomial_pairs_are_ordered():
    assert sorted(_pairs(None, None, TYPES, "binomial")) == sorted((f, t) for f in TYPES for t in TYPES)
    assert _pairs("A", ["B", "C"], TYPES, "binomial") == [("A", "B"), ("A", "C")]
    assert _pairs(None, "A", TYPES, "binomial") == [("A", "A"), ("B", "A"), ("C", "A")]


@pytest.mark.parametrize("family", ["poisson", "binomial"])
def test_single_pair_is_fitted_as_given(family):
    assert _pairs("B", "A", TYPES, family) == [("B", "A")]


def test_binomial_directions_differ_and_match_single_pair_fits():
    cells = _cells()
    out = spicy_glm(cells, condition="response", family="binomial", k=10).results.set_index(["from", "to"])
    assert len(out) == 9
    assert abs(out.loc[("A", "B"), "log_odds_ratio"] - out.loc[("B", "A"), "log_odds_ratio"]) > 1e-3
    for f, t in [("A", "B"), ("B", "A")]:
        one = spicy_glm(cells, condition="response", family="binomial", k=10, from_=f, to=t).results.iloc[0]
        assert out.loc[(f, t), "log_odds_ratio"] == pytest.approx(one["log_odds_ratio"], rel=1e-12)
        assert out.loc[(f, t), "p_value"] == pytest.approx(one["p_value"], rel=1e-12)


def test_poisson_is_direction_invariant():
    cells = _cells()
    ab = spicy_glm(cells, condition="response", r=60, from_="A", to="B").results.iloc[0]
    ba = spicy_glm(cells, condition="response", r=60, from_="B", to="A").results.iloc[0]
    assert ab["log_rate_ratio"] == pytest.approx(ba["log_rate_ratio"], rel=1e-12)
    assert ab["p_value"] == pytest.approx(ba["p_value"], rel=1e-10)
