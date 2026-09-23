"""ref= chooses the reference condition, as spicyR's spicyGLM(ref=)."""

import numpy as np
import pandas as pd
import pytest

from spicyglm import spicy_glm

CASES = {"poisson": dict(r=60), "binomial": dict(family="binomial", k=10)}
EFFECT = {"poisson": "log_rate_ratio", "binomial": "log_odds_ratio"}


def _cells(seed=1, n_img=8, per_type=60):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_img):
        n = per_type * 3
        rows.append(pd.DataFrame({"x": rng.uniform(0, 1000, n), "y": rng.uniform(0, 1000, n),
                                  "cellType": np.repeat(["A", "B", "C"], per_type), "imageID": f"i{i}",
                                  "response": "beta" if i % 2 else "alpha"}))
    return pd.concat(rows, ignore_index=True)


def _by_pair(res):
    return res.set_index(["from", "to"]).sort_index()


@pytest.mark.parametrize("case", CASES)
def test_ref_default_level_is_a_no_op(case):
    cells = _cells()
    default = spicy_glm(cells, condition="response", **CASES[case])
    explicit = spicy_glm(cells, condition="response", ref="alpha", **CASES[case])
    pd.testing.assert_frame_equal(default.results, explicit.results)


@pytest.mark.parametrize("case", CASES)
def test_ref_other_level_swaps_the_contrast(case):
    cells = _cells()
    a = _by_pair(spicy_glm(cells, condition="response", **CASES[case]).results)
    b = _by_pair(spicy_glm(cells, condition="response", ref="beta", **CASES[case]).results)
    assert (b["condition_ref"] == "beta").all() and (b["condition_comp"] == "alpha").all()
    np.testing.assert_allclose(b["coef_ref"], a["coef_comp"], rtol=1e-12)
    np.testing.assert_allclose(b["coef_comp"], a["coef_ref"], rtol=1e-12)
    np.testing.assert_allclose(b[EFFECT[case]], -a[EFFECT[case]], rtol=1e-12)
    np.testing.assert_allclose(b["p_value"], a["p_value"], rtol=1e-10)


def test_ref_overrides_category_order():
    cells = _cells()
    cells["response"] = pd.Categorical(cells["response"], categories=["beta", "alpha"])
    out = spicy_glm(cells, condition="response", r=60, ref="alpha").results
    assert (out["condition_ref"] == "alpha").all()


def test_ref_swap_leaves_patient_diagnostics_unchanged():
    cells = _cells()
    a = spicy_glm(cells, condition="response", r=60, compute_diagnostics=True).diagnostics["patient"]
    b = spicy_glm(cells, condition="response", r=60, compute_diagnostics=True, ref="beta").diagnostics["patient"]
    m = a.merge(b, on=["from", "to", "patient_id"], suffixes=("_a", "_b"))
    assert len(m) == len(a) == len(b)
    assert (m["group_a"] == m["group_b"]).all()
    # delta_i is the shift in the patient's own group coefficient, so it does not depend on the reference
    for col in ["l_i", "influence_i", "delta_i"]:
        np.testing.assert_allclose(m[col + "_b"], m[col + "_a"], rtol=1e-10, err_msg=col)


def test_ref_not_a_level_raises():
    with pytest.raises(ValueError, match="not a level of `condition`"):
        spicy_glm(_cells(), condition="response", r=60, ref="gamma")
