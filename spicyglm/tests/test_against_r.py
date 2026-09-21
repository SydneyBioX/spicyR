"""End-to-end comparison with the R implementation (fixtures from make_fixtures.R)."""

from pathlib import Path

import pandas as pd
import pytest

from r_compare import compare_results
from spicyglm import spicy_glm

REF = Path(__file__).parent / "r_reference"
CASES = {
    "subject_convex_firth": dict(subject="subject", window="convex", estimator="firth", r=40),
    "image_rectangle_firth": dict(subject=None, window="rectangle", estimator="firth", r=30),
    "subject_convex_mle": dict(subject="subject", window="convex", estimator="mle", r=40),
    "binomial_subject_firth": dict(subject="subject", estimator="firth", family="binomial", k=10),
    "binomial_image_mle": dict(subject=None, estimator="mle", family="binomial", k=6),
    "naive_subject_convex_firth": dict(subject="subject", window="convex", estimator="firth", r=40,
                                       cr2_method="naive"),
    "naive_image_rectangle_mle": dict(subject=None, window="rectangle", estimator="mle", r=30, cr2_method="naive"),
    "naive_binomial_subject_firth": dict(subject="subject", estimator="firth", family="binomial", k=10,
                                         cr2_method="naive"),
}

pytestmark = pytest.mark.skipif(not (REF / "cells.csv").exists(),
                                reason="run tests/r_reference/make_fixtures.R first")


@pytest.fixture(scope="module")
def cells():
    return pd.read_csv(REF / "cells.csv")


@pytest.mark.parametrize("n_jobs", [1, 3])
@pytest.mark.parametrize("case", CASES)
def test_matches_r(cells, case, n_jobs):
    out = spicy_glm(cells, condition="condition", n_jobs=n_jobs, **CASES[case])
    compare_results(out, REF, case, CASES[case].get("family", "poisson"))
