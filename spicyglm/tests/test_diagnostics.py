"""Diagnostics against the R implementation, plus the exact-sum identities."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from r_compare import DIAG_TABLES, compare_diagnostic_table
from spicyglm import spicy_glm
from spicyglm.diagnostics import wilson_interval

REF = Path(__file__).parent / "r_reference"
CASES = {
    "subject_convex_firth": dict(subject="subject", window="convex", r=40),
    "image_rectangle_firth": dict(subject=None, window="rectangle", r=30),
}


@pytest.fixture(scope="module")
def outputs():
    if not (REF / "cells.csv").exists():
        pytest.skip("run tests/r_reference/make_fixtures.R first")
    cells = pd.read_csv(REF / "cells.csv")
    return {case: spicy_glm(cells, condition="condition", compute_diagnostics=True, **kw)
            for case, kw in CASES.items()}


@pytest.mark.parametrize("table", DIAG_TABLES)
@pytest.mark.parametrize("case", CASES)
def test_tables_match_r(outputs, case, table):
    compare_diagnostic_table(outputs[case], REF, case, table)


def test_share_identities(outputs):
    d = outputs["subject_convex_firth"].diagnostics
    patient, image = d["patient"], d["image"]
    # influence and leverage shares sum to one (Props. 14-16)
    np.testing.assert_allclose(patient.groupby(["from", "to"])["influence_i"].sum(), 1.0)
    np.testing.assert_allclose(patient.groupby(["from", "to", "group"])["l_i"].sum(), 1.0)
    np.testing.assert_allclose(image.groupby(["from", "to", "patient_id"])["l_ij"].sum(), 1.0)
    np.testing.assert_allclose(image.groupby(["from", "to"])["influence_ij"].sum(), 1.0)


def test_mle_disables_diagnostics():
    cells = pd.DataFrame({"imageID": ["a"] * 4 + ["b"] * 4 + ["c"] * 4 + ["d"] * 4,
                          "condition": ["x"] * 8 + ["y"] * 8,
                          "cellType": ["A", "B"] * 8,
                          "x": np.arange(16.0), "y": np.arange(16.0) % 3})
    with pytest.warns(UserWarning):
        out = spicy_glm(cells, condition="condition", r=3, estimator="mle", compute_diagnostics=True)
    assert out.diagnostics is None


def test_wilson_matches_binom():
    # binom::binom.wilson(c(0, 1, 3, 5), c(5, 5, 7, 5))
    lower, upper = wilson_interval([0, 1, 3, 5], [5, 5, 7, 5])
    np.testing.assert_allclose(lower, [0.0, 0.03622411, 0.1582199, 0.5655175], atol=1e-7)
    np.testing.assert_allclose(upper, [0.4344825, 0.6244654, 0.7495416, 1.0], atol=1e-7)
