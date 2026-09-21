"""Compare spicy_glm output with CSVs written by tests/r_reference/run_spicyglm.R."""

from pathlib import Path

import numpy as np
import pandas as pd

R_EFFECT = {"poisson": ("log_rate_ratio", "logRateRatio"), "binomial": ("log_odds_ratio", "logOddsRatio")}
# Binomial reference fits are refitted to epsilon = 1e-14 (tighten_binomial_fits
# in run_spicyglm.R), so both families are compared at the same tolerance.
TOL = {"poisson": dict(rtol=1e-7, atol=1e-12), "binomial": dict(rtol=1e-7, atol=1e-12)}
DIAG_TABLES = {
    "pair": lambda d: d["pair"],
    "patient": lambda d: d["patient"],
    "image": lambda d: d["image"],
    "cross_patient": lambda d: d["cross_pair"]["patient"],
    "cross_image": lambda d: d["cross_pair"]["image"],
}


def compare_results(out, r_dir, name, family):
    """Assert fitted pairs, estimates, p-values and skip reasons match R.

    R's buildGLM() takes the reference condition from the first image with data
    for a pair (dplyr::bind_rows merges factor levels in order of appearance),
    so on some pairs R's reference is the second level. spicy_glm always uses
    the first level; on those pairs the R coefficients are swapped and the log
    ratio negated before comparing. Returns (number of pairs, pairs where R
    flipped the reference, pairs whose R fit diverged).
    """
    r_dir = Path(r_dir)
    r_res = pd.read_csv(r_dir / f"results_{name}.csv")
    r_res = r_res.rename(columns=lambda c: c if c in ("from", "to") else "R_" + c)
    merged = out.results.merge(r_res, on=["from", "to"], how="outer", indicator=True)
    unmatched = merged.loc[merged["_merge"] != "both", ["from", "to", "_merge"]]
    assert unmatched.empty, f"pairs fitted by only one implementation:\n{unmatched}"

    flipped = (merged["R_conditionRef"] != merged["condition_ref"]).to_numpy()
    assert (merged.loc[flipped, "R_conditionRef"] == merged.loc[flipped, "condition_comp"]).all()
    py_log, r_log = R_EFFECT[family]
    r_coef_ref = np.where(flipped, merged["R_coef_comp"], merged["R_coef_ref"])
    r_coef_comp = np.where(flipped, merged["R_coef_ref"], merged["R_coef_comp"])
    r_log_effect = np.where(flipped, -merged["R_" + r_log], merged["R_" + r_log])

    # brglmFit can fail to converge on separated binomial pairs; buildGLM() does
    # not check fit$converged and reports coefficients near -1e15. Those pairs
    # are excluded, and since their spurious p-values shift every BH adjustment,
    # p_adj is then checked against BH of our own p-values instead of R's.
    diverged = (np.abs(r_coef_ref) > 1e6) | (np.abs(r_coef_comp) > 1e6)
    ok = ~diverged
    for col, expected in [("coef_ref", r_coef_ref), ("coef_comp", r_coef_comp), (py_log, r_log_effect),
                          ("p_value", merged["R_p.value"])]:
        np.testing.assert_allclose(merged.loc[ok, col], np.asarray(expected)[ok], err_msg=col, **TOL[family])
    if diverged.any():
        from scipy.stats import false_discovery_control
        np.testing.assert_allclose(merged["p_adj"], false_discovery_control(merged["p_value"], method="bh"))
    else:
        np.testing.assert_allclose(merged["p_adj"], merged["R_p.adj"], err_msg="p_adj", **TOL[family])
    assert (merged["mle_would_skip"] == merged["R_mle_would_skip"]).all(), "mle_would_skip"

    skip_file = r_dir / f"skipped_{name}.csv"
    r_skip = pd.read_csv(skip_file) if skip_file.stat().st_size > 1 else pd.DataFrame(columns=["from", "to", "reason"])
    py_skip = out.skipped.set_index(["from", "to"])["reason"].sort_index()
    r_skip = r_skip.set_index(["from", "to"])["reason"].sort_index()
    assert py_skip.to_dict() == r_skip.to_dict(), "skip reasons differ"
    flipped_pairs = set(map(tuple, merged.loc[flipped, ["from", "to"]].to_numpy()))
    diverged_pairs = set(map(tuple, merged.loc[diverged, ["from", "to"]].to_numpy()))
    return len(merged), flipped_pairs, diverged_pairs


def compare_diagnostic_table(out, r_dir, name, table, flipped_pairs=frozenset()):
    """Assert one diagnostics table matches R column for column, row for row.

    e_i = w_g s_i with w = (-1/S_1, 1/S_2), so its sign (and e_ij's) follows the
    reference condition; it is negated for pairs where R flipped the reference.
    """
    py = DIAG_TABLES[table](out.diagnostics).reset_index(drop=True)
    r = pd.read_csv(Path(r_dir) / f"diag_{table}_{name}.csv")
    assert list(py.columns) == list(r.columns), f"{table}: columns differ"
    assert len(py) == len(r), f"{table}: {len(py)} rows vs {len(r)} in R"
    if table.startswith("cross"):
        # Rows are sorted by wilson_lower_influence. With no flagged pairs its true
        # value is 0, but R's arithmetic leaves noise around 1e-18 that varies with
        # n, so R's order among those ties is arbitrary: check our order is
        # non-increasing, then compare rows matched by id.
        assert (np.diff(py["wilson_lower_influence"].to_numpy()) <= 1e-12).all(), f"{table}: not sorted"
        keys = [c for c in ("patient_id", "image_id") if c in r.columns]
        py = py.sort_values(keys, kind="stable").reset_index(drop=True)
        r = r.assign(**{c: r[c].astype(str) for c in keys}).sort_values(keys, kind="stable").reset_index(drop=True)
    if flipped_pairs and {"from", "to"} <= set(r.columns):
        in_flipped = np.array([(f, t) in flipped_pairs for f, t in zip(r["from"], r["to"])])
        for col in ("e_i", "e_ij"):
            if col in r.columns:
                r[col] = np.where(in_flipped, -r[col], r[col])
    for col in r.columns:
        if pd.api.types.is_numeric_dtype(r[col]) and not pd.api.types.is_bool_dtype(r[col]):
            np.testing.assert_allclose(py[col].to_numpy(dtype=float), r[col].to_numpy(dtype=float),
                                       rtol=1e-7, atol=1e-12, err_msg=f"{table}.{col}")
        else:
            assert (py[col].astype(str) == r[col].astype(str)).all(), f"{table}.{col}"
    return len(r)
