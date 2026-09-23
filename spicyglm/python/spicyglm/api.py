"""Python front end for spicyGLM: pair enumeration, skip rules and result tables.

Numerical work (window areas, neighbour search, fitting, CR2 and degrees of
freedom) happens in the C++ core, ``spicyglm._core``.
"""

import warnings
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from itertools import combinations, product
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import stats

from . import _core
from . import diagnostics as diag

EFFECT_COLUMNS = {"poisson": ("log_rate_ratio", "rate_ratio"), "binomial": ("log_odds_ratio", "odds_ratio")}
SKIP_COLUMNS = ["from", "to", "reason", "message"]


def result_columns(family):
    log_effect, effect = EFFECT_COLUMNS[family]
    return ["from", "to", "condition_ref", "condition_comp", "coef_ref", "coef_comp", log_effect, effect,
            "p_value", "estimator", "family", "mle_would_skip", "mle_skip_reason", "p_adj"]


@dataclass
class SpicyGLMResult:
    results: pd.DataFrame
    skipped: pd.DataFrame
    # None unless compute_diagnostics; otherwise keys "pair", "patient", "image"
    # and "cross_pair" (a dict with "patient" and "image")
    diagnostics: dict | None = None


def spicy_glm(cells, condition, r=None, subject=None, image_id="imageID", cell_type="cellType",
              spatial_coords=("x", "y"), from_=None, to=None, window="convex", family="poisson",
              k=None, estimator="firth", cr2_method="fast", compute_diagnostics=False, top_percent=0.05,
              n_jobs=1):
    """Test for a change in co-localisation of cell-type pairs between two conditions.

    ``family="poisson"``: for each reference cell (``from_``) the number of
    target cells (``to``) within radius ``r`` is Poisson with a log-density
    offset. ``family="binomial"``: the number of target cells among its ``k``
    nearest neighbours is Binomial with a logit background-proportion offset.
    Either way there is one coefficient per condition, and the log rate (odds)
    ratio is tested with the CR2 cluster-robust variance (clustered by
    ``subject``, or by image when ``subject`` is missing or one-to-one with
    images) and Satterthwaite degrees of freedom, or with the model-based
    variance (``cr2_method="naive"``).

    Parameters
    ----------
    cells : pandas.DataFrame
        One row per cell.
    condition : str
        Column with the image-level condition. Exactly two levels are required;
        the first level (category order, otherwise sorted) is the reference.
    r : float
        Neighbourhood radius (poisson only).
    from_, to : str or list of str, optional
        Cell types. A single ``from_`` and ``to`` fits that pair only. Otherwise,
        for ``family="poisson"`` (direction-invariant) every pair among the given
        (or all) cell types is fitted, one direction per unordered pair plus
        self-pairs; for ``family="binomial"`` (directional: from->to and to->from
        differ) every ordered pair in ``from_`` x ``to`` is fitted, where an
        omitted side means all cell types.
    window : {"convex", "rectangle"}
        Observation window used for each image's area (poisson only).
    family : {"poisson", "binomial"}
    k : int
        Number of nearest neighbours (binomial only).
    estimator : {"firth", "mle"}
    cr2_method : {"fast", "naive"}
        ``"fast"``: closed-form CR2 variance, t-test with Satterthwaite degrees
        of freedom. ``"naive"``: model-based variance ignoring clustering,
        z-test; pairs are not skipped for having a single cluster in a condition.
    compute_diagnostics : bool
        Add per-patient and per-image leverage, influence and point-estimate
        shift, their within-pair percentile ranks, and cross-pair flagging.
        Requires ``family="poisson"``, ``estimator="firth"`` and ``cr2_method="fast"``.
    top_percent : float
        Fraction of the within-pair percentile rank counted as flagged.
    n_jobs : int
        Threads used to fit pairs concurrently. Memory grows with n_jobs only by
        the per-pair working data; the cell index is shared.
    """
    if estimator not in ("firth", "mle"):
        raise ValueError("estimator must be 'firth' or 'mle'")
    if cr2_method not in ("fast", "naive"):
        raise ValueError("cr2_method must be 'fast' or 'naive'")
    if int(n_jobs) != n_jobs or n_jobs < 1:
        raise ValueError("n_jobs must be a positive integer")
    if family == "poisson":
        if r is None or not r > 0:
            raise ValueError("family='poisson' requires a positive radius r")
    elif family == "binomial":
        if k is None or int(k) != k or k < 1:
            raise ValueError("family='binomial' requires k, a positive integer")
        k = int(k)
    else:
        raise ValueError("family must be 'poisson' or 'binomial'")
    if compute_diagnostics and (family != "poisson" or estimator != "firth" or cr2_method != "fast"):
        warnings.warn("compute_diagnostics requires family='poisson', estimator='firth' and "
                      "cr2_method='fast'; diagnostics are not computed.")
        compute_diagnostics = False
    x_col, y_col = spatial_coords
    needed = [condition, image_id, cell_type, x_col, y_col] + ([subject] if subject else [])
    missing = [c for c in needed if c not in cells.columns]
    if missing:
        raise KeyError(f"columns not found in cells: {missing}")

    # sort cells by image so each image is a contiguous block for the core
    image_codes, image_labels = pd.factorize(cells[image_id], sort=True)
    order = np.argsort(image_codes, kind="stable")
    df = cells.iloc[order].reset_index(drop=True)
    image_codes = image_codes[order]

    levels = _condition_levels(cells[condition])
    if len(levels) != 2:
        raise ValueError(f"spicyGLM compares exactly two conditions; found {len(levels)}: {list(levels)}")
    if (df.groupby(image_codes)[condition].nunique() > 1).any():
        raise ValueError(f"'{condition}' must be constant within each image")
    if subject and (cells.groupby(subject)[condition].nunique() > 1).any():
        raise ValueError(f"each subject must belong to a single '{condition}' level")
    image_condition = df.groupby(image_codes, sort=True)[condition].first().to_numpy()

    one_to_one = subject is None or cells[subject].nunique() == cells[image_id].nunique()
    if one_to_one:
        image_cluster, cluster_labels = np.arange(len(image_labels)), image_labels
    else:
        image_cluster, cluster_labels = pd.factorize(df.groupby(image_codes, sort=True)[subject].first())

    type_labels = pd.unique(cells[cell_type])  # first-appearance order, as R's unique()
    type_codes = pd.Categorical(df[cell_type], categories=type_labels).codes.astype(np.int32)
    image_offsets = np.searchsorted(image_codes, np.arange(len(image_labels) + 1)).astype(np.int32)
    ctx = SimpleNamespace(
        family=family, k=k, estimator=estimator, cr2_method=cr2_method, levels=levels,
        image_group=pd.Categorical(image_condition, categories=levels).codes,
        image_cluster=image_cluster,
        type_index={t: i for i, t in enumerate(type_labels)},
        data=_core.Dataset(df[x_col].to_numpy(dtype=float), df[y_col].to_numpy(dtype=float), type_codes,
                           image_offsets, len(type_labels)),
    )
    ctx.presence = _presence(type_codes, image_codes, ctx.image_group, len(type_labels))
    if family == "poisson":
        ctx.areas = ctx.data.image_areas(window)
        ctx.data.build_radius_index(r)
    else:
        ctx.data.build_knn(k, n_jobs)

    pairs = _pairs(from_, to, list(type_labels), family)
    for f, t in pairs:
        if f not in ctx.type_index or t not in ctx.type_index:
            raise KeyError(f"cell type not found: {f if f not in ctx.type_index else t}")
    if n_jobs == 1:
        outcomes = [_fit_pair(ctx, f, t, compute_diagnostics) for f, t in pairs]
    else:
        with ThreadPoolExecutor(max_workers=n_jobs) as pool:
            outcomes = list(pool.map(lambda p: _fit_pair(ctx, p[0], p[1], compute_diagnostics), pairs))

    fitted, skipped, pair_diag = [], [], []
    for (f, t), outcome in zip(pairs, outcomes):
        if "reason" in outcome:
            skipped.append(outcome)
            continue
        fit = outcome.pop("_fit")
        fitted.append(outcome)
        if compute_diagnostics:
            pair_diag.append(diag.pair_tables(fit, f, t, levels, cluster_labels, image_labels))

    results = pd.DataFrame(fitted, columns=result_columns(family)[:-1])
    results["p_adj"] = _bh(results["p_value"].to_numpy(dtype=float))
    results = results.sort_values("p_adj", kind="stable").reset_index(drop=True)
    diagnostics = _assemble_diagnostics(pair_diag, top_percent) if pair_diag else None
    return SpicyGLMResult(results=results, skipped=pd.DataFrame(skipped, columns=SKIP_COLUMNS),
                          diagnostics=diagnostics)


def _assemble_diagnostics(pair_diag, top_percent):
    patient = diag.relative_diagnostics(pd.concat([p for _, p, _ in pair_diag], ignore_index=True), "patient")
    image = diag.relative_diagnostics(pd.concat([i for _, _, i in pair_diag], ignore_index=True), "image")
    return {
        "pair": pd.DataFrame([s for s, _, _ in pair_diag]),
        "patient": patient,
        "image": image,
        "cross_pair": {
            "patient": diag.cross_pair_diagnostics(patient, top_percent),
            "image": diag.cross_pair_diagnostics(image, top_percent),
        },
    }


def _fit_pair(ctx, f, t, compute_diagnostics):
    from_code, to_code = ctx.type_index[f], ctx.type_index[t]
    if ctx.family == "poisson":
        md = ctx.data.poisson_model_data(ctx.areas, from_code, to_code)
    else:
        md = ctx.data.binomial_model_data(from_code, to_code)
    group = ctx.image_group[md["image"]]
    present = set(np.unique(group))
    if len(present) < 2:
        missing = [g for g in range(2) if g not in present]
        return _skip(f, t, *_diagnose_missing(f, t, missing, ctx))

    n = md["n"]
    mle_skip_reason, message = _boundary(f, t, n, group, ctx)
    if ctx.estimator == "mle" and mle_skip_reason:
        return _skip(f, t, mle_skip_reason, message)

    cluster = ctx.image_cluster[md["image"]]
    for g in range(2 if ctx.cr2_method == "fast" else 0):  # as R: only CR2 needs two clusters per group
        if np.unique(cluster[group == g]).size < 2:
            return _skip(f, t, "one_patient_per_group",
                         f"Skipping pair {f}__{t}: condition '{ctx.levels[g]}' has fewer than two "
                         f"clusters with data for this pair; CR2 needs at least two per group.")

    if ctx.family == "poisson":
        fit = _core.fit_pair_poisson(cluster, md["image"], group, n, md["density"], ctx.estimator,
                                     ctx.cr2_method, compute_diagnostics)
    else:
        fit = _core.fit_pair_binomial(cluster, md["image"], group, n, ctx.k, md["p0"], ctx.estimator,
                                      ctx.cr2_method)
    b_ref, b_comp = fit["beta"]
    log_effect = b_comp - b_ref
    t_stat = log_effect / np.sqrt(fit["v_hat"])
    if ctx.cr2_method == "naive":
        p_value = 2 * stats.norm.sf(abs(t_stat))
    else:
        p_value = 2 * stats.t.sf(abs(t_stat), fit["df"])
    log_col, col = EFFECT_COLUMNS[ctx.family]
    return {
        "_fit": fit,
        "from": f, "to": t, "condition_ref": ctx.levels[0], "condition_comp": ctx.levels[1],
        "coef_ref": b_ref, "coef_comp": b_comp, log_col: log_effect, col: np.exp(log_effect),
        "p_value": p_value, "estimator": ctx.estimator, "family": ctx.family,
        "mle_would_skip": mle_skip_reason is not None, "mle_skip_reason": mle_skip_reason,
    }


def _boundary(f, t, n, group, ctx):
    """Conditions where the MLE is infinite: all counts zero, or (binomial) all k."""
    floor = [ctx.levels[g] for g in range(2) if not np.any(n[group == g] > 0)]
    ceil = [ctx.levels[g] for g in range(2) if ctx.family == "binomial" and np.all(n[group == g] == ctx.k)]
    boundary = floor + [c for c in ceil if c not in floor]
    if not boundary:
        return None, None
    if len(boundary) == 2:
        reason = "all_zero" if not ceil else "all_max" if not floor else "all_boundary"
    else:
        reason = "one_condition_zero" if floor else "one_condition_max"
    return reason, (f"Skipping pair {f}__{t}: neighbour counts sit at a separation boundary in condition(s) "
                    f"{boundary} (reason '{reason}'), so the MLE is not finite. Use estimator='firth'.")


def _condition_levels(col):
    if isinstance(col.dtype, pd.CategoricalDtype):
        return [c for c in col.cat.categories if (col == c).any()]
    return sorted(col.dropna().unique())


def _unique_types(group):
    types = []
    for ct in ([group] if isinstance(group, str) else (group or [])):
        if ct not in types:
            types.append(ct)
    return types


def _pairs(from_, to, all_types, family):
    if isinstance(from_, str) and isinstance(to, str):
        return [(from_, to)]
    if family == "binomial":
        # the kNN effect is directional (A->B != B->A), so fit every ordered pair
        return list(product(_unique_types(from_) if from_ is not None else all_types,
                            _unique_types(to) if to is not None else all_types))
    if from_ is not None or to is not None:
        types = _unique_types(from_) + [ct for ct in _unique_types(to) if ct not in _unique_types(from_)]
    else:
        types = all_types
    return list(combinations(types, 2)) + [(ct, ct) for ct in types]


def _presence(type_codes, image_codes, image_group, n_types):
    """present[type, group] = number of images in the group containing the type."""
    has = np.zeros((n_types, image_group.size), dtype=bool)
    has[type_codes, image_codes] = True
    return np.stack([has[:, image_group == g].sum(axis=1) for g in range(2)], axis=1)


def _diagnose_missing(f, t, missing_groups, ctx):
    reasons, texts = [], []
    for g in missing_groups:
        level = ctx.levels[g]
        f_in, t_in = ctx.presence[ctx.type_index[f], g] > 0, ctx.presence[ctx.type_index[t], g] > 0
        if not f_in and not t_in:
            reasons.append("both_absent")
            texts.append(f"condition '{level}' has no images containing '{f}' or '{t}'.")
        elif not (f_in and t_in):
            reasons.append("type_absent")
            texts.append(f"condition '{level}' has no images containing '{t if f_in else f}'.")
        else:
            reasons.append("no_cooccurrence")
            texts.append(f"condition '{level}' has no image containing both '{f}' and '{t}'.")
    reason = reasons[0] if len(set(reasons)) == 1 else "mixed"
    return reason, f"Skipping pair {f}__{t}: " + " ".join(texts)


def _skip(f, t, reason, message):
    return {"from": f, "to": t, "reason": reason, "message": message}


def _bh(p):
    """Benjamini-Hochberg adjustment ignoring NaNs, as R's p.adjust(method = 'fdr')."""
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        out[ok] = stats.false_discovery_control(p[ok], method="bh")
    return out
