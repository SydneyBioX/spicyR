"""Diagnostic tables built from the per-pair quantities computed in C++.

Mirrors spicyR's assembleDiagnosticsTables(), computePairSummary(),
computeRelativeDiagnostics() and crossPairDiagnostics() (Section 11 of
spicyClub_Supplementary_Math.pdf).
"""

import numpy as np
import pandas as pd
from scipy import stats

PATIENT_COLUMNS = ["patient_id", "group", "n_i", "T_i", "S_g", "l_i", "raw_residual_sum",
                   "adjusted_residual_sum", "e_i", "influence_i", "y_i", "d_i", "delta_i"]
IMAGE_COLUMNS = ["patient_id", "image_id", "group", "n_ij", "density_ij", "l_ij", "l_ij_group_share",
                 "raw_residual_sum_ij", "adjusted_residual_sum_ij", "e_ij", "e_ij_share_within_patient",
                 "influence_ij"]


def pair_tables(fit, f, t, levels, cluster_labels, image_labels):
    """Pair summary row plus patient and image tables for one fitted pair."""
    p = fit["patient"]
    patient = pd.DataFrame({k: p[k] for k in PATIENT_COLUMNS[2:] if k != "S_g"})
    patient["patient_id"] = [str(cluster_labels[c]) for c in p["cluster_id"]]
    patient["group"] = [levels[g] for g in p["group"]]
    patient["S_g"] = np.asarray(p["S_g"])[p["group"]]
    patient = patient[PATIENT_COLUMNS].sort_values("patient_id", kind="stable").reset_index(drop=True)

    im = fit["image"]
    image = pd.DataFrame({k: im[k] for k in IMAGE_COLUMNS[3:]})
    image["patient_id"] = [str(cluster_labels[p["cluster_id"][c]]) for c in im["cluster"]]
    image["image_id"] = [str(image_labels[i]) for i in im["image_id"]]
    image["group"] = [levels[p["group"][c]] for c in im["cluster"]]
    image = image[IMAGE_COLUMNS].sort_values(["patient_id", "image_id"], kind="stable").reset_index(drop=True)

    abs_delta = patient["delta_i"].abs()
    summary = {
        "from": f, "to": t,
        "n_patients": len(patient),
        "nu": fit["df"],
        "max_influence": patient["influence_i"].max(),
        "patient_with_max_influence": patient["patient_id"][patient["influence_i"].idxmax()],
        "max_abs_delta_logRR": abs_delta.max(),
        "patient_with_max_delta_logRR": (patient["patient_id"][abs_delta.idxmax()]
                                         if abs_delta.notna().any() else None),
    }
    for tbl in (patient, image):
        tbl["from"], tbl["to"] = f, t
    return summary, patient, image


def _percent_rank(s):
    """dplyr::percent_rank: (min rank - 1) / (non-missing count - 1)."""
    return (s.rank(method="min") - 1) / (s.notna().sum() - 1)


def relative_diagnostics(table, level):
    out = table.copy()
    pair = ["from", "to"]
    if level == "patient":
        by_group = out.groupby(pair + ["group"], sort=False)
        out["n_patients_per_group"] = by_group["patient_id"].transform("size")
        out["rel_l_i"] = out["l_i"] * out["n_patients_per_group"]
        out["percentile_rank_leverage"] = out.groupby(pair + ["group"], sort=False)["rel_l_i"].transform(_percent_rank)
        by_pair = out.groupby(pair, sort=False)
        out["n_patients_in_pair"] = by_pair["patient_id"].transform("size")
        out["rel_influence_i"] = out["influence_i"] * out["n_patients_in_pair"]
        out["percentile_rank_influence"] = out.groupby(pair, sort=False)["rel_influence_i"].transform(_percent_rank)
        out["abs_delta"] = out["delta_i"].abs()
        out["percentile_rank_delta"] = out.groupby(pair, sort=False)["abs_delta"].transform(_percent_rank)
        return out[["patient_id", "group", "from", "to", "n_patients_per_group", "n_patients_in_pair",
                    "l_i", "influence_i", "delta_i", "rel_l_i", "percentile_rank_leverage",
                    "rel_influence_i", "percentile_rank_influence", "percentile_rank_delta",
                    "n_i", "T_i", "S_g", "raw_residual_sum", "adjusted_residual_sum", "e_i", "y_i", "d_i"]]

    by_patient = out.groupby(pair + ["patient_id"], sort=False)
    out["n_images_for_patient"] = by_patient["image_id"].transform("size")
    out["rel_l_ij"] = out["l_ij"] * out["n_images_for_patient"]
    out["percentile_rank_leverage"] = out.groupby(pair + ["patient_id"], sort=False)["rel_l_ij"].transform(_percent_rank)
    by_pair = out.groupby(pair, sort=False)
    out["n_images_in_pair"] = by_pair["image_id"].transform("size")
    out["rel_influence_ij"] = out["influence_ij"] * out["n_images_in_pair"]
    out["percentile_rank_influence"] = out.groupby(pair, sort=False)["rel_influence_ij"].transform(_percent_rank)
    return out[["patient_id", "image_id", "group", "from", "to", "n_images_for_patient", "n_images_in_pair",
                "l_ij", "l_ij_group_share", "influence_ij", "rel_l_ij", "percentile_rank_leverage",
                "rel_influence_ij", "percentile_rank_influence", "density_ij", "raw_residual_sum_ij",
                "adjusted_residual_sum_ij", "e_ij", "e_ij_share_within_patient"]]


def wilson_interval(x, n, conf_level=0.95):
    """Wilson score interval, as binom::binom.wilson (no clipping)."""
    x, n = np.asarray(x, dtype=float), np.asarray(n, dtype=float)
    z = stats.norm.ppf(1 - (1 - conf_level) / 2)
    z2 = z * z
    p = x / n
    centre = p + 0.5 * z2 / n
    half = z * np.sqrt((p * (1 - p) + 0.25 * z2 / n) / n)
    return (centre - half) / (1 + z2 / n), (centre + half) / (1 + z2 / n)


def cross_pair_diagnostics(relative, top_percent=0.05):
    image_level = "image_id" in relative.columns
    keys = ["patient_id", "image_id"] if image_level else ["patient_id"]
    if (relative.groupby(keys)["group"].nunique() > 1).any():
        raise ValueError("a patient/image's condition group varies across cell-type pairs")

    rel = relative.copy()
    pathways = ["influence", "leverage"] + ([] if image_level else ["delta"])
    for pw in pathways:
        rel[f"flagged_{pw}"] = rel[f"percentile_rank_{pw}"] >= 1 - top_percent

    grouped = rel.groupby(keys, sort=True)
    out = grouped.agg(group=("group", "first"), n_pairs_present=("group", "size"))
    if image_level:
        for new, col in [("mean_influence_share_ij", "e_ij_share_within_patient"),
                         ("mean_influence_ij", "influence_ij"), ("mean_l_ij", "l_ij"),
                         ("mean_l_ij_group_share", "l_ij_group_share")]:
            out[new] = grouped[col].mean()
    for pw in pathways:
        out[f"n_pairs_flagged_{pw}"] = grouped[f"flagged_{pw}"].sum()
    for pw in pathways:
        out[f"prop_flagged_{pw}"] = out[f"n_pairs_flagged_{pw}"] / out["n_pairs_present"]
    for pw in pathways:
        out[f"wilson_lower_{pw}"], out[f"wilson_upper_{pw}"] = wilson_interval(
            out[f"n_pairs_flagged_{pw}"], out["n_pairs_present"])
    out = out.reset_index()

    if image_level:
        out = out[["patient_id", "image_id", "group", "n_pairs_present", "mean_influence_share_ij",
                   "mean_influence_ij", "mean_l_ij", "mean_l_ij_group_share",
                   "n_pairs_flagged_influence", "prop_flagged_influence",
                   "wilson_lower_influence", "wilson_upper_influence",
                   "n_pairs_flagged_leverage", "prop_flagged_leverage",
                   "wilson_lower_leverage", "wilson_upper_leverage"]]
    return out.sort_values("wilson_lower_influence", ascending=False, kind="stable").reset_index(drop=True)
