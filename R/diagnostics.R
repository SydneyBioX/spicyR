# Diagnostic tables built from the per-pair quantities computed in C++.
# A port of python/spicyglm/diagnostics.py, which in turn mirrors spicyR's
# assembleDiagnosticsTables(), computePairSummary(), computeRelativeDiagnostics()
# and crossPairDiagnostics().

PATIENT_COLUMNS <- c("patient_id", "group", "n_i", "T_i", "S_g", "l_i", "raw_residual_sum",
                     "adjusted_residual_sum", "e_i", "influence_i", "y_i", "d_i", "delta_i")
IMAGE_COLUMNS <- c("patient_id", "image_id", "group", "n_ij", "density_ij", "l_ij",
                   "l_ij_group_share", "raw_residual_sum_ij", "adjusted_residual_sum_ij",
                   "e_ij", "e_ij_share_within_patient", "influence_ij")

# Pair summary row plus patient and image tables for one fitted pair.
pair_tables <- function(fit, f, t, levels_, cluster_labels, image_labels) {
  p <- fit$patient
  patient <- as.data.frame(p[setdiff(PATIENT_COLUMNS[-(1:2)], "S_g")], stringsAsFactors = FALSE)
  patient$patient_id <- as.character(cluster_labels[p$cluster_id + 1L])
  patient$group <- levels_[p$group + 1L]
  patient$S_g <- p$S_g[p$group + 1L]
  patient <- patient[order(patient$patient_id, method = "radix"), PATIENT_COLUMNS, drop = FALSE]
  rownames(patient) <- NULL

  im <- fit$image
  image <- as.data.frame(im[IMAGE_COLUMNS[-(1:3)]], stringsAsFactors = FALSE)
  image$patient_id <- as.character(cluster_labels[p$cluster_id[im$cluster + 1L] + 1L])
  image$image_id <- as.character(image_labels[im$image_id + 1L])
  image$group <- levels_[p$group[im$cluster + 1L] + 1L]
  image <- image[order(image$patient_id, image$image_id, method = "radix"),
                 IMAGE_COLUMNS, drop = FALSE]
  rownames(image) <- NULL

  abs_delta <- abs(patient$delta_i)
  i_infl <- which.max(patient$influence_i)
  i_delta <- if (any(!is.na(abs_delta))) which.max(abs_delta) else integer(0)
  summary <- data.frame(
    from = f, to = t, n_patients = nrow(patient), nu = fit$df,
    max_influence = max(patient$influence_i, na.rm = TRUE),
    patient_with_max_influence = patient$patient_id[i_infl],
    max_abs_delta_logRR = if (length(i_delta)) abs_delta[i_delta] else NA_real_,
    patient_with_max_delta_logRR = if (length(i_delta)) patient$patient_id[i_delta] else NA_character_,
    stringsAsFactors = FALSE)

  patient$from <- f; patient$to <- t
  image$from <- f;   image$to <- t
  list(summary = summary, patient = patient, image = image)
}

# dplyr::percent_rank: (min rank - 1) / (non-missing count - 1)
percent_rank <- function(s) {
  n <- sum(!is.na(s))
  (rank(s, ties.method = "min", na.last = "keep") - 1) / (n - 1)
}

# transform within groups defined by the columns `by`, preserving row order
group_transform <- function(df, by, fun) {
  key <- do.call(paste, c(unname(df[by]), sep = "\r"))
  out <- rep(NA_real_, nrow(df))
  for (idx in split(seq_len(nrow(df)), factor(key, levels = unique(key)))) out[idx] <- fun(idx)
  out
}

relative_diagnostics <- function(table, level) {
  out <- table
  pair <- c("from", "to")
  if (level == "patient") {
    out$n_patients_per_group <- group_transform(out, c(pair, "group"), length)
    out$rel_l_i <- out$l_i * out$n_patients_per_group
    out$percentile_rank_leverage <- group_transform(out, c(pair, "group"),
                                                    function(i) percent_rank(out$rel_l_i[i]))
    out$n_patients_in_pair <- group_transform(out, pair, length)
    out$rel_influence_i <- out$influence_i * out$n_patients_in_pair
    out$percentile_rank_influence <- group_transform(out, pair,
                                                     function(i) percent_rank(out$rel_influence_i[i]))
    out$abs_delta <- abs(out$delta_i)
    out$percentile_rank_delta <- group_transform(out, pair,
                                                 function(i) percent_rank(out$abs_delta[i]))
    return(out[, c("patient_id", "group", "from", "to", "n_patients_per_group", "n_patients_in_pair",
                   "l_i", "influence_i", "delta_i", "rel_l_i", "percentile_rank_leverage",
                   "rel_influence_i", "percentile_rank_influence", "percentile_rank_delta",
                   "n_i", "T_i", "S_g", "raw_residual_sum", "adjusted_residual_sum", "e_i",
                   "y_i", "d_i"), drop = FALSE])
  }
  out$n_images_for_patient <- group_transform(out, c(pair, "patient_id"), length)
  out$rel_l_ij <- out$l_ij * out$n_images_for_patient
  out$percentile_rank_leverage <- group_transform(out, c(pair, "patient_id"),
                                                  function(i) percent_rank(out$rel_l_ij[i]))
  out$n_images_in_pair <- group_transform(out, pair, length)
  out$rel_influence_ij <- out$influence_ij * out$n_images_in_pair
  out$percentile_rank_influence <- group_transform(out, pair,
                                                   function(i) percent_rank(out$rel_influence_ij[i]))
  out[, c("patient_id", "image_id", "group", "from", "to", "n_images_for_patient",
          "n_images_in_pair", "l_ij", "l_ij_group_share", "influence_ij", "rel_l_ij",
          "percentile_rank_leverage", "rel_influence_ij", "percentile_rank_influence",
          "density_ij", "raw_residual_sum_ij", "adjusted_residual_sum_ij", "e_ij",
          "e_ij_share_within_patient"), drop = FALSE]
}

# Wilson score interval, as binom::binom.wilson (no clipping)
wilson_interval <- function(x, n, conf_level = 0.95) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  z2 <- z * z
  p <- x / n
  centre <- p + 0.5 * z2 / n
  half <- z * sqrt((p * (1 - p) + 0.25 * z2 / n) / n)
  list(lower = (centre - half) / (1 + z2 / n), upper = (centre + half) / (1 + z2 / n))
}

cross_pair_diagnostics <- function(relative, top_percent = 0.05) {
  image_level <- "image_id" %in% names(relative)
  keys <- if (image_level) c("patient_id", "image_id") else "patient_id"
  key <- do.call(paste, c(unname(relative[keys]), sep = "\r"))
  if (any(tapply(relative$group, key, function(g) length(unique(g))) > 1))
    stop("a patient/image's condition group varies across cell-type pairs")

  rel <- relative
  pathways <- c("influence", "leverage", if (!image_level) "delta")
  for (pw in pathways)
    rel[[paste0("flagged_", pw)]] <- rel[[paste0("percentile_rank_", pw)]] >= 1 - top_percent

  idx <- split(seq_len(nrow(rel)), key)
  idx <- idx[order(names(idx), method = "radix")]        # groupby(sort=True)
  first <- vapply(idx, function(i) i[1], integer(1))
  out <- rel[first, keys, drop = FALSE]
  out$group <- rel$group[first]
  out$n_pairs_present <- vapply(idx, length, integer(1))
  if (image_level) {
    for (nc in list(c("mean_influence_share_ij", "e_ij_share_within_patient"),
                    c("mean_influence_ij", "influence_ij"),
                    c("mean_l_ij", "l_ij"),
                    c("mean_l_ij_group_share", "l_ij_group_share")))
      out[[nc[1]]] <- vapply(idx, function(i) mean(rel[[nc[2]]][i]), numeric(1))
  }
  for (pw in pathways)
    out[[paste0("n_pairs_flagged_", pw)]] <-
      vapply(idx, function(i) sum(rel[[paste0("flagged_", pw)]][i]), numeric(1))
  for (pw in pathways)
    out[[paste0("prop_flagged_", pw)]] <-
      out[[paste0("n_pairs_flagged_", pw)]] / out$n_pairs_present
  for (pw in pathways) {
    w <- wilson_interval(out[[paste0("n_pairs_flagged_", pw)]], out$n_pairs_present)
    out[[paste0("wilson_lower_", pw)]] <- w$lower
    out[[paste0("wilson_upper_", pw)]] <- w$upper
  }
  if (image_level)
    out <- out[, c("patient_id", "image_id", "group", "n_pairs_present", "mean_influence_share_ij",
                   "mean_influence_ij", "mean_l_ij", "mean_l_ij_group_share",
                   "n_pairs_flagged_influence", "prop_flagged_influence",
                   "wilson_lower_influence", "wilson_upper_influence",
                   "n_pairs_flagged_leverage", "prop_flagged_leverage",
                   "wilson_lower_leverage", "wilson_upper_leverage"), drop = FALSE]
  out <- out[order(-out$wilson_lower_influence, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}

assemble_diagnostics <- function(pair_diag, top_percent) {
  patient <- relative_diagnostics(
    do.call(rbind, lapply(pair_diag, `[[`, "patient")), "patient")
  image <- relative_diagnostics(
    do.call(rbind, lapply(pair_diag, `[[`, "image")), "image")
  list(pair = do.call(rbind, lapply(pair_diag, `[[`, "summary")),
       patient = patient, image = image,
       cross_pair = list(patient = cross_pair_diagnostics(patient, top_percent),
                         image = cross_pair_diagnostics(image, top_percent)))
}
