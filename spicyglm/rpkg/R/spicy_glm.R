# R front end for the spicyGLM core: pair enumeration, skip rules and result
# tables. A port of python/spicyglm/api.py; the numeric work happens in C++.

EFFECT_COLUMNS <- list(poisson = c("log_rate_ratio", "rate_ratio"),
                       binomial = c("log_odds_ratio", "odds_ratio"))
SKIP_COLUMNS <- c("from", "to", "reason", "message")

result_columns <- function(family) {
  c("from", "to", "condition_ref", "condition_comp", "coef_ref", "coef_comp",
    EFFECT_COLUMNS[[family]], "p_value", "estimator", "family",
    "mle_would_skip", "mle_skip_reason", "p_adj")
}

#' Test for a change in co-localisation of cell-type pairs between two conditions
#'
#' For \code{family = "poisson"} the number of target cells within radius \code{r}
#' of each reference cell is Poisson with a log-density offset; for
#' \code{family = "binomial"} the number of target cells among its \code{k}
#' nearest neighbours is Binomial with a logit background-proportion offset.
#' Either way there is one coefficient per condition, and the log rate (odds)
#' ratio is tested with the CR2 cluster-robust variance, clustered by
#' \code{subject} or by image when \code{subject} is missing or one-to-one with
#' images, and Satterthwaite degrees of freedom; \code{cr2_method = "naive"}
#' uses the model-based variance and a z-test instead.
#'
#' @param cells A data frame with one row per cell.
#' @param condition Column holding the image-level condition. Exactly two levels
#'   are required; the first level (factor order, otherwise sorted) is the reference.
#' @param r Neighbourhood radius (poisson only).
#' @param subject Optional column holding the clustering unit.
#' @param image_id,cell_type Column names.
#' @param spatial_coords Length-two character vector naming the coordinate columns.
#' @param from,to Cell types. A single \code{from} and \code{to} fits that pair
#'   only. Otherwise, for \code{family = "poisson"} (direction-invariant) every
#'   pair among the given (or all) cell types is fitted, one direction per
#'   unordered pair plus self-pairs; for \code{family = "binomial"} (directional:
#'   from -> to and to -> from differ) every ordered pair in \code{from} x
#'   \code{to} is fitted, where an omitted side means all cell types.
#' @param window Observation window for each image's area, "convex" or "rectangle".
#' @param family "poisson" or "binomial".
#' @param k Number of nearest neighbours (binomial only).
#' @param estimator "firth" or "mle".
#' @param cr2_method "fast" or "naive".
#' @param compute_diagnostics Add per-patient and per-image leverage, influence and
#'   point-estimate shift, their within-pair percentile ranks, and cross-pair
#'   flagging. Requires \code{family = "poisson"}, \code{estimator = "firth"}
#'   and \code{cr2_method = "fast"}.
#' @param top_percent Fraction of the within-pair percentile rank counted as flagged.
#' @param n_threads Threads used to build the neighbour lists.
#'
#' @return A list with \code{results} (one row per fitted pair),
#'   \code{skipped} (pairs that could not be fitted, with a reason code) and,
#'   when \code{compute_diagnostics} is set, \code{diagnostics} with elements
#'   "pair", "patient", "image" and "cross_pair".
#' @export
spicy_glm <- function(cells, condition, r = NULL, subject = NULL,
                      image_id = "imageID", cell_type = "cellType",
                      spatial_coords = c("x", "y"), from = NULL, to = NULL,
                      window = "convex", family = c("poisson", "binomial"), k = NULL,
                      estimator = c("firth", "mle"), cr2_method = c("fast", "naive"),
                      compute_diagnostics = FALSE, top_percent = 0.05, n_threads = 1) {
  family <- match.arg(family)
  estimator <- match.arg(estimator)
  cr2_method <- match.arg(cr2_method)
  if (!is.data.frame(cells)) stop("`cells` must be a data frame with one row per cell.")
  if (length(n_threads) != 1 || n_threads < 1 || n_threads != as.integer(n_threads))
    stop("`n_threads` must be a positive integer.")
  if (family == "poisson") {
    if (is.null(r) || !is.finite(r) || r <= 0)
      stop("family = 'poisson' requires a positive radius `r`.")
  } else {
    if (is.null(k) || length(k) != 1 || !is.finite(k) || k < 1 || k != as.integer(k))
      stop("family = 'binomial' requires `k`, a single positive integer.")
    k <- as.integer(k)
    if (!is.null(r)) message("family = 'binomial' uses `k` nearest neighbours; `r` is ignored.")
  }
  if (compute_diagnostics && (family != "poisson" || estimator != "firth" || cr2_method != "fast")) {
    warning("compute_diagnostics requires family = 'poisson', estimator = 'firth' and ",
            "cr2_method = 'fast'; diagnostics are not computed.", call. = FALSE)
    compute_diagnostics <- FALSE
  }
  x_col <- spatial_coords[1]; y_col <- spatial_coords[2]
  needed <- c(condition, image_id, cell_type, x_col, y_col, if (!is.null(subject)) subject)
  missing_cols <- setdiff(needed, names(cells))
  if (length(missing_cols))
    stop("columns not found in `cells`: ", paste(missing_cols, collapse = ", "))

  # sort cells by image so each image is a contiguous block for the core
  image_chr <- as.character(cells[[image_id]])
  image_labels <- sort(unique(image_chr))
  image_codes <- match(image_chr, image_labels) - 1L
  ord <- order(image_codes, method = "radix")
  df <- cells[ord, , drop = FALSE]
  image_codes <- image_codes[ord]
  n_images <- length(image_labels)

  levels_ <- condition_levels(cells[[condition]])
  if (length(levels_) != 2L)
    stop("spicyGLM compares exactly two conditions; found ", length(levels_), ": ",
         paste(levels_, collapse = ", "))
  cond_chr <- as.character(df[[condition]])
  per_image_cond <- split(cond_chr, image_codes)
  if (any(vapply(per_image_cond, function(z) length(unique(z)) > 1L, logical(1))))
    stop("'", condition, "' must be constant within each image.")
  if (!is.null(subject)) {
    bysub <- split(as.character(cells[[condition]]), as.character(cells[[subject]]))
    if (any(vapply(bysub, function(z) length(unique(z)) > 1L, logical(1))))
      stop("each subject must belong to a single '", condition, "' level.")
  }
  image_condition <- vapply(per_image_cond[as.character(seq_len(n_images) - 1L)],
                            function(z) z[1], character(1))
  image_group <- match(image_condition, levels_) - 1L

  one_to_one <- is.null(subject) ||
    length(unique(as.character(cells[[subject]]))) == length(unique(image_chr))
  if (one_to_one) {
    image_cluster <- seq_len(n_images) - 1L
    cluster_labels <- image_labels
  } else {
    sub_chr <- as.character(df[[subject]])
    per_image_sub <- vapply(split(sub_chr, image_codes)[as.character(seq_len(n_images) - 1L)],
                            function(z) z[1], character(1))
    cluster_labels <- unique(per_image_sub)            # first-appearance order, as pd.factorize
    image_cluster <- match(per_image_sub, cluster_labels) - 1L
  }

  type_labels <- unique(as.character(cells[[cell_type]]))   # first-appearance order, as R's unique()
  type_codes <- as.integer(match(as.character(df[[cell_type]]), type_labels) - 1L)
  image_offsets <- as.integer(c(0L, cumsum(tabulate(image_codes + 1L, nbins = n_images))))

  data <- dataset_create(as.numeric(df[[x_col]]), as.numeric(df[[y_col]]),
                         type_codes, image_offsets, length(type_labels))
  presence <- presence_table(type_codes, image_codes, image_group, length(type_labels))
  areas <- NULL
  if (family == "poisson") {
    areas <- dataset_image_areas(data, window)
    dataset_build_radius_index(data, r)
  } else {
    dataset_build_knn(data, k, as.integer(n_threads))
  }

  pairs <- enumerate_pairs(from, to, type_labels, family)
  unknown <- setdiff(unique(unlist(pairs)), type_labels)
  if (length(unknown)) stop("cell type not found: ", paste(unknown, collapse = ", "))

  ctx <- list(family = family, k = k, estimator = estimator, cr2_method = cr2_method,
              levels = levels_, image_group = image_group, image_cluster = image_cluster,
              type_index = stats::setNames(seq_along(type_labels) - 1L, type_labels),
              data = data, areas = areas, presence = presence)

  outcomes <- lapply(pairs, function(p) fit_one_pair(ctx, p[1], p[2], compute_diagnostics))
  is_skip <- vapply(outcomes, function(o) !is.null(o$reason), logical(1))

  pair_diag <- NULL
  if (compute_diagnostics && any(!is_skip))
    pair_diag <- lapply(outcomes[!is_skip], function(o)
      pair_tables(o$.fit, o$from, o$to, levels_, cluster_labels, image_labels))
  for (i in which(!is_skip)) outcomes[[i]]$.fit <- NULL

  skipped <- if (any(is_skip)) {
    do.call(rbind, lapply(outcomes[is_skip], function(o)
      data.frame(from = o$from, to = o$to, reason = o$reason, message = o$message,
                 stringsAsFactors = FALSE)))
  } else {
    stats::setNames(data.frame(matrix(character(0), 0, 4), stringsAsFactors = FALSE), SKIP_COLUMNS)
  }

  cols <- result_columns(family)
  results <- if (any(!is_skip)) {
    do.call(rbind, lapply(outcomes[!is_skip], function(o)
      as.data.frame(o[cols[-length(cols)]], stringsAsFactors = FALSE)))
  } else {
    stats::setNames(data.frame(matrix(numeric(0), 0, length(cols) - 1)), cols[-length(cols)])
  }
  if (nrow(results)) {
    results$p_adj <- stats::p.adjust(results$p_value, method = "BH")
    results <- results[order(results$p_adj, method = "radix"), , drop = FALSE]
    rownames(results) <- NULL
  } else {
    results$p_adj <- numeric(0)
  }
  list(results = results, skipped = skipped,
       diagnostics = if (!is.null(pair_diag)) assemble_diagnostics(pair_diag, top_percent))
}

condition_levels <- function(col) {
  if (is.factor(col)) {
    present <- unique(as.character(col[!is.na(col)]))
    return(levels(col)[levels(col) %in% present])
  }
  sort(unique(as.character(col[!is.na(col)])))
}

enumerate_pairs <- function(from, to, all_types, family) {
  if (is.character(from) && length(from) == 1L && is.character(to) && length(to) == 1L)
    return(list(c(from, to)))
  # the kNN effect is directional (A->B != B->A), so fit every ordered pair
  if (family == "binomial") {
    from_types <- if (is.null(from)) all_types else unique(from)
    to_types <- if (is.null(to)) all_types else unique(to)
    return(unlist(lapply(from_types, function(f) lapply(to_types, function(t) c(f, t))),
                  recursive = FALSE))
  }
  types <- if (!is.null(from) || !is.null(to)) unique(c(from, to)) else all_types
  cross <- if (length(types) >= 2L)
    utils::combn(types, 2L, simplify = FALSE) else list()
  c(cross, lapply(types, function(t) c(t, t)))
}

# present[type, group] = number of images in the group containing the type
presence_table <- function(type_codes, image_codes, image_group, n_types) {
  has <- matrix(FALSE, n_types, length(image_group))
  has[cbind(type_codes + 1L, image_codes + 1L)] <- TRUE
  cbind(rowSums(has[, image_group == 0L, drop = FALSE]),
        rowSums(has[, image_group == 1L, drop = FALSE]))
}

skip_pair <- function(f, t, reason, message) {
  list(from = f, to = t, reason = reason, message = message)
}

fit_one_pair <- function(ctx, f, t, compute_diagnostics = FALSE) {
  from_code <- ctx$type_index[[f]]; to_code <- ctx$type_index[[t]]
  md <- if (ctx$family == "poisson")
    dataset_poisson_model_data(ctx$data, ctx$areas, from_code, to_code)
  else
    dataset_binomial_model_data(ctx$data, from_code, to_code)

  group <- ctx$image_group[md$image + 1L]
  present <- unique(group)
  if (length(present) < 2L) {
    d <- diagnose_missing(f, t, setdiff(0:1, present), ctx)
    return(skip_pair(f, t, d$reason, d$message))
  }

  n <- md$n
  b <- boundary_reason(f, t, n, group, ctx)
  if (ctx$estimator == "mle" && !is.null(b$reason))
    return(skip_pair(f, t, b$reason, b$message))

  cluster <- ctx$image_cluster[md$image + 1L]
  if (ctx$cr2_method == "fast") {
    for (g in 0:1) {
      if (length(unique(cluster[group == g])) < 2L)
        return(skip_pair(f, t, "one_patient_per_group", paste0(
          "Skipping pair ", f, "__", t, ": condition '", ctx$levels[g + 1L],
          "' has fewer than two clusters with data for this pair; CR2 needs at least two per group.")))
    }
  }

  fit <- if (ctx$family == "poisson")
    fit_pair_poisson_cpp(cluster, md$image, group, n, md$density, ctx$estimator, ctx$cr2_method,
                         compute_diagnostics)
  else
    fit_pair_binomial_cpp(cluster, md$image, group, n, ctx$k, md$p0, ctx$estimator, ctx$cr2_method)

  b_ref <- fit$beta[1]; b_comp <- fit$beta[2]
  log_effect <- b_comp - b_ref
  t_stat <- log_effect / sqrt(fit$v_hat)
  p_value <- if (ctx$cr2_method == "naive")
    2 * stats::pnorm(-abs(t_stat)) else 2 * stats::pt(-abs(t_stat), fit$df)

  eff <- EFFECT_COLUMNS[[ctx$family]]
  out <- list(from = f, to = t, condition_ref = ctx$levels[1], condition_comp = ctx$levels[2],
              coef_ref = b_ref, coef_comp = b_comp,
              log_effect = log_effect, effect = exp(log_effect),
              p_value = p_value, estimator = ctx$estimator, family = ctx$family,
              mle_would_skip = !is.null(b$reason),
              mle_skip_reason = if (is.null(b$reason)) NA_character_ else b$reason)
  if (compute_diagnostics) out$.fit <- fit
  names(out)[names(out) == "log_effect"] <- eff[1]
  names(out)[names(out) == "effect"] <- eff[2]
  out
}

# Conditions where the MLE is infinite: all counts zero, or (binomial) all k.
boundary_reason <- function(f, t, n, group, ctx) {
  floor_g <- ctx$levels[vapply(0:1, function(g) !any(n[group == g] > 0), logical(1))]
  ceil_g <- if (ctx$family == "binomial")
    ctx$levels[vapply(0:1, function(g) all(n[group == g] == ctx$k), logical(1))] else character(0)
  boundary <- c(floor_g, setdiff(ceil_g, floor_g))
  if (!length(boundary)) return(list(reason = NULL, message = NULL))
  reason <- if (length(boundary) == 2L) {
    if (!length(ceil_g)) "all_zero" else if (!length(floor_g)) "all_max" else "all_boundary"
  } else if (length(floor_g)) "one_condition_zero" else "one_condition_max"
  list(reason = reason, message = paste0(
    "Skipping pair ", f, "__", t, ": neighbour counts sit at a separation boundary in condition(s) ",
    paste(boundary, collapse = ", "), " (reason '", reason,
    "'), so the MLE is not finite. Use estimator = 'firth'."))
}

diagnose_missing <- function(f, t, missing_groups, ctx) {
  reasons <- character(0); texts <- character(0)
  for (g in missing_groups) {
    level <- ctx$levels[g + 1L]
    f_in <- ctx$presence[ctx$type_index[[f]] + 1L, g + 1L] > 0
    t_in <- ctx$presence[ctx$type_index[[t]] + 1L, g + 1L] > 0
    if (!f_in && !t_in) {
      reasons <- c(reasons, "both_absent")
      texts <- c(texts, paste0("condition '", level, "' has no images containing '", f, "' or '", t, "'."))
    } else if (!(f_in && t_in)) {
      reasons <- c(reasons, "type_absent")
      texts <- c(texts, paste0("condition '", level, "' has no images containing '",
                               if (f_in) t else f, "'."))
    } else {
      reasons <- c(reasons, "no_cooccurrence")
      texts <- c(texts, paste0("condition '", level, "' has no image containing both '", f,
                               "' and '", t, "'."))
    }
  }
  list(reason = if (length(unique(reasons)) == 1L) reasons[1] else "mixed",
       message = paste0("Skipping pair ", f, "__", t, ": ", paste(texts, collapse = " ")))
}
