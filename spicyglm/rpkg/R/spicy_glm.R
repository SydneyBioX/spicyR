# R front end for the spicyGLM core: pair enumeration, skip rules and result
# tables. A port of python/spicyglm/api.py; the numeric work happens in C++.

EFFECT_COLUMNS <- list(poisson = c("log_rate_ratio", "rate_ratio"),
                       binomial = c("log_odds_ratio", "odds_ratio"))
SKIP_COLUMNS <- c("from", "to", "reason", "message")

result_columns <- function(family, effect = "ratio") {
  c("from", "to", "condition_ref", "condition_comp", "coef_ref", "coef_comp",
    if (effect == "excess") "excess_difference" else EFFECT_COLUMNS[[family]], "p_value", "estimator", "family",
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
#'   are required; the reference is \code{ref} if given, otherwise the first level
#'   (factor order, otherwise sorted).
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
#' @param window Observation window of each image, "convex" (the cells' convex
#'   hull) or "rectangle" (their bounding box): its area sets the offset, and in
#'   the inhomogeneous design it also bounds the intensity discs and the
#'   translation edge correction.
#' @param family "poisson" or "binomial".
#' @param k Number of nearest neighbours (binomial only).
#' @param sigma Radius of the disc kernel for the inhomogeneous cross-K design
#'   (poisson only), in the units of \code{spatial_coords}. \code{NULL} (default)
#'   fits the homogeneous design, whose null is complete spatial randomness over
#'   each image. When supplied, each cell type's intensity at each of its cells
#'   is the number of other cells of the type within \code{sigma} over the area
#'   of that disc inside the window, and every reference--target pair within
#'   \code{r} is weighted by the renormalised inverse intensities of its two
#'   cells times a translation edge correction, so the group coefficients are
#'   the pooled inhomogeneous cross-K function over \code{pi * r^2} and the
#'   effect is the log ratio of intensity-adjusted co-localisation. Structure at
#'   scales up to about \code{sigma} is treated as trend rather than
#'   interaction, so choose \code{sigma} several times \code{r}. As spicyR's
#'   \code{spicyGLM(sigma =)}.
#' @param min_lambda Floor for the estimated intensities, as a fraction of each
#'   image's average intensity for the type (inhomogeneous design only).
#' @param edge_correct Apply the edge correction: the translation correction
#'   to each neighbour pair in the inhomogeneous design, or in the Kontextual
#'   Poisson design the correction of each cell's context count for the part of
#'   its disc outside the window (as Statial's Kontextual). No effect otherwise.
#' @param parent Cell types forming the context (parent) population for
#'   Kontextual (Ameen et al.; Statial's \code{Kontextual()}): the co-localisation
#'   of \code{from} with \code{to} is judged against where the context lies, so
#'   the null is that \code{to} cells are a random subset of the context cells.
#'   The \code{to} type must belong to \code{parent}. For \code{family =
#'   "poisson"}, each reference cell's \code{to} neighbours within \code{r} are
#'   weighted by its context intensity over theirs, against an offset of the
#'   \code{to} share of the context times its context cells within \code{r}; a
#'   reference cell's weight is its context intensity, so cells with no context
#'   cell within \code{r} drop out. For \code{family = "binomial"}, the context
#'   cells among each reference cell's \code{k} nearest neighbours are the
#'   trials and the \code{to} cells among them the successes, against the
#'   \code{to} share of the context. Pairs are directional: every ordered pair
#'   in \code{from} x \code{to} with \code{to} in \code{parent} is fitted, an
#'   omitted \code{to} meaning all of \code{parent}. Self-pairs leave the cell
#'   itself out of the counts. \code{NULL} (default) is the ordinary design.
#' @param null For family = "poisson": the null the offset expresses. "csr"
#'   (default): complete spatial randomness in the window, n_to / area * pi r^2.
#'   "random_labelling": condition on every cell's location and on the REF
#'   cells; each REF cell's offset is its number of candidate neighbours
#'   within r (the non-REF cells, or all other cells for a self-pair) times the
#'   TARGET share of the candidates. This removes the image-wide shift in
#'   co-localisation that arises when tissue fills a window unevenly. Pairs are
#'   then directional, as for the binomial design.
#' @param frailty NULL (default) means TRUE for \code{effect = "excess"} and FALSE otherwise.
#'   TRUE fits the cell-level GLM with an image (subject) frailty and
#'   the random-labelling within-image variance: each image's cells get the
#'   prior weight a_u / phi_i, where phi_i is the exact random-labelling
#'   variance of the image's count over its Poisson (binomial) variance and
#'   a_u = 1 / (1 + tau2 J_u) interpolates between count weighting and equal
#'   weighting of subjects. tau2 is estimated per pair (Paule-Mandel), and the
#'   test uses the CR2 variance under this working covariance with
#'   Satterthwaite df. Available for every design: for Kontextual (`parent`)
#'   phi is the random-labelling variance within the context; for the
#'   inhomogeneous design (`sigma`) it is the variance under an inhomogeneous
#'   Poisson null for the TARGET (Campbell's theorem). Per-pair tau2, df and SE
#'   are returned in `$frailty`.
#' @param label_clustering With frailty = TRUE: inflate each image's phi for
#'   TARGET cells that cluster among themselves, which random labelling ignores.
#'   A spatial HAC (Conley) estimate of Var(O) (Bartlett weights, bandwidth 2r,
#'   residuals of the TARGET indicator within the image) is divided by the
#'   random-labelling variance, and the ratio is pooled to one factor per
#'   (image, TARGET type), the median over the fitted REF types, floored at 1. Not
#'   used for the inhomogeneous design, whose phi comes from a Poisson-process null.
#' @param moderate With frailty = TRUE: shrink each pair's CR2 variance toward
#'   its frailty-model variance at a common tau2 (the median of the per-pair
#'   estimates) by empirical Bayes (Smyth's
#'   hyperparameter estimates, with no trend: the model variance already
#'   carries the dependence on counts). The test uses the moderated variance
#'   on d + d0 degrees of freedom. Needs several pairs. `$frailty` then holds
#'   the prior df d0 and scale s0 as attributes.
#' @param effect "excess" (default): the additive effect under random labelling. For
#'   a directional pair REF -> TARGET, each group coefficient is the mean number
#'   of extra REF cells among a TARGET cell's neighbours (within \code{r}, or among
#'   its \code{k} nearest) beyond random labelling of the non-REF cells; for the
#'   radius graph that is lambda_REF (K - K under random labelling), the
#'   KAMP-adjusted cross-K in neighbour units. Recruiting a fraction f of the
#'   TARGET cells next to REF cells gives f whatever the density or abundance of
#'   either type, whereas a ratio is capped by composition (the observed/expected
#'   ratio cannot exceed 1 / the share of candidates near a REF cell) and so
#'   shrinks in dense images. Each image's count has its exact random-labelling
#'   mean and variance; the fit is the frailty GEE (frailty = TRUE; tau2 by
#'   Paule-Mandel, on the scale of neighbours per cell squared) or working
#'   independence (frailty = FALSE), with the closed-form CR2 variance and
#'   Satterthwaite df. \code{label_clustering} applies. Returns
#'   \code{excess_difference} (comparison minus reference) and per-pair tau2, df
#'   and SE in \code{$frailty}. \code{moderate = TRUE} shrinks each pair's CR2
#'   variance toward its working-model variance as for the ratio (one prior df, no
#'   trend), which matters when a few high-information units dominate and the CR2
#'   df are small. Not available with \code{sigma}, \code{parent} or
#'   \code{test = "moderated"} (when \code{effect} is not given, these use "ratio"); \code{null}
#'   and \code{estimator} are not used. "ratio": the effect is a log rate (odds) ratio, as described above.
#' @param variance For \code{effect = "excess"}: the variance of the difference between
#'   conditions. "cr2" (default): the CR2
#'   variance on Satterthwaite df, calibrated with few units in simulations and with shuffled labels on real
#'   data. "hartung_knapp": the working-model variance at the pair's tau2,
#'   scaled up by the Pearson dispersion when that exceeds 1 (modified Hartung-Knapp) and never
#'   below the CR2 variance, tested on n_units - 2 df. This is exact for normal unit summaries
#'   whose variances are known up to a common scale, and it keeps inference informative when a
#'   condition has only a few patients, where the CR2 Satterthwaite df collapse, but it can be
#'   anti-conservative for rare types. With
#'   \code{moderate = TRUE} the CR2 variance is moderated and this argument is ignored.
#' @param test "glm" (default): the cell-level GLM with CR2 variance described
#'   above. "moderated": an image-level test. Each image is summarised by its
#'   log observed / expected ratio (Poisson designs) or log odds ratio against
#'   its background (binomial designs), images are averaged within
#'   \code{subject}, and the two conditions are compared by a two-sample
#'   t-test whose variance is moderated across all fitted pairs by empirical
#'   Bayes with a trend on the pair's mean log count (limma's \code{eBayes(trend
#'   = TRUE)}). It weights units equally rather than by cell counts, so one
#'   large or unusual image cannot dominate, and the moderation borrows degrees
#'   of freedom across pairs. Needs several pairs; diagnostics are not
#'   available.
#' @param density_adjust For \code{test = "moderated"}: remove the dependence of
#'   the image summaries on the two cell types' abundances (log counts per
#'   image). The two slopes are shared by all fitted pairs and estimated
#'   within each pair and condition (an analysis of covariance with a common
#'   slope), so that a difference in abundance between conditions is not
#'   reported as a change in co-localisation.
#' @param estimator "firth" or "mle".
#' @param cr2_method "fast" or "naive".
#' @param compute_diagnostics Add per-patient and per-image leverage, influence and
#'   point-estimate shift, their within-pair percentile ranks, and cross-pair
#'   flagging. Requires \code{family = "poisson"}, \code{estimator = "firth"}
#'   and \code{cr2_method = "fast"}.
#' @param top_percent Fraction of the within-pair percentile rank counted as flagged.
#' @param n_threads Threads used to build the neighbour lists.
#' @param ref Reference level of \code{condition}, as spicyR's \code{spicyGLM(ref =)}.
#'   Must be one of the levels present; the other level is the comparison.
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
                      sigma = NULL, min_lambda = 0.05, edge_correct = TRUE, parent = NULL,
                      null = c("csr", "random_labelling"), frailty = NULL, moderate = FALSE,
                      label_clustering = TRUE, effect = c("excess", "ratio"), variance = c("cr2", "hartung_knapp"),
                      test = c("glm", "moderated"), density_adjust = FALSE,
                      estimator = c("firth", "mle"), cr2_method = c("fast", "naive"),
                      compute_diagnostics = FALSE, top_percent = 0.05, n_threads = 1,
                      ref = NULL) {
  family <- match.arg(family)
  null <- match.arg(null)
  test <- match.arg(test)
  # excess is the default; the designs that exist only for the ratio fall back to it when effect is not given
  if (missing(effect) && (!is.null(sigma) || !is.null(parent) || test == "moderated")) effect <- "ratio"
  effect <- match.arg(effect)
  variance <- match.arg(variance)
  if (effect == "excess") {
    if (!is.null(sigma) || !is.null(parent))
      stop("effect = 'excess' uses the random-labelling null; it is not available with `sigma` or `parent`.")
    if (test == "moderated")
      stop("effect = 'excess' is a cell-level model; it is not available with test = 'moderated'.")
  }
  # frailty defaults to TRUE for effect = "excess" (patient heterogeneity in the working covariance)
  # and to FALSE for the rate-ratio GLM, as before
  if (is.null(frailty)) frailty <- identical(match.arg(effect), "excess")
  if (!is.logical(frailty) || length(frailty) != 1 || is.na(frailty))
    stop("`frailty` must be TRUE or FALSE.")
  if (null == "random_labelling" && family != "poisson")
    stop("null = 'random_labelling' applies to family = 'poisson'; the k-nearest-neighbour ",
         "design already conditions on each cell's neighbourhood.")
  if (null == "random_labelling" && (!is.null(sigma) || !is.null(parent)))
    stop("null = 'random_labelling' replaces `sigma` and `parent`; supply only one adjustment.")
  if (!is.logical(moderate) || length(moderate) != 1 || is.na(moderate))
    stop("`moderate` must be TRUE or FALSE.")
  if (moderate && !frailty && effect[1] != "excess") stop("moderate = TRUE needs frailty = TRUE.")
  if (!is.logical(label_clustering) || length(label_clustering) != 1 || is.na(label_clustering))
    stop("`label_clustering` must be TRUE or FALSE.")
  if (frailty && test == "moderated")
    stop("frailty = TRUE is a cell-level model; use it with test = 'glm'.")
  if (frailty && effect == "ratio" && (estimator[1] != "firth" || cr2_method[1] == "naive"))
    stop("frailty = TRUE uses the Firth-adjusted estimate and the CR2 variance.")
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
  if (!is.null(sigma)) {
    if (family == "binomial") {
      message("`sigma` applies to family = 'poisson' only (the k-nearest-neighbour design ",
              "already adapts to local cell density); ignoring it.")
      sigma <- NULL
    } else {
      if (!is.numeric(sigma) || length(sigma) != 1 || !is.finite(sigma) || sigma <= 0)
        stop("`sigma` must be a single positive number: the radius of the disc kernel ",
             "used to estimate each cell type's intensity.")
      if (!is.numeric(min_lambda) || length(min_lambda) != 1 || !is.finite(min_lambda) ||
          min_lambda <= 0)
        stop("`min_lambda` must be a single positive number.")
      if (!is.logical(edge_correct) || length(edge_correct) != 1 || is.na(edge_correct))
        stop("`edge_correct` must be TRUE or FALSE.")
      if (sigma <= r)
        message("`sigma` (", sigma, ") is no larger than `r` (", r, "): intensity smoothing at ",
                "the scale of the test radius absorbs the co-localisation being tested. ",
                "Consider `sigma` several times `r`.")
    }
  }
  if (!is.null(parent)) {
    if (!is.character(parent) || !length(parent) || anyNA(parent))
      stop("`parent` must be a character vector of cell types (the Kontextual context).")
    if (!is.null(sigma))
      stop("`sigma` and `parent` are alternative adjustments; supply one of them.")
    parent <- unique(parent)
  }
  if (compute_diagnostics && (test == "moderated" || frailty || null != "csr" || effect == "excess")) {
    warning("compute_diagnostics is available for the default CSR GLM only; diagnostics are not computed.",
            call. = FALSE)
    compute_diagnostics <- FALSE
  }
  if (compute_diagnostics && test == "moderated") {
    warning("compute_diagnostics is not available for test = 'moderated'; diagnostics are not computed.",
            call. = FALSE)
    compute_diagnostics <- FALSE
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
  if (!is.null(ref)) {
    if (length(ref) != 1L || !as.character(ref) %in% levels_)
      stop("`ref` = \"", paste(ref, collapse = ", "), "\" is not a level of `condition`. Available levels: ",
           paste(levels_, collapse = ", "), call. = FALSE)
    levels_ <- c(as.character(ref), setdiff(levels_, as.character(ref)))
  }
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
  if (!is.null(parent)) {
    unknown <- setdiff(parent, type_labels)
    if (length(unknown)) stop("`parent` cell type not found: ", paste(unknown, collapse = ", "))
  }
  if (family == "poisson") {
    areas <- dataset_image_areas(data, window)
    dataset_build_radius_index(data, r)
    if (!is.null(sigma)) dataset_build_intensity(data, sigma, min_lambda, window)
  } else {
    dataset_build_knn(data, k, as.integer(n_threads))
  }
  if (!is.null(parent))
    dataset_build_context(data, as.integer(match(parent, type_labels) - 1L), window, edge_correct)

  pairs <- enumerate_pairs(from, to, type_labels, if (null == "random_labelling" || effect == "excess") "binomial" else family, parent)
  unknown <- setdiff(unique(unlist(pairs)), type_labels)
  if (length(unknown)) stop("cell type not found: ", paste(unknown, collapse = ", "))
  if (!is.null(parent)) {
    outside <- setdiff(unique(vapply(pairs, `[`, character(1), 2L)), parent)
    if (length(outside))
      stop("Kontextual needs every `to` cell type inside `parent` (the context); not in it: ",
           paste(outside, collapse = ", "))
  }

  frailty_data <- NULL
  if (effect == "excess") {
    knn <- family == "binomial"
    frailty_data <- list(totals = dataset_pair_neighbour_totals(data, knn, length(type_labels)),
                         out_sq_totals = dataset_pair_neighbour_out_sq_totals(data, knn, length(type_labels)),
                         counts = unclass(table(factor(image_codes, levels = seq_len(n_images) - 1L),
                                                factor(type_codes, levels = seq_along(type_labels) - 1L))))
  } else if (frailty && is.null(sigma) && is.null(parent)) {
    knn <- family == "binomial"
    frailty_data <- list(totals = dataset_pair_neighbour_totals(data, knn, length(type_labels)),
                         sq_totals = dataset_pair_neighbour_sq_totals(data, knn, length(type_labels)),
                         counts = unclass(table(factor(image_codes, levels = seq_len(n_images) - 1L),
                                                factor(type_codes, levels = seq_along(type_labels) - 1L))))
  }
  if (frailty && is.null(frailty_data))
    frailty_data <- list(counts = unclass(table(factor(image_codes, levels = seq_len(n_images) - 1L),
                                                factor(if (is.null(parent)) type_codes else type_codes,
                                                       levels = seq_along(type_labels) - 1L))))
  ctx <- c(frailty_data, list(null = null, frailty = frailty, moderate = moderate, effect = effect, variance = variance))
  ctx <- c(ctx, list(family = family, r = r, k = k, sigma = sigma, edge_correct = edge_correct, parent = parent,
              estimator = estimator, cr2_method = cr2_method,
              levels = levels_, image_group = image_group, image_cluster = image_cluster,
              type_index = stats::setNames(seq_along(type_labels) - 1L, type_labels),
              data = data, areas = areas, presence = presence))

  if ((frailty || effect == "excess") && label_clustering && is.null(sigma)) {
    if (family == "binomial") {
      # the HAC bandwidth needs a distance: twice the typical k-NN radius
      n_img_cells <- tabulate(image_codes + 1L, nbins = n_images)
      scale <- stats::median(sqrt(k / (pi * n_img_cells / dataset_image_areas(data, window))))
      dataset_build_radius_index(data, scale)
      ctx$hac_h <- 2 * scale
    } else ctx$hac_h <- 2 * r
    if (is.null(ctx$totals)) {
      if (family == "poisson" && is.null(parent)) ctx$pair_totals <- dataset_pair_neighbour_totals(data, FALSE, length(type_labels))
      else if (family == "binomial") ctx$pair_totals <- dataset_pair_neighbour_totals(data, TRUE, length(type_labels))
    } else ctx$pair_totals <- ctx$totals
    ctx$phi_inflation <- frailty_label_clustering(ctx, pairs)
  }
  outcomes <- if (test == "moderated") {
    moderated_tests(ctx, pairs, density_adjust,
                    table(factor(image_codes, levels = seq_len(n_images) - 1L),
                          factor(type_codes, levels = seq_along(type_labels) - 1L)))
  } else if (effect == "excess") {
    lapply(pairs, function(p) excess_outcome(ctx, p[1], p[2]))
  } else {
    lapply(pairs, function(p) fit_one_pair(ctx, p[1], p[2], compute_diagnostics))
  }
  frailty_prior <- NULL
  if (moderate && effect == "excess") { outcomes <- excess_moderate(outcomes, ctx); frailty_prior <- attr(outcomes, "frailty_prior") }
  else if (frailty && moderate) { outcomes <- frailty_moderate(outcomes, ctx); frailty_prior <- attr(outcomes, "frailty_prior") }
  moderation <- attr(outcomes, "moderation")
  is_skip <- vapply(outcomes, function(o) !is.null(o$reason), logical(1))

  pair_diag <- NULL
  if (compute_diagnostics && any(!is_skip))
    pair_diag <- lapply(outcomes[!is_skip], function(o)
      pair_tables(o$.fit, o$from, o$to, levels_, cluster_labels, image_labels))
  frailty_table <- NULL
  if ((frailty || effect == "excess") && any(!is_skip))
    frailty_table <- do.call(rbind, lapply(outcomes[!is_skip], function(o)
      data.frame(from = o$from, to = o$to, tau2 = o$.frailty[["tau2"]], df = o$.frailty[["df"]],
                 se = o$.frailty[["se"]], stringsAsFactors = FALSE)))
  for (i in which(!is_skip)) { outcomes[[i]]$.fit <- NULL; outcomes[[i]]$.frailty <- NULL; outcomes[[i]]$.images <- NULL }

  skipped <- if (any(is_skip)) {
    do.call(rbind, lapply(outcomes[is_skip], function(o)
      data.frame(from = o$from, to = o$to, reason = o$reason, message = o$message,
                 stringsAsFactors = FALSE)))
  } else {
    stats::setNames(data.frame(matrix(character(0), 0, 4), stringsAsFactors = FALSE), SKIP_COLUMNS)
  }

  cols <- result_columns(family, effect)
  results <- if (any(!is_skip)) {
    do.call(rbind, lapply(outcomes[!is_skip], function(o)
      as.data.frame(o[cols[-length(cols)]], stringsAsFactors = FALSE)))
  } else {
    stats::setNames(data.frame(matrix(numeric(0), 0, length(cols) - 1)), cols[-length(cols)])
  }
  if (nrow(results)) {
    results$p_adj <- stats::p.adjust(results$p_value, method = "BH")
    results <- results[order(results$p_value, method = "radix"), , drop = FALSE]
    rownames(results) <- NULL
  } else {
    results$p_adj <- numeric(0)
  }
  out <- list(results = results, skipped = skipped,
              diagnostics = if (!is.null(pair_diag)) assemble_diagnostics(pair_diag, top_percent))
  if (!is.null(parent)) out$parent <- parent
  if (!is.null(moderation)) out$moderation <- moderation
  if (!is.null(frailty_table)) { if (!is.null(frailty_prior)) attr(frailty_table, "prior") <- frailty_prior; out$frailty <- frailty_table }
  out
}

condition_levels <- function(col) {
  if (is.factor(col)) {
    present <- unique(as.character(col[!is.na(col)]))
    return(levels(col)[levels(col) %in% present])
  }
  sort(unique(as.character(col[!is.na(col)])))
}

enumerate_pairs <- function(from, to, all_types, family, parent = NULL) {
  if (is.character(from) && length(from) == 1L && is.character(to) && length(to) == 1L)
    return(list(c(from, to)))
  # Kontextual is directional too, and needs `to` in the context
  if (!is.null(parent)) {
    from_types <- if (is.null(from)) all_types else unique(from)
    to_types <- if (is.null(to)) parent else unique(to)
    return(unlist(lapply(from_types, function(f) lapply(to_types, function(t) c(f, t))),
                  recursive = FALSE))
  }
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

pair_model_data <- function(ctx, f, t) {
  from_code <- ctx$type_index[[f]]; to_code <- ctx$type_index[[t]]
  if (ctx$family == "binomial" && !is.null(ctx$parent))
    dataset_kontextual_binomial_model_data(ctx$data, from_code, to_code)
  else if (ctx$family == "binomial")
    dataset_binomial_model_data(ctx$data, from_code, to_code)
  else if (!is.null(ctx$parent))
    dataset_kontextual_model_data(ctx$data, from_code, to_code)
  else if (ctx$null == "random_labelling") {
    # REF cells with no candidate neighbour have offset 0 and carry no information
    md <- dataset_rl_model_data(ctx$data, from_code, to_code)
    keep <- md$density > 0
    lapply(md, `[`, keep)
  }
  else if (!is.null(ctx$sigma))
    dataset_inhom_model_data(ctx$data, ctx$areas, from_code, to_code, ctx$edge_correct)
  else
    dataset_poisson_model_data(ctx$data, ctx$areas, from_code, to_code)
}

fit_one_pair <- function(ctx, f, t, compute_diagnostics = FALSE) {
  md <- pair_model_data(ctx, f, t)

  group <- ctx$image_group[md$image + 1L]
  present <- unique(group)
  if (length(present) < 2L) {
    d <- diagnose_missing(f, t, setdiff(0:1, present), ctx)
    return(skip_pair(f, t, d$reason, d$message))
  }

  n <- md$n
  trials <- if (!is.null(md$trials)) md$trials else ctx$k
  b <- boundary_reason(f, t, n, group, ctx, trials)
  if (ctx$estimator == "mle" && !is.null(b$reason))
    return(skip_pair(f, t, b$reason, b$message))

  cluster <- ctx$image_cluster[md$image + 1L]
  if (ctx$frailty) return(frailty_outcome(md, ctx, f, t, b))
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
  else if (!is.null(md$trials))
    fit_pair_binomial_trials_cpp(cluster, md$image, group, n, md$trials, md$p0, ctx$estimator,
                                 ctx$cr2_method)
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

# Conditions where the MLE is infinite: all counts zero, or (binomial) all trials.
boundary_reason <- function(f, t, n, group, ctx, trials = ctx$k) {
  floor_g <- ctx$levels[vapply(0:1, function(g) !any(n[group == g] > 0), logical(1))]
  full <- n == trials
  ceil_g <- if (ctx$family == "binomial")
    ctx$levels[vapply(0:1, function(g) all(full[group == g]), logical(1))] else character(0)
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
