#' Test for changes in the co-localisation of cell types between conditions
#'
#' `spicy()` tests, for every ordered pair of cell types *from* → *to*, whether the co-localisation of
#' the two types differs between conditions, or is associated with survival. The direction follows the
#' spatial-statistics convention for cross-type statistics: *from* is the type whose neighbourhoods are
#' examined, *to* the type counted in them.
#'
#' **`method = "cell"` (the default, spicyR Cell).** For each `to` cell, the number of `from` cells
#' within radius `r` is compared with its exact expectation under random labelling of the observed
#' cells. The effect is the **excess**: the number of extra `to` cells within `r` of each `from` cell. Images are
#' combined within patients and patients within conditions by a frailty GEE, and the difference
#' between conditions is tested with a CR2 cluster-robust variance on Satterthwaite degrees of
#' freedom, with **patients (`subject`) as the units**. Results for each pair also include the
#' difference at equal availability of the `from` type (adjusted for its share of all cells), which
#' guards against changes in abundance being read as changes in attraction.
#'
#' **Changed in version 2.0:** up to spicyR 1.99.0, `spicy_glm(effect = "excess")` counted extra `from`
#' cells around each `to` cell. `spicy()` now uses the convention above.
#'
#' **`method = "image"` (the original spicyR test).** A per-image L-function summary of each pair is
#' compared between conditions with a weighted linear model, or a mixed model when `subject` is given
#' (Canete et al. 2022).
#'
#' @param cells A `data.frame`, `SingleCellExperiment` or `SpatialExperiment` with one row (column) per
#'   cell.
#' @param condition The column of the image-level condition: two groups, or a `survival::Surv` column
#'   for a survival outcome.
#' @param subject The column of the patient (unit) of each image. Images of one patient are combined;
#'   if omitted, every image is its own patient.
#' @param covariates Image- or patient-level columns to adjust for (cell method: added to the design of
#'   the excess; survival: added to the null Cox model).
#' @param imageID,cellType,spatialCoords Column names of the image, cell type and coordinates.
#' @param r Radius (or radii) of the neighbourhood, in the units of the coordinates. Cell method: one
#'   radius, or several to be combined by `combine`. Image method: the radii of the L function.
#' @param from,to Cell types to test (all ordered pairs by default).
#' @param method `"cell"` (spicyR Cell, the default) or `"image"` (the original spicyR test).
#' @param k Cell method: use the `k` nearest neighbours instead of a radius.
#' @param combine Cell method with several radii: `"maxT"` (max-T with the sandwich correlation across
#'   radii) or `"cauchy"` (Cauchy combination).
#' @param availability Cell method: also report the difference at equal availability of the `to`
#'   type, i.e. adjusted for its share of all cells (`adjusted_*` columns).
#' @param variance Cell method: `"cr2"` (CR2 on Satterthwaite df, the default) or `"hartung_knapp"`
#'   (for very few patients: the model-based variance floored at CR2, on m - 2 df).
#' @param frailty,labelClustering Cell method: the patient frailty and the label-clustering inflation
#'   of the within-image variance (both on by default).
#' @param ref Cell method: the reference level of `condition`.
#' @param cores Number of threads (cell method) or cores (image method).
#' @param ... Arguments of the image method: `sigma`, `alternateResult`, `minLambda`, `weights`,
#'   `weightsByPair`, `weightFactor`, `weightZThreshold`, `window`, `window.length`, `edgeCorrect`,
#'   `includeZeroCells`, `verbose`, `BPPARAM`. Supplying `alternateResult` selects the image method.
#' @return A `SpicyResults` object. `topPairs()`, `signifPlot()`, `spicyBoxPlot()` and `bind()` work
#'   for both methods. For the cell method, `$cellResults` holds the full table: the excess in each
#'   condition, the difference, its standard error, df, p-value and BH-adjusted p-value, the frailty
#'   variance, and the availability-adjusted difference and p-value.
#' @references Canete NP et al. (2022). spicyR: spatial analysis of in situ cytometry data in R.
#'   Bioinformatics 38(11), 3099-3105.
#' @examples
#' data("diabetesData")
#' # spicyR Cell: patients ("case") are the units
#' # extra Th and beta cells within 50 units of each Tc cell
#' res <- spicy(diabetesData, condition = "stage", subject = "case", r = 50,
#'              from = "Tc", to = c("Th", "beta"))
#' topPairs(res)
#' res$cellResults
#'
#' # the original image-level test
#' resImage <- spicy(diabetesData, condition = "stage", subject = "case",
#'                   from = "Tc", to = "Th", method = "image")
#' topPairs(resImage)
#' @aliases spicy spicy,spicy-method
#' @export
spicy <- function(cells,
                  condition,
                  subject = NULL,
                  covariates = NULL,
                  imageID = "imageID",
                  cellType = "cellType",
                  spatialCoords = c("x", "y"),
                  r = NULL,
                  from = NULL,
                  to = NULL,
                  method = c("cell", "image"),
                  k = NULL,
                  combine = c("maxT", "cauchy"),
                  availability = TRUE,
                  variance = c("cr2", "hartung_knapp"),
                  frailty = TRUE,
                  labelClustering = TRUE,
                  ref = NULL,
                  cores = 1,
                  ...) {
  dots <- list(...)
  if (missing(method) && (!is.null(dots$alternateResult) || !is.null(dots$sigma))) {
    method <- "image"
    message("`alternateResult` / `sigma` given: using method = \"image\" (the original spicyR test).")
  }
  method <- match.arg(method)
  if (method == "image") {
    return(do.call(.spicyImage, c(list(cells = cells, condition = condition, subject = subject, covariates = covariates,
                                       imageID = imageID, cellType = cellType, spatialCoords = spatialCoords,
                                       r = r, from = from, to = to, cores = cores), dots)))
  }
  if (length(dots)) stop("arguments not used by method = \"cell\": ", paste(names(dots), collapse = ", "),
                         ". They belong to method = \"image\".", call. = FALSE)
  combine <- match.arg(combine); variance <- match.arg(variance)
  if (is.null(r) && is.null(k)) r <- 50
  if (!is.null(k) && length(r)) { message("`k` given: using the k nearest neighbours; `r` is ignored."); r <- NULL }

  cells <- .format_data(cells, imageID, cellType, spatialCoords, FALSE)
  if (!is.null(subject) && !subject %in% names(cells)) stop("`subject` column not found.", call. = FALSE)
  survival <- inherits(cells[[condition]], "Surv")
  types <- unique(as.character(cells$cellType))
  bad <- setdiff(c(from, to), types)
  if (length(bad)) stop("cell type not found: ", paste(bad, collapse = ", "), call. = FALSE)
  # the user's from -> to: extra `to` cells within r of each `from` cell (the spatial-statistics convention).
  # The core works in (counted, centre) order, so the internal pairs are reversed; .cell_results swaps back.
  pairs <- lapply(enumerate_pairs(from, to, types, "binomial"), rev)
  if (survival) {
    sv <- cells[[condition]]; cells$.time <- sv[, 1]; cells$.event <- sv[, 2]; cells[[condition]] <- NULL
  }
  ctx <- .cell_context(cells, if (survival) NULL else condition, subject, "imageID", "cellType", c("x", "y"),
                       ref = ref, survival = survival)
  pheno <- ctx$df[ctx$first, , drop = FALSE]

  Z_extra <- NULL
  if (!is.null(covariates)) {
    miss <- setdiff(covariates, names(pheno))
    if (length(miss)) stop("covariates not found: ", paste(miss, collapse = ", "), call. = FALSE)
    Z_extra <- stats::model.matrix(stats::reformulate(covariates), pheno)[, -1, drop = FALSE]
    Z_extra <- sweep(Z_extra, 2, colMeans(Z_extra))
  }

  radii <- if (is.null(k)) sort(unique(r)) else NA
  if (length(radii) > 1L && !survival && length(.cell_levels(cells[[condition]])) > 2L)
    stop("several radii are supported for two conditions; give one `r`.", call. = FALSE)
  if (survival) {
    res <- .cell_survival(ctx, pairs, radii, k, pheno, covariates, labelClustering, cores)
  } else {
    per_r <- lapply(radii, function(rr) {
      g <- .cell_graph(ctx, pairs, r = if (is.na(rr)) NULL else rr, k = k, label_clustering = labelClustering,
                       n_threads = cores)
      lapply(pairs, function(p) .cell_pair_test(ctx, g, p[1], p[2], frailty, variance, availability, Z_extra))
    })
    res <- .cell_combine(per_r, radii, ctx, availability, covariates, combine)
  }
  .cell_results(res, ctx, pheno, condition, subject, survival, radii, k)
}

## Several radii: per-radius tables, and one row per pair with the combined p-value.
.cell_combine <- function(per_r, radii, ctx, availability, covariates, combine) {
  tabs <- lapply(per_r, .cell_table, ctx = ctx, availability = availability, covariates = covariates)
  if (length(radii) == 1L) return(list(table = tabs[[1]], fits = per_r[[1]]))
  long <- do.call(rbind, Map(function(t, rr) if (!is.null(t)) cbind(r = rr, t), tabs, radii))
  keys <- unique(long[, c("from", "to")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    f <- keys$from[i]; t <- keys$to[i]
    tests <- lapply(per_r, function(fs) { o <- fs[[which(vapply(fs, function(z) z$from == f && z$to == t, TRUE))]]
      if (o$ok) o$test else NULL })
    ok <- !vapply(tests, is.null, TRUE)
    if (!any(ok)) return(NULL)
    tt <- vapply(tests[ok], function(z) z$difference / z$se, 0); df <- vapply(tests[ok], `[[`, 0, "df")
    pv <- vapply(tests[ok], `[[`, 0, "p")
    if (combine == "maxT") {
      infl <- do.call(cbind, lapply(tests[ok], `[[`, "influence"))
      mt <- stats_max_t(infl, tt, df); best <- mt$best; p <- mt$p
    } else { best <- which.min(pv); p <- stats_cauchy(pv) }
    z <- tests[ok][[best]]
    data.frame(from = f, to = t, r = radii[ok][best], excess_ref = z$coef_ref, excess_comp = z$coef_comp,
               excess_difference = z$difference, se = z$se, df = z$df, p_value = p, tau2 = z$tau2,
               p_value_best_radius = pv[best], stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows); tab$p_adj <- stats::p.adjust(tab$p_value, "BH")
  if (availability) {
    adj <- tapply(long$adjusted_p_value, paste(long$from, long$to, sep = "__"), function(p) stats_cauchy(p[is.finite(p)]))
    tab$adjusted_p_value <- as.numeric(adj[paste(tab$from, tab$to, sep = "__")])
    tab$adjusted_p_adj <- stats::p.adjust(tab$adjusted_p_value, "BH")
  }
  rownames(tab) <- paste(tab$from, tab$to, sep = "__")
  list(table = tab, fits = per_r[[which.min(abs(radii - stats::median(radii)))]], radius_table = long)
}

## Survival: score test and the hazard ratio of the shrunken excess, per pair.
.cell_survival <- function(ctx, pairs, radii, k, pheno, covariates, labelClustering, cores) {
  unit_first <- match(seq_along(ctx$unit_labels) - 1L, ctx$image_unit)
  time <- pheno$.time[unit_first]; event <- as.integer(pheno$.event[unit_first])
  if (any(tapply(pheno$.time, ctx$image_unit, function(z) length(unique(z))) > 1L))
    stop("the survival outcome must be constant within each subject.", call. = FALSE)
  W <- if (is.null(covariates)) matrix(0, length(time), 0) else
    stats::model.matrix(stats::reformulate(covariates), pheno[unit_first, , drop = FALSE])[, -1, drop = FALSE]
  ok_u <- is.finite(time) & !is.na(event)
  null <- stats_cox_fit(time[ok_u], event[ok_u], W[ok_u, , drop = FALSE])
  if (!null$ok) stop("the null Cox model did not converge.", call. = FALSE)
  M <- rep(NA_real_, length(time)); M[ok_u] <- null$martingale
  rr <- radii[1]
  if (length(radii) > 1L) message("survival uses one radius; using r = ", rr, ".")
  g <- .cell_graph(ctx, pairs, r = if (is.na(rr)) NULL else rr, k = k, label_clustering = labelClustering, n_threads = cores)
  fits <- lapply(pairs, function(p) {
    rows <- .cell_rows(ctx, g, p[1], p[2])
    s <- stats_survival_test(rows, rows$unit, length(ctx$unit_labels), M, time, event)
    list(from = p[1], to = p[2], ok = s$ok, reason = s$reason, rows = rows, surv = s) })
  ok <- vapply(fits, `[[`, TRUE, "ok")
  tab <- do.call(rbind, lapply(fits[ok], function(o) { s <- o$surv
    data.frame(from = o$from, to = o$to, score_coefficient = s$score_coef, score_se = s$score_se, score_df = s$score_df,
               p_value = s$score_p, hazard_ratio_sd = s$hr_sd, log_hr_sd = s$log_hr_sd, log_hr_se = s$hr_se,
               hr_p_value = s$hr_p, log_hr_per_cell = s$log_hr_unit, tau2 = s$tau2, stringsAsFactors = FALSE) }))
  if (!is.null(tab)) { tab$p_adj <- stats::p.adjust(tab$p_value, "BH"); rownames(tab) <- paste(tab$from, tab$to, sep = "__") }
  list(table = tab, fits = fits)
}

## Assemble a SpicyResults object that topPairs(), signifPlot(), spicyBoxPlot() and bind() understand.
.cell_swap <- function(d) {
  if (is.null(d)) return(d)
  f <- d$from; d$from <- d$to; d$to <- f
  rownames(d) <- if (!is.null(d$level)) paste(d$from, d$to, d$level, sep = "__") else if (anyDuplicated(paste(d$from, d$to))) NULL else paste(d$from, d$to, sep = "__")
  d
}

.cell_results <- function(res, ctx, pheno, condition, subject, survival, radii, k) {
  tab <- .cell_swap(res$table); res$radius_table <- .cell_swap(res$radius_table)
  if (is.null(tab)) stop("no pair could be tested (each condition needs at least two patients with both cell types).",
                         call. = FALSE)
  out <- list(method = "cell", cellResults = tab)
  if (!survival && length(ctx$levels) > 2L) {
    # wide matrices: one column per level against the reference, as the image method's model terms
    key <- paste(tab$from, tab$to, sep = "__"); labels <- unique(key)
    wide <- function(col) { d <- data.frame(row.names = labels)
      d[["(Intercept)"]] <- if (col == "excess_difference") tab$excess_ref[match(labels, key)] else NA_real_
      for (l in ctx$levels[-1]) { k <- tab$level == l; d[[paste0("condition", l)]] <- tab[[col]][k][match(labels, key[k])] }
      d }
    out$coefficient <- wide("excess_difference"); out$p.value <- wide("p_value"); out$se <- wide("se"); out$df <- wide("df")
    out$statistic <- out$coefficient; out$statistic[, -1] <- out$coefficient[, -1] / out$se[, -1]; out$statistic[, 1] <- NA
    out$condition <- factor(ctx$levels[ctx$image_group + 1L], levels = ctx$levels)
    tab <- tab[match(labels, key), ]
  } else {
  labels <- rownames(tab)
  term <- if (survival) "condition" else paste0("condition", ctx$levels[2])
  mk <- function(a, b) { d <- data.frame(a, b); names(d) <- c("(Intercept)", term); rownames(d) <- labels; d }
  if (survival) {
    out$coefficient <- mk(NA_real_, tab$log_hr_sd); out$p.value <- mk(NA_real_, tab$p_value)
    out$se <- mk(NA_real_, tab$score_se); out$statistic <- mk(NA_real_, tab$score_coefficient / tab$score_se)
    out$df <- mk(NA_real_, tab$score_df)
    out$survivalOutcome <- survival::Surv(pheno$.time, pheno$.event)
  } else {
    out$coefficient <- mk(tab$excess_ref, tab$excess_difference); out$p.value <- mk(NA_real_, tab$p_value)
    out$se <- mk(NA_real_, tab$se); out$statistic <- mk(NA_real_, tab$excess_difference / tab$se)
    out$df <- mk(NA_real_, tab$df)
    out$condition <- factor(ctx$levels[ctx$image_group + 1L], levels = ctx$levels)
  }
  }
  labels <- paste(tab$from, tab$to, sep = "__")
  if (!is.null(res$radius_table)) out$radiusResults <- res$radius_table
  out$comparisons <- data.frame(from = tab$from, to = tab$to, labels = labels)
  pa <- .cell_image_excess(res$fits, ctx)
  names(pa) <- vapply(strsplit(names(pa), "__", fixed = TRUE), function(z) paste(z[2], z[1], sep = "__"), "")
  out$pairwiseAssoc <- pa[labels]
  out$imageIDs <- ctx$image_labels
  out$imageID <- ctx$image_labels
  if (!is.null(subject)) out$subject <- as.character(pheno[[subject]])
  out$nCells <- table(ctx$df$imageID, ctx$df$cellType)
  out$alternateResult <- FALSE
  out$r <- if (is.null(k)) radii else NULL
  out$k <- k
  methods::new("SpicyResults", out)
}
