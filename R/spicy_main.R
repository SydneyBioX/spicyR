#' Test for changes in the co-localisation of cell types between conditions
#'
#' `spicy()` tests, for every ordered pair of cell types `from` → `to`, whether the co-localisation of
#' the two types differs between conditions, or is associated with survival. A pair asks whether `to`
#' cells are placed near `from` cells more than other cells are.
#'
#' **`method = "cell"` (the default).** For each `to` cell, the number of `from` cells within radius `r`
#' is compared with its exact expectation if the `to` cells were a random choice among the cells that
#' are not `from` cells in the same image. The effect is the **excess**: the number of extra `from` cells
#' within `r` of each `to` cell. Images are combined within patients and patients within conditions by
#' a frailty GEE, and the difference between conditions is tested with a CR2 cluster-robust variance on
#' Satterthwaite degrees of freedom, with **patients (`subject`) as the units**. By default the
#' difference is adjusted for how common the `from` type is in each image (the log of its share of all
#' cells) and for any `covariates`, so that a change in abundance alone does not appear as a change in
#' co-localisation. The unadjusted test is reported alongside (`unadjusted_*` columns).
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
#'   the excess, and the effect of each is reported as `<column>_effect` and `<column>_p_value`;
#'   survival: added to the null Cox model).
#' @param imageID,cellType,spatialCoords Column names of the image, cell type and coordinates.
#' @param r Radius (or radii) of the neighbourhood, in the units of the coordinates. Cell method: one
#'   radius (default 50), or several to be combined by `combine`. Image method: the radii of the L
#'   function (default 20, 50 and 100).
#' @param from,to Cell types to test (all ordered pairs by default).
#' @param method `"cell"` (spicyR Cell, the default) or `"image"` (the original spicyR test).
#' @param k Cell method: use the `k` nearest neighbours instead of a radius.
#' @param combine Cell method with several radii: `"maxT"` (max-T with the sandwich correlation across
#'   radii) or `"cauchy"` (Cauchy combination).
#' @param adjustAbundance Cell method: adjust the test for the log share of the `from` type in each image
#'   (default `TRUE`). Its effect is reported as `abundance_effect`. `FALSE` gives the test without it.
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
#'   condition (at the average abundance and covariates), the difference, its standard error, df,
#'   p-value and BH-adjusted p-value, the frailty variance, what the test was adjusted for
#'   (`adjusted_for`), the effect and p-value of each adjustment, and the unadjusted test
#'   (`unadjusted_difference`, `unadjusted_p_value`, `unadjusted_p_adj`).
#' @references Canete NP et al. (2022). spicyR: spatial analysis of in situ cytometry data in R.
#'   Bioinformatics 38(11), 3099-3105. \doi{10.1093/bioinformatics/btac268}
#'
#'   Bell RM, McCaffrey DF (2002). Bias reduction in standard errors for linear regression with multi-stage
#'   samples. Survey Methodology 28(2), 169-181.
#'
#'   Pustejovsky JE, Tipton E (2018). Small-sample methods for cluster-robust variance estimation and hypothesis
#'   testing in fixed effects models. Journal of Business & Economic Statistics 36(4), 672-683.
#'   \doi{10.1080/07350015.2016.1247004}
#'
#'   Paule RC, Mandel J (1982). Consensus values and weighting factors. Journal of Research of the National
#'   Bureau of Standards 87(5), 377-385. \doi{10.6028/jres.087.022}
#'
#'   Hartung J, Knapp G (2001). A refined method for the meta-analysis of controlled clinical trials with binary
#'   outcome. Statistics in Medicine 20(24), 3875-3889. \doi{10.1002/sim.1009}
#'
#'   Liu Y, Xie J (2020). Cauchy combination test: a powerful test with analytic p-value calculation under
#'   arbitrary dependency structures. Journal of the American Statistical Association 115(529), 393-402.
#'   \doi{10.1080/01621459.2018.1554485}
#' @examples
#' data("diabetesData")
#' # spicyR Cell: patients ("case") are the units
#' # extra Tc cells within 50 units of each Th cell and of each beta cell
#' res <- spicy(diabetesData, condition = "stage", subject = "case", r = 50,
#'              from = "Tc", to = c("Th", "beta"))
#' topPairs(res)
#' res$cellResults
#'
#' # the original image-level test
#' resImage <- spicy(diabetesData, condition = "stage", subject = "case",
#'                   from = "Tc", to = "Th", method = "image")
#' topPairs(resImage)
#' @aliases spicy
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
                  adjustAbundance = TRUE,
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
  # from -> to: extra `from` cells within r of each `to` cell, beyond the `to` cells being a random subset of the
  # cells that are not `from` cells. This is the core's own (counted, centre) order, so pairs pass through as given.
  pairs <- enumerate_pairs(from, to, types, "binomial")
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
    # one row per image, NA where a covariate is missing (those images are left out of the covariate model)
    Z_extra <- stats::model.matrix(stats::reformulate(covariates), stats::model.frame(stats::reformulate(covariates), pheno,
                                   na.action = stats::na.pass))[, -1, drop = FALSE]
    Z_extra <- sweep(Z_extra, 2, colMeans(Z_extra, na.rm = TRUE))
  }

  radii <- if (is.null(k)) sort(unique(r)) else NA
  if (length(radii) > 1L && !survival && length(.cell_levels(cells[[condition]])) > 2L)
    stop("several radii are supported for two conditions; give one `r`.", call. = FALSE)
  if (survival) {
    res <- .cell_survival(ctx, pairs, radii, k, pheno, covariates, labelClustering, cores, adjustAbundance)
  } else {
    per_r <- lapply(radii, function(rr) {
      g <- .cell_graph(ctx, pairs, r = if (is.na(rr)) NULL else rr, k = k, label_clustering = labelClustering,
                       n_threads = cores)
      lapply(pairs, function(p) .cell_pair_test(ctx, g, p[1], p[2], frailty, variance, adjustAbundance, Z_extra))
    })
    res <- .cell_combine(per_r, radii, ctx, adjustAbundance || !is.null(covariates), combine)
  }
  .cell_results(res, ctx, pheno, condition, subject, survival, radii, k)
}

## Several radii: per-radius tables, and one row per pair with the combined p-value (the main test and,
## when adjusted, the unadjusted test). The other columns are those of the radius chosen by max-T.
.cell_combine <- function(per_r, radii, ctx, adjusted, combine) {
  tabs <- lapply(per_r, .cell_table, ctx = ctx, adjusted = adjusted)
  if (length(radii) == 1L) return(list(table = tabs[[1]], fits = per_r[[1]]))
  long <- do.call(rbind, Map(function(t, rr) if (!is.null(t)) cbind(r = rr, t), tabs, radii))
  keys <- unique(long[, c("from", "to")])
  comb <- function(tests) {
    ok <- !vapply(tests, is.null, TRUE)
    if (!any(ok)) return(NULL)
    tt <- vapply(tests[ok], function(z) z$difference / z$se, 0); df <- vapply(tests[ok], `[[`, 0, "df")
    pv <- vapply(tests[ok], `[[`, 0, "p")
    if (combine == "maxT") {
      mt <- stats_max_t(do.call(cbind, lapply(tests[ok], `[[`, "influence")), tt, df)
      list(best = which(ok)[mt$best], p = mt$p, p_best = pv[mt$best])
    } else { b <- which.min(pv); list(best = which(ok)[b], p = stats_cauchy(pv), p_best = pv[b]) }
  }
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    f <- keys$from[i]; t <- keys$to[i]
    fit_at <- lapply(per_r, function(fs) { o <- fs[[which(vapply(fs, function(z) z$from == f && z$to == t, TRUE))]]
      if (o$ok) o else NULL })
    m <- comb(lapply(fit_at, function(o) o$test))
    if (is.null(m)) return(NULL)
    row <- long[long$from == f & long$to == t & long$r == radii[m$best], , drop = FALSE]
    row$p_value <- m$p; row$p_value_best_radius <- m$p_best
    if (adjusted) {
      u <- comb(lapply(fit_at, function(o) o$unadjusted))
      row$unadjusted_p_value <- if (is.null(u)) NA_real_ else u$p
    }
    row
  })
  tab <- do.call(rbind, rows); tab$p_adj <- stats::p.adjust(tab$p_value, "BH")
  if (adjusted) tab$unadjusted_p_adj <- stats::p.adjust(tab$unadjusted_p_value, "BH")
  for (col in c("unadjusted_difference", "unadjusted_se", "unadjusted_df")) tab[[col]] <- NULL
  rownames(tab) <- paste(tab$from, tab$to, sep = "__")
  list(table = .cell_order_columns(tab), fits = per_r[[which.min(abs(radii - stats::median(radii)))]], radius_table = long)
}

## Survival: score test and the hazard ratio of the shrunken excess, per pair.
.cell_survival <- function(ctx, pairs, radii, k, pheno, covariates, labelClustering, cores, adjust) {
  sv <- .cell_survival_setup(ctx, pheno, covariates)
  rr <- radii[1]
  if (length(radii) > 1L) message("survival uses one radius; using r = ", rr, ".")
  g <- .cell_graph(ctx, pairs, r = if (is.na(rr)) NULL else rr, k = k, label_clustering = labelClustering, n_threads = cores)
  fits <- lapply(pairs, function(p) .cell_survival_test(ctx, .cell_rows(ctx, g, p[1], p[2]), p[1], p[2], sv, adjust, covariates))
  list(table = .cell_survival_table(fits, adjust, covariates), fits = fits)
}

## The null Cox model (covariates only) and its martingale residuals, one per patient.
.cell_survival_setup <- function(ctx, pheno, covariates) {
  unit_first <- match(seq_along(ctx$unit_labels) - 1L, ctx$image_unit)
  time <- pheno$.time[unit_first]; event <- as.integer(pheno$.event[unit_first])
  if (any(tapply(pheno$.time, ctx$image_unit, function(z) length(unique(z))) > 1L))
    stop("the survival outcome must be constant within each subject.", call. = FALSE)
  W <- if (is.null(covariates)) matrix(0, length(time), 0) else {
    pu <- pheno[unit_first, , drop = FALSE]
    stats::model.matrix(stats::reformulate(covariates), stats::model.frame(stats::reformulate(covariates), pu,
                        na.action = stats::na.pass))[, -1, drop = FALSE] }
  ok_u <- is.finite(time) & !is.na(event) & stats::complete.cases(W)
  null <- stats_cox_fit(time[ok_u], event[ok_u], W[ok_u, , drop = FALSE])
  if (!null$ok) stop("the null Cox model did not converge.", call. = FALSE)
  M <- rep(NA_real_, length(time)); M[ok_u] <- null$martingale
  list(M = M, time = time, event = event)
}

## The survival test of one pair on given image rows.
.cell_survival_test <- function(ctx, rows, f, t, sv, adjust, covariates) {
  m <- length(ctx$unit_labels)
  u <- stats_survival_test(rows, rows$unit, m, sv$M, sv$time, sv$event, numeric(0))
  x <- if (adjust) .cell_share(ctx, rows, f) else numeric(0)
  s <- if (length(x) && stats::var(x) > 0) stats_survival_test(rows, rows$unit, m, sv$M, sv$time, sv$event, x) else u
  list(from = f, to = t, ok = s$ok, reason = s$reason, rows = rows, surv = s, unadjusted_surv = u,
       adjusted_for = paste(c(if (!identical(s, u)) "abundance", if (!is.null(covariates)) "covariates"), collapse = "+"))
}

.cell_survival_table <- function(fits, adjust, covariates) {
  ok <- vapply(fits, `[[`, TRUE, "ok")
  tab <- do.call(rbind, lapply(fits[ok], function(o) { s <- o$surv
    row <- data.frame(from = o$from, to = o$to, score_coefficient = s$score_coef, score_se = s$score_se, score_df = s$score_df,
               p_value = s$score_p, hazard_ratio_sd = s$hr_sd, log_hr_sd = s$log_hr_sd, log_hr_se = s$hr_se,
               hr_p_value = s$hr_p, log_hr_per_cell = s$log_hr_unit, tau2 = s$tau2,
               adjusted_for = if (nzchar(o$adjusted_for)) o$adjusted_for else "none",
               unadjusted_p_value = if (isTRUE(o$unadjusted_surv$ok)) o$unadjusted_surv$score_p else NA_real_,
               unadjusted_hazard_ratio_sd = if (isTRUE(o$unadjusted_surv$ok)) o$unadjusted_surv$hr_sd else NA_real_,
               stringsAsFactors = FALSE)
    if (!is.null(o$parent)) row <- cbind(row[1:2], parent = o$parent, row[-(1:2)])
    row }))
  if (!is.null(tab)) {
    tab$p_adj <- stats::p.adjust(tab$p_value, "BH"); tab$unadjusted_p_adj <- stats::p.adjust(tab$unadjusted_p_value, "BH")
    if (!adjust && is.null(covariates)) tab[c("adjusted_for", grep("^unadjusted_", names(tab), value = TRUE))] <- NULL
    rownames(tab) <- .cell_labels(tab)
  }
  tab
}

## Assemble a SpicyResults object that topPairs(), signifPlot(), spicyBoxPlot() and bind() understand.
.cell_results <- function(res, ctx, pheno, condition, subject, survival, radii, k) {
  tab <- res$table
  num <- vapply(tab, is.double, TRUE); tab[num] <- lapply(tab[num], function(z) { z[is.nan(z)] <- NA_real_; z })
  if (is.null(tab)) stop("no pair could be tested (each condition needs at least two patients with both cell types).",
                         call. = FALSE)
  out <- list(method = "cell", cellResults = tab)
  if (!survival && length(ctx$levels) > 2L) {
    # wide matrices: one column per level against the reference, as the image method's model terms
    key <- .cell_labels(tab); labels <- unique(key)
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
  labels <- .cell_labels(tab)
  if (!is.null(res$radius_table)) out$radiusResults <- res$radius_table
  out$comparisons <- data.frame(from = tab$from, to = tab$to, labels = labels)
  if (!is.null(tab$parent)) { out$comparisons$parent <- tab$parent; out$isKontextual <- TRUE }
  out$pairwiseAssoc <- .cell_image_excess(res$fits, ctx)[labels]
  out$imageWeights <- .cell_image_weight(res$fits, ctx)[labels]
  out$imageIDs <- ctx$image_labels
  out$imageID <- ctx$image_labels
  if (!is.null(subject)) out$subject <- as.character(pheno[[subject]])
  out$nCells <- table(ctx$df$imageID, ctx$df$cellType)
  out$alternateResult <- FALSE
  out$r <- if (is.null(k)) radii else NULL
  out$k <- k
  methods::new("SpicyResults", out)
}
