## The engine of Statial::kontextualTest(): Statial's Kontextual statistic, tested between conditions with spicyR's
## random-labelling null, frailty GEE and CR2 variance. Documented in Statial. Written with AI assistance (Claude, Anthropic), directed by the authors; see NEWS.

#' Engine of Statial's Kontextual test
#'
#' The computation behind `Statial::kontextualTest()`; use that function, whose help page describes the method.
#' For each triple of `parentDf` it compares Statial's Kontextual statistic in each image with its exact
#' expectation and variance when the `to` cells are a random choice among the parent's cells, and combines the
#' images as [spicy()] does (frailty GEE, each condition's CR2 variance, enough effective patients per condition,
#' abundance adjustment).
#'
#' @param cells A SingleCellExperiment, SpatialExperiment or data frame with a row per cell.
#' @param parentDf A data frame of triples with columns `from`, `to`, `parent` (a list column of cell types
#'   containing `to`) and optionally `parent_name`, as made by `Statial::parentCombinations()`.
#' @param condition The column with each image's condition, or a `Surv` column of survival outcomes.
#' @param subject,covariates,r,from,to,imageID,cellType,spatialCoords,adjustAbundance,variance,frailty,ref
#'   As for [spicy()]; `r` is one radius.
#' @param labelClustering Whether to multiply each image's variance by a factor for clustering of the `to` labels.
#'   Off by default: the factor puts most of the weight on a few images, and without it the test stays at or below
#'   its nominal level when the conditions are shuffled between patients.
#' @param edgeCorrect Correct the parent densities for the part of each disc outside the image's window.
#' @param window The window of each image for the edge correction: `"convex"` or `"rectangle"`.
#' @param cellWeights Optional cell weights, a vector or a soft membership matrix as for [spicy()]: with a matrix,
#'   `from` and `to` (and the `from`, `to` of `parentDf`) name its columns, while `parent` still lists cell types. A `from` cell enters with its
#'   weight, the score of a candidate counts with its `to` weight, and the null shuffles the `to` weights over the
#'   candidates. Pairs of a label with itself are not tested, and `labelClustering` is not available.
#'
#' @return A [SpicyResults][SpicyResults-class] object with one row per triple (`from`, `to`, `parent`).
#'
#' @examples
#' data("diabetesData")
#' parentDf <- data.frame(from = "alpha", to = c("Tc", "Th"), parent_name = "tcells")
#' parentDf$parent <- list(c("Tc", "Th"), c("Tc", "Th"))
#' res <- kontextualEngine(diabetesData, parentDf, condition = "stage", subject = "case", r = 50)
#'
#' @keywords internal
#' @export
kontextualEngine <- function(cells,
                             parentDf,
                             condition,
                             subject = NULL,
                             covariates = NULL,
                             r = 50,
                             from = NULL,
                             to = NULL,
                             imageID = "imageID",
                             cellType = "cellType",
                             spatialCoords = c("x", "y"),
                             adjustAbundance = TRUE,
                             variance = c("cr2", "hartung_knapp"),
                             frailty = TRUE,
                             labelClustering = FALSE,
                             edgeCorrect = TRUE,
                             window = c("convex", "rectangle"),
                             ref = NULL,
                             cellWeights = NULL) {
  variance <- match.arg(variance); window <- match.arg(window)
  if (length(r) != 1L) stop("give one radius `r`.", call. = FALSE)
  weighted <- !is.null(cellWeights)
  if (weighted && labelClustering) stop("`labelClustering` is not available with `cellWeights`.", call. = FALSE)
  # the weights follow the input order of the cells
  cells <- .format_data(cells, imageID, cellType, spatialCoords, FALSE, sort = !weighted)
  if (weighted) memb <- .membership(cellWeights, nrow(cells), cells$cellType)
  if (!is.null(subject) && !subject %in% names(cells)) stop("`subject` column not found.", call. = FALSE)
  survival <- inherits(cells[[condition]], "Surv")
  types <- unique(as.character(cells$cellType))
  trip <- if (weighted) .kontextual_triples_weighted(parentDf, types, memb$labels, from, to)
          else .kontextual_triples(parentDf, types, from, to)
  if (survival) {
    sv <- cells[[condition]]; cells$.time <- sv[, 1]; cells$.event <- sv[, 2]; cells[[condition]] <- NULL
  }
  ctx <- .cell_context(cells, if (survival) NULL else condition, subject, "imageID", "cellType", c("x", "y"), ref = ref,
                       survival = survival, dataset = TRUE)
  pheno <- ctx$pheno
  Z_extra <- NULL
  if (!is.null(covariates) && !survival) {
    miss <- setdiff(covariates, names(pheno))
    if (length(miss)) stop("covariates not found: ", paste(miss, collapse = ", "), call. = FALSE)
    Z_extra <- stats::model.matrix(stats::reformulate(covariates), stats::model.frame(stats::reformulate(covariates), pheno,
                                   na.action = stats::na.pass))[, -1, drop = FALSE]
    Z_extra <- sweep(Z_extra, 2, colMeans(Z_extra, na.rm = TRUE))
  }

  sv <- if (survival) .cell_survival_setup(ctx, pheno, covariates)
  dataset_build_radius_index(ctx$data, r)
  if (weighted)
    return(.kontextual_weighted(ctx, memb, trip, r, window, edgeCorrect, condition, subject, survival, sv,
                                frailty, variance, adjustAbundance, Z_extra, covariates, pheno))
  code <- function(z) match(z, ctx$type_labels) - 1L
  fits <- list()
  for (pn in unique(trip$parent_name)) {
    tr <- trip[trip$parent_name == pn, , drop = FALSE]
    dataset_build_context(ctx$data, code(tr$parent[[1]]), window, edgeCorrect)
    sums <- lapply(seq_len(nrow(tr)), function(i) dataset_kontextual_sums(ctx$data, code(tr$from[i]), code(tr$to[i])))
    # psi of a `to` type: the median over every cell type as `from` (the `to` type itself included, as in spicy()),
    # whichever triples were requested, so a triple's result does not depend on the others asked for
    psi <- if (labelClustering) {
      pp <- expand.grid(from = ctx$type_labels, to = unique(tr$to), stringsAsFactors = FALSE)
      raw <- do.call(rbind, lapply(seq_len(nrow(pp)), function(i) {
        j <- which(tr$from == pp$from[i] & tr$to == pp$to[i])
        (if (length(j)) sums[[j[1]]] else dataset_kontextual_sums(ctx$data, code(pp$from[i]), code(pp$to[i])))[6, ] }))
      stats_kontextual_clustering(ctx$data, code(pp$from), code(pp$to), raw, ctx$counts, 2 * r)
    } else matrix(numeric(0), 0, 0)
    fits <- c(fits, lapply(seq_len(nrow(tr)), function(i) {
      rows <- stats_kontextual_image_rows(sums[[i]], ctx$counts, code(tr$from[i]), code(tr$to[i]), psi)
      rows$unit <- ctx$image_unit[rows$img + 1L]
      if (!survival) rows$group <- ctx$image_group[rows$img + 1L]
      o <- if (survival) .cell_survival_test(ctx, rows, tr$from[i], tr$to[i], sv, adjustAbundance, covariates)
           else .cell_rows_test(ctx, rows, tr$from[i], tr$to[i], frailty, variance, adjustAbundance, Z_extra)
      o$parent <- pn
      o }))
  }
  tab <- if (survival) .cell_survival_table(fits, adjustAbundance, covariates)
         else .cell_table(fits, ctx, adjustAbundance || !is.null(covariates))
  out <- .cell_results(list(table = tab, fits = fits), ctx, pheno, condition, subject, survival, r, NULL, effect = "kontextual")
  out$kontextual <- list(parents = stats::setNames(lapply(unique(trip$parent_name), function(pn)
    trip$parent[[match(pn, trip$parent_name)]]), unique(trip$parent_name)), edgeCorrect = edgeCorrect, window = window)
  out
}

## The (from, to, parent) triples of a parentDf (Statial::parentCombinations), restricted to the cell types present
## and to `from` / `to`.
.kontextual_triples <- function(parentDf, types, from, to) {
  need <- c("from", "to", "parent")
  if (!all(need %in% names(parentDf))) stop("`parentDf` needs columns from, to and parent.", call. = FALSE)
  d <- data.frame(from = as.character(parentDf$from), to = as.character(parentDf$to), stringsAsFactors = FALSE)
  d$parent <- lapply(parentDf$parent, function(p) intersect(as.character(unlist(p)), types))
  d$parent_name <- if (!is.null(parentDf$parent_name)) as.character(parentDf$parent_name)
                   else vapply(d$parent, function(p) paste(sort(p), collapse = "+"), "")
  keyed <- vapply(d$parent, function(p) paste(sort(p), collapse = "\r"), "")
  if (any(tapply(keyed, d$parent_name, function(z) length(unique(z))) > 1L))
    stop("each parent_name must name one parent population.", call. = FALSE)
  keep <- d$from %in% types & d$to %in% types & vapply(seq_len(nrow(d)), function(i) d$to[i] %in% d$parent[[i]], TRUE)
  if (!is.null(from)) keep <- keep & d$from %in% from
  if (!is.null(to)) keep <- keep & d$to %in% to
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) stop("no triple to test: each `to` must be a cell type of the data in its parent.", call. = FALSE)
  d
}

## The triples of a parentDf with soft membership: `from` and `to` are labels (columns of the weights), `parent` lists
## cell types (restricted to those present); pairs of a label with itself are left out.
.kontextual_triples_weighted <- function(parentDf, types, labels, from, to) {
  need <- c("from", "to", "parent")
  if (!all(need %in% names(parentDf))) stop("`parentDf` needs columns from, to and parent.", call. = FALSE)
  d <- data.frame(from = as.character(parentDf$from), to = as.character(parentDf$to), stringsAsFactors = FALSE)
  d$parent <- lapply(parentDf$parent, function(p) intersect(as.character(unlist(p)), types))
  d$parent_name <- if (!is.null(parentDf$parent_name)) as.character(parentDf$parent_name)
                   else vapply(d$parent, function(p) paste(sort(p), collapse = "+"), "")
  bad <- setdiff(c(d$from, d$to), labels)
  if (length(bad)) stop("not a column of `cellWeights`: ", paste(unique(bad), collapse = ", "), call. = FALSE)
  keep <- d$from != d$to & lengths(d$parent) > 0
  if (!is.null(from)) keep <- keep & d$from %in% from
  if (!is.null(to)) keep <- keep & d$to %in% to
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) stop("no triple to test.", call. = FALSE)
  d
}

## kontextualEngine() with soft membership: per parent, the context, then each triple's weighted Kontextual sums
## (Dataset::kontextual_sums_weighted), its image rows and the test, as for cell types.
.kontextual_weighted <- function(ctx, memb, trip, r, window, edgeCorrect, condition, subject, survival, sv, frailty,
                                 variance, adjustAbundance, Z_extra, covariates, pheno) {
  code <- function(z) match(z, ctx$type_labels) - 1L
  wctx <- .weighted_context(ctx, memb)
  col <- function(l) { a <- numeric(memb$nCells); sel <- memb$j == match(l, memb$labels); a[memb$i[sel]] <- memb$x[sel]; a[ctx$ord] }
  fits <- list()
  for (pn in unique(trip$parent_name)) {
    tr <- trip[trip$parent_name == pn, , drop = FALSE]
    dataset_build_context(ctx$data, code(tr$parent[[1]]), window, edgeCorrect)
    fits <- c(fits, lapply(seq_len(nrow(tr)), function(i) {
      s <- dataset_kontextual_sums_weighted(ctx$data, col(tr$from[i]), col(tr$to[i]), FALSE)
      rows <- stats_kontextual_image_rows_weighted(s)
      rows$unit <- ctx$image_unit[rows$img + 1L]
      if (!survival) rows$group <- ctx$image_group[rows$img + 1L]
      o <- if (survival) .cell_survival_test(wctx, rows, tr$from[i], tr$to[i], sv, adjustAbundance, covariates)
           else .cell_rows_test(wctx, rows, tr$from[i], tr$to[i], frailty, variance, adjustAbundance, Z_extra)
      o$parent <- pn
      o }))
  }
  tab <- if (survival) .cell_survival_table(fits, adjustAbundance, covariates)
         else .cell_table(fits, ctx, adjustAbundance || !is.null(covariates))
  out <- .cell_results(list(table = tab, fits = fits), ctx, pheno, condition, subject, survival, r, NULL, effect = "kontextual")
  out$kontextual <- list(parents = stats::setNames(lapply(unique(trip$parent_name), function(pn)
    trip$parent[[match(pn, trip$parent_name)]]), unique(trip$parent_name)), edgeCorrect = edgeCorrect, window = window)
  out
}
