## The engine of Statial::kontextualTest(): Statial's Kontextual statistic, tested between conditions with spicyR's
## random-labelling null, frailty GEE and CR2 variance. Documented in Statial. Written with AI assistance (Claude, Anthropic), directed by the authors; see NEWS.

#' Engine of Statial's Kontextual test
#'
#' The computation behind `Statial::kontextualTest()`; use that function, whose help page describes the method.
#' For each triple of `parentDf` it compares Statial's Kontextual statistic in each image with its exact
#' expectation and variance when the `to` cells are a random choice among the parent's cells, and combines the
#' images as [spicy()] does (frailty GEE, CR2 variance, label-clustering factor, abundance adjustment).
#'
#' @param cells A SingleCellExperiment, SpatialExperiment or data frame with a row per cell.
#' @param parentDf A data frame of triples with columns `from`, `to`, `parent` (a list column of cell types
#'   containing `to`) and optionally `parent_name`, as made by `Statial::parentCombinations()`.
#' @param condition The column with each image's condition, or a `Surv` column of survival outcomes.
#' @param subject,covariates,r,from,to,imageID,cellType,spatialCoords,adjustAbundance,variance,frailty,labelClustering,ref
#'   As for [spicy()]; `r` is one radius.
#' @param edgeCorrect Correct the parent densities for the part of each disc outside the image's window.
#' @param window The window of each image for the edge correction: `"convex"` or `"rectangle"`.
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
                             labelClustering = TRUE,
                             edgeCorrect = TRUE,
                             window = c("convex", "rectangle"),
                             ref = NULL) {
  variance <- match.arg(variance); window <- match.arg(window)
  if (length(r) != 1L) stop("give one radius `r`.", call. = FALSE)
  cells <- .format_data(cells, imageID, cellType, spatialCoords, FALSE)
  if (!is.null(subject) && !subject %in% names(cells)) stop("`subject` column not found.", call. = FALSE)
  survival <- inherits(cells[[condition]], "Surv")
  types <- unique(as.character(cells$cellType))
  trip <- .kontextual_triples(parentDf, types, from, to)
  if (survival) {
    sv <- cells[[condition]]; cells$.time <- sv[, 1]; cells$.event <- sv[, 2]; cells[[condition]] <- NULL
  }
  ctx <- .cell_context(cells, if (survival) NULL else condition, subject, "imageID", "cellType", c("x", "y"), ref = ref,
                       survival = survival)
  pheno <- ctx$df[ctx$first, , drop = FALSE]
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
  out <- .cell_results(list(table = tab, fits = fits), ctx, pheno, condition, subject, survival, r, NULL)
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
