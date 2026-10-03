## The Kontextual test: Statial's Kontextual statistic, tested between conditions with spicyR's random-labelling
## null, frailty GEE and CR2 variance. Written with AI assistance (Claude, Anthropic), directed by the authors; see NEWS.

#' Test for changes in Kontextual relationships between conditions
#'
#' `kontextualTest()` tests, for every triple `from` → `to` within `parent`, whether the Kontextual relationship
#' differs between conditions. It is the engine behind `Statial::kontextualTest()`, which most users call.
#'
#' Kontextual (Ameen et al. 2025) asks whether `from` cells are near `to` cells more than they are near the other
#' cells of a parent population that contains `to`. Its statistic weights each `from`-`to` pair within `r` by
#' \eqn{e_a \lambda_P(a) / \lambda_P(b)}, where \eqn{\lambda_P(x)} is the density of parent cells within `r` of
#' `x` (the cell itself included) and \eqn{e_a} the edge correction of the `from` cell; the per-image value is
#' the sum of the weights over the image's `from` exposure \eqn{\sum_a \lambda_P(a)}. `Statial::Kontextual()`
#' returns this as an L-function.
#'
#' The test does not compare these values with zero. In each image their exact expectation and variance if the
#' `to` cells were a random choice among the parent's cells (that are not `from` cells) are computed, so a change
#' in tissue composition between conditions does not by itself give a difference. The **Kontextual excess** of an
#' image is its weighted count beyond that expectation, per unit of `from` exposure. Images are combined within
#' patients and patients within conditions as in [spicy()]: a frailty GEE with a CR2 cluster-robust variance on
#' the patients, a label-clustering factor over the parent's cells, and by default an adjustment for the
#' abundance of `from` (its log share of the image's cells) and for `covariates`.
#'
#' @param cells A SingleCellExperiment, SpatialExperiment or data frame with a row per cell.
#' @param parentDf A data frame of triples with columns `from`, `to` and `parent` (a list column of cell types,
#'   containing `to`), and optionally `parent_name`, as made by `Statial::parentCombinations()`.
#' @param parentList Alternatively, a named list of parent populations: every cell type is tested as `from`
#'   against every type of each parent as `to` (`from` = `to` included).
#' @param condition The column with each image's condition (constant within a patient), or a `Surv` column of
#'   survival outcomes: then the test is a score test of the Kontextual excess against the hazard, as [spicy()].
#' @param subject The column identifying patients, when patients have several images.
#' @param covariates Columns of image- or patient-level covariates to adjust for.
#' @param r The radius.
#' @param from,to Optional: keep only the triples with these `from` / `to` cell types.
#' @param imageID,cellType,spatialCoords The columns of the image, cell type and coordinates.
#' @param adjustAbundance Adjust for the abundance of `from` (as [spicy()]).
#' @param variance `"cr2"` (default) or `"hartung_knapp"`.
#' @param frailty Use the frailty (between-patient) variance component.
#' @param labelClustering Allow for clustering of the `to` labels within the parent (the label-clustering factor).
#' @param edgeCorrect Correct the parent densities for the part of each disc outside the image's window.
#' @param window The window of each image for the edge correction: `"convex"` (as Statial) or `"rectangle"`.
#' @param ref The reference level of `condition`.
#'
#' @return A [SpicyResults][SpicyResults-class] object with one row per triple (`from`, `to`, `parent`):
#'   [topPairs()], [signifPlot()] and [spicyBoxPlot()] work on it. The excess is on the scale of Kontextual's
#'   K function (an area, in the squared units of the coordinates): an image's Kontextual L value is
#'   \eqn{\sqrt{K/\pi} - r}, with \eqn{K = \pi r^2} expected under random labelling within the parent.
#'
#' @references Ameen F, Robertson N, Lin DM, Ghazanfar S, Patrick E (2025). Kontextual reframes analysis of
#'   spatial omics data and reveals consistent cell relationships across images. *Cell Reports Methods*.
#'
#' @examples
#' data("diabetesData")
#' res <- kontextualTest(diabetesData, parentList = list(immune = c("Tc", "Th", "macrophage")),
#'                       condition = "stage", subject = "case", r = 50, from = "alpha")
#' topPairs(res)
#'
#' @export
kontextualTest <- function(cells,
                           parentDf = NULL,
                           parentList = NULL,
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
  trip <- .kontextual_triples(parentDf, parentList, types, from, to)
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
    psi <- if (labelClustering)
      stats_kontextual_clustering(ctx$data, code(tr$from), code(tr$to), do.call(rbind, lapply(sums, function(s) s[6, ])),
                                  ctx$counts, 2 * r)
    else matrix(numeric(0), 0, 0)
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

## The (from, to, parent) triples of a parentDf (Statial::parentCombinations) or a named list of parents,
## restricted to the cell types present and to `from` / `to`.
.kontextual_triples <- function(parentDf, parentList, types, from, to) {
  if (is.null(parentDf) == is.null(parentList)) stop("give one of `parentDf` and `parentList`.", call. = FALSE)
  if (!is.null(parentList)) {
    if (is.null(names(parentList)) || any(!nzchar(names(parentList)))) stop("`parentList` must be a named list.", call. = FALSE)
    parentDf <- do.call(rbind, lapply(names(parentList), function(pn) {
      P <- intersect(parentList[[pn]], types)
      g <- expand.grid(from = types, to = P, stringsAsFactors = FALSE)
      if (!nrow(g)) return(NULL)
      g$parent <- rep(list(P), nrow(g)); g$parent_name <- pn; g }))
    if (is.null(parentDf)) stop("no cell type of `parentList` is in the data.", call. = FALSE)
  }
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
