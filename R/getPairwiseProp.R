#' Observed/expected proportion ratio for a fixed-k neighbourhood
#'
#' Computes a per-image, per-cell-type-pair spatial localisation statistic as a
#' drop-in alternative to the fixed-radius L-function statistic returned by
#' \code{\link{getPairwise}}. Instead of a search radius \code{r}, each reference
#' cell's neighbourhood is its \code{k} nearest neighbours (Euclidean distance,
#' any cell type, the cell itself excluded).
#'
#' For an ordered pair (\code{from} = reference type A, \code{to} = target type
#' B), the statistic for image \eqn{j} is the plain observed-over-expected ratio
#' \deqn{v_j = \hat\pi_j / p_{0,j},}
#' where \eqn{\hat\pi_j} is the mean over reference cells of the fraction of their
#' \code{k} nearest neighbours that are type B, and \eqn{p_{0,j}} is the image's
#' background proportion of type-B cells (B count / total cell count).
#'
#' \eqn{v_j = 1} is the null (complete spatial randomness); \eqn{v_j > 1}
#' indicates attraction/enrichment and \eqn{0 \le v_j < 1} depletion. No
#' variance-stabilising transform (arcsine, log) is applied: \eqn{v_j} is bounded
#' below at 0 and right-skewed, a rougher fit to the linear model's Gaussian-noise
#' assumption than a transformed statistic, but chosen for direct
#' interpretability ("observed is 1.4x expected").
#'
#' The returned matrix is shaped identically to \code{getPairwise()}'s output
#' (images in rows, \code{from__to} cell-type pairs in columns, same order), so it
#' can be passed straight to \code{spicy()} as \code{alternateResult} to run the
#' standard spicyR weighting and linear (mixed) model pipeline on the ratio.
#'
#' \strong{Weights.} \code{spicy()}'s variance-weighting model (\code{weights =
#' TRUE}, the default) is calibrated to the L-function's numeric scale. It works
#' for this ratio when the dataset has enough cross-pair variance heterogeneity
#' (e.g. several cell types, some rare) to anchor the fit, but on a dataset with
#' few cell types or weak spatial structure it can degenerate and return an
#' all-\code{NA} p-value table. If that happens, refit with \code{weights =
#' FALSE}.
#'
#' No edge correction is applied: \eqn{\hat\pi_j} estimates neighbourhood
#' composition rather than a count within an area, so a boundary cell's truncated
#' neighbourhood carries the same null label distribution as an interior cell's.
#' \code{k} controls only the variance of \eqn{\hat\pi_j}, not the location of its
#' null value.
#'
#' @param cells A SingleCellExperiment, SpatialExperiment or data.frame.
#' @param imageID The name of the imageID column if a column name is provided.
#' @param cellType The name of the cellType column if a column name is provided.
#' @param spatialCoords The names of the spatial coordinates columns if provided.
#' @param k Number of nearest neighbours to search per reference cell.
#' @param from The reference cell types. Defaults to all cell types.
#' @param to The target cell types. Defaults to all cell types.
#' @param cores Number of cores or a BiocParallelParam object for parallel
#'     processing over images.
#' @param includeZeroCells
#'     If FALSE (default), image-pairs where the reference or target cell type is
#'     absent are returned as NA. If TRUE, image-pairs whose reference type is
#'     absent (but whose target type is present, so \eqn{p_{0,j} > 0}) are floored
#'     at \eqn{\hat\pi_j = 0}, i.e. \eqn{v_j = 0}. Pairs whose target type is
#'     entirely absent (\eqn{p_{0,j} = 0}) stay NA regardless, since the ratio is
#'     undefined.
#' @param BPPARAM A BiocParallelParam object. Overrides \code{cores}.
#'
#' @return A matrix with one row per image and one column per \code{from__to}
#'     cell-type pair, containing the observed/expected proportion ratio.
#'
#' @examples
#' data("diabetesData")
#' propAssoc <- getPairwiseProp(diabetesData, k = 15)
#'
#' # Use as the spicyR response, reusing the standard model:
#' \dontrun{
#' spicy(diabetesData,
#'   condition = "stage", subject = "case",
#'   alternateResult = propAssoc
#' )
#' }
#'
#' @export
#' @importFrom spatstat.geom nnwhich ppp owin
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam SerialParam
getPairwiseProp <- function(
    cells,
    imageID = "imageID",
    cellType = "cellType",
    spatialCoords = c("x", "y"),
    k = 15,
    from = NULL,
    to = NULL,
    cores = 1,
    includeZeroCells = FALSE,
    BPPARAM = NULL) {
  if (length(k) != 1 || !is.finite(k) || k < 1) {
    stop("`k` must be a single positive integer.")
  }
  k <- as.integer(k)

  if (is(cells, "SummarizedExperiment") || is.data.frame(cells)) {
    cells <- .format_data(cells, imageID, cellType, spatialCoords, FALSE)
  }

  if (is.null(BPPARAM)) {
    if (cores > 1 && .Platform$OS.type != "windows") {
      BPPARAM <- BiocParallel::MulticoreParam(workers = cores)
    } else if (cores > 1) {
      BPPARAM <- BiocParallel::SnowParam(workers = cores)
    } else {
      BPPARAM <- BiocParallel::SerialParam()
    }
  }

  cells2 <- getCellSummary(cells, bind = FALSE)

  # Match spicy()'s own label construction (order of first appearance), so the
  # returned matrix drops into spicy(..., alternateResult=) without reordering.
  allTypes <- as.character(unique(getCellType(cells)))
  if (is.null(from)) from <- allTypes
  if (is.null(to)) to <- allTypes

  pairwiseVals <- BiocParallel::bplapply(cells2,
    propPair,
    k = k,
    from = from,
    to = to,
    includeZeroCells = includeZeroCells,
    BPPARAM = BPPARAM
  )

  do.call("rbind", pairwiseVals)
}


#' @importFrom spatstat.geom nnwhich ppp owin
propPair <- function(data,
                     k = 15,
                     from = NULL,
                     to = NULL,
                     includeZeroCells = FALSE) {
  x <- as.numeric(data$x)
  y <- as.numeric(data$y)
  cellTypeLevels <- levels(data$cellType)
  types <- as.character(data$cellType)
  nAll <- length(types)

  if (is.null(from)) from <- cellTypeLevels
  if (is.null(to)) to <- cellTypeLevels

  m1 <- rep(from, times = length(to))
  m2 <- rep(to, each = length(from))
  labels <- paste(m1, m2, sep = "__")
  assoc <- rep(NA_real_, length(labels))
  names(assoc) <- labels

  if (nAll <= k) {
    return(assoc)
  }

  pp <- spatstat.geom::ppp(x, y,
    window = spatstat.geom::owin(range(x), range(y)),
    check = FALSE
  )
  nnIdx <- matrix(spatstat.geom::nnwhich(pp, k = seq_len(k)), ncol = k)
  nnType <- matrix(types[nnIdx], ncol = k)

  for (idx in seq_along(labels)) {
    A <- m1[idx]
    B <- m2[idx]
    idxFrom <- which(types == A)
    nB <- sum(types == B)
    p0 <- nB / nAll

    degenerate <- length(idxFrom) == 0 || nB == 0 || nB == nAll

    if (degenerate) {
      if (includeZeroCells && nB > 0) {
        assoc[idx] <- 0
      }
      next
    }

    piHat <- mean(rowSums(nnType[idxFrom, , drop = FALSE] == B)) / k
    assoc[idx] <- piHat / p0
  }

  assoc
}
