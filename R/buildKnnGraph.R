#' Build a k-nearest-neighbour graph between cells
#'
#' Finds the `k` nearest neighbours of every cell within its own image and stores them as a
#' `colPair` of `cells`, ready for [convPairs()]. Distances are Euclidean, a cell is not its own
#' neighbour, and ties are broken as in `spatstat.geom::nnwhich()`. The graph is the one
#' `imcRtools::buildSpatialGraph(type = "knn")` builds, computed by the C++ core that [spicy()]
#' uses for `k`.
#'
#' @param cells A `SingleCellExperiment` or `SpatialExperiment`.
#' @param k The number of neighbours of each cell. Cells in images with at most `k` cells get no
#'   neighbours.
#' @param imageID The `colData` column holding the image of each cell.
#' @param spatialCoords The `colData` columns holding the coordinates. A `SpatialExperiment`'s
#'   own `spatialCoords()` are used unless these name `colData` columns.
#' @param name The name of the `colPair` the graph is stored in.
#' @param cores The number of threads; images are processed in parallel.
#'
#' @return `cells`, with the graph in `colPair(cells, name)`: one edge from each cell to each of
#'   its neighbours.
#' @export
#'
#' @examples
#' data("diabetesData")
#' diabetesData <- buildKnnGraph(diabetesData, k = 10)
#' head(SingleCellExperiment::colPair(diabetesData, "knn_interaction_graph"))
buildKnnGraph <- function(cells,
                          k = 20,
                          imageID = "imageID",
                          spatialCoords = c("x", "y"),
                          name = "knn_interaction_graph",
                          cores = 1) {
  .need("SingleCellExperiment", "for buildKnnGraph()")
  if (!.is_class(cells, "SingleCellExperiment")) {
    stop("`cells` must be a SingleCellExperiment or SpatialExperiment.", call. = FALSE)
  }
  if (length(k) != 1L || !is.numeric(k) || k < 1 || k != round(k)) {
    stop("`k` must be a positive integer.", call. = FALSE)
  }
  cd <- SummarizedExperiment::colData(cells)
  if (!imageID %in% colnames(cd)) stop("`imageID` (", imageID, ") is not a colData column.", call. = FALSE)
  if (all(spatialCoords %in% colnames(cd))) {
    xy <- cbind(cd[[spatialCoords[1]]], cd[[spatialCoords[2]]])
  } else if (.is_class(cells, "SpatialExperiment")) {
    xy <- SpatialExperiment::spatialCoords(cells)[, 1:2, drop = FALSE]
  } else {
    stop("`spatialCoords` (", paste(spatialCoords, collapse = ", "), ") are not colData columns.", call. = FALSE)
  }

  # the core wants the cells grouped by image
  img <- match(as.character(cd[[imageID]]), unique(as.character(cd[[imageID]])))
  ord <- order(img, method = "radix")
  offsets <- as.integer(c(0L, cumsum(tabulate(img, nbins = max(img, 0L)))))
  nb <- knn_rows(as.numeric(xy[ord, 1]), as.numeric(xy[ord, 2]), offsets, as.integer(k), as.integer(cores))

  # edges from each cell to its neighbours, as indices into the original cell order
  from <- rep(ord, each = k)
  keep <- nb >= 0L
  SingleCellExperiment::colPair(cells, name) <- S4Vectors::SelfHits(
    from[keep], ord[nb[keep] + 1L], nnode = ncol(cells)
  )
  cells
}
