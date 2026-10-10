#' Converts colPairs object into an abundance matrix based on number of nearby
#' interactions for every cell type.
#'
#' @param cells
#'   A SingleCellExperiment that contains objects in the colPairs slot.
#' @param colPair
#'   The name of the object in the colPairs slot for which the dataframe is
#'   constructed from.
#' @param imageID The image ID if using SingleCellExperiment.
#' @param cellType The cell type if using SingleCellExperiment.
#'
#' @return Matrix of abundances
#' @export
#'
#' @examples
#' data("diabetesData")
#' images <- c("A09", "A11", "A16", "A17")
#' diabetesData <- diabetesData[
#'   , SummarizedExperiment::colData(diabetesData)$imageID %in% images
#' ]
#'
#' diabetesData <- buildKnnGraph(diabetesData, k = 20)
#'
#' pairAbundances <- convPairs(diabetesData,
#'   colPair = "knn_interaction_graph"
#' )
#'
#' @export
convPairs <- function(cells,
                      colPair,
                      imageID = "imageID",
                      cellType = "cellType") {
  .need("SingleCellExperiment", "for convPairs()")
  hits <- SingleCellExperiment::colPair(cells, colPair)
  from <- S4Vectors::from(hits)
  to <- S4Vectors::to(hits)
  cd <- as.data.frame(SummarizedExperiment::colData(cells))
  img <- cd[[imageID]]
  type <- cd[[cellType]]

  # order of groups as dplyr::group_by() sorts them: factor levels, else C-locale order
  sortKey <- function(v) {
    if (is.factor(v)) as.integer(v) else match(v, sort(unique(v), method = "radix"))
  }

  # number of `to` cells of each type next to `from` cells of each type, per image
  edges <- data.frame(imageID = img[from], cellType_from = type[from], cellType_to = type[to])
  key <- paste(sortKey(edges$imageID), sortKey(edges$cellType_from), sortKey(edges$cellType_to), sep = "\r")
  first <- !duplicated(key)
  groups <- edges[first, , drop = FALSE]
  groups$n_close <- as.vector(table(factor(key, levels = key[first])))
  groups <- groups[order(sortKey(groups$imageID), sortKey(groups$cellType_from), sortKey(groups$cellType_to),
                         method = "radix"), , drop = FALSE]

  # divided by the number of `from` cells in the image
  nType <- table(paste(img, type, sep = "\r"))
  groups$association <- groups$n_close /
    as.vector(nType[paste(groups$imageID, groups$cellType_from, sep = "\r")])
  groups$test <- paste(groups$cellType_from, groups$cellType_to, sep = "__")

  # wide: one row per image, one column per pair, absent pairs 0
  rows <- unique(as.character(groups$imageID))
  tests <- unique(groups$test)
  m <- matrix(0, length(rows), length(tests), dimnames = list(rows, tests))
  m[cbind(match(as.character(groups$imageID), rows), match(groups$test, tests))] <- groups$association
  all_pairs <- as.data.frame(m, check.names = FALSE)

  # Hot fix for spicy input when no cell type interactions exist for a pairwise
  # relation.
  vector <- unique(type)

  pairwise_vector <- c()

  for (i in vector) {
    for (j in vector) {
      pairwise_vector <- c(pairwise_vector, paste(i, j, sep = "__"))
    }
  }

  tmp <- setdiff(pairwise_vector, colnames(all_pairs))
  df <- data.frame(matrix(0, nrow = nrow(all_pairs), ncol = length(tmp)))
  colnames(df) <- tmp

  all_pairs <- cbind(all_pairs, df)

  return(all_pairs)
}
