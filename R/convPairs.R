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
#' data("diabetesData_SCE")
#' 
#' diabetesData_SPE <- SpatialExperiment::SpatialExperiment(diabetesData_SCE,
#'   colData = SingleCellExperiment::colData(diabetesData_SCE))
#' SpatialExperiment::spatialCoords(diabetesData_SPE) <- data.frame(
#'   SingleCellExperiment::colData(diabetesData_SPE)$x, 
#'   SingleCellExperiment::colData(diabetesData_SPE)$y) |> 
#'     as.matrix()
#'     
#' SpatialExperiment::spatialCoordsNames(diabetesData_SPE) <- c("x", "y")
#' 
#' diabetesData_SPE <- imcRtools::buildSpatialGraph(diabetesData_SPE, 
#'   img_id = "imageID", 
#'   type = "knn", 
#'   k = 20, 
#'   coords = c("x", "y"))
#' 
#' pairAbundances <- convPairs(diabetesData_SPE,
#'   colPair = "knn_interaction_graph")
#' 
#' @export
#' @importFrom SingleCellExperiment colPairs colData
#' @importFrom tibble rownames_to_column column_to_rownames
convPairs <- function(cells,
                      colPair,
                      cellType = "cellType",
                      imageID = "imageID") {
  all_pairs <- SingleCellExperiment::colPair(cells, colPair) |>
    dplyr::as_tibble() |>
    # join the `from` cellType
    dplyr::left_join(
      SummarizedExperiment::colData(cells) |>
        dplyr::as_tibble() |>
        dplyr::select(imageID, cellType) |>
        tibble::rownames_to_column() |>
        dplyr::mutate(rowname = as.integer(rowname)),
      by = c("from" = "rowname")
    ) |>
    dplyr::rename(cellType_from = cellType) |>
    # join the `to` cellType
    dplyr::left_join(
      SummarizedExperiment::colData(cells) |>
        dplyr::as_tibble() |>
        dplyr::select(cellType) |>
        tibble::rownames_to_column() |>
        dplyr::mutate(rowname = as.integer(rowname)),
      by = c("to" = "rowname")
    ) |>
    dplyr::rename(cellType_to = cellType) |>
    # count the the number of `to` cellType associated with `from` cellTypes
    dplyr::mutate(one = 1L) |>
    dplyr::group_by(imageID, cellType_from, cellType_to) |>
    dplyr::summarise(n_close = sum(one), .groups = "drop") |>
    # join the total number of each cellType (within sample) into the dataframe
    dplyr::left_join(
      SummarizedExperiment::colData(cells) |>
        dplyr::as_tibble() |>
        dplyr::select(imageID, cellType) |>
        dplyr::group_by(imageID, cellType) |>
        dplyr::count(),
      by = c("cellType_from" = cellType, "imageID" = imageID)
    ) |>
    # calculate the association
    dplyr::mutate(association = n_close / n) |>
    dplyr::select(-n_close, -n) |>
    # wrangle the data into the correct format for spicy
    dplyr::mutate(
      test = paste(cellType_from, cellType_to, sep = "__")
    ) |>
    dplyr::select(imageID, test, association) |>
    tidyr::pivot_wider(
      names_from = test, values_from = association, values_fill = 0
    ) |>
    tibble::column_to_rownames(imageID)
  
  
  # Hot fix for spicy input when no cell type interactions exist for a pairwise 
  # relation.
  vector <- cells$cellType |> unique()
  
  pairwise_vector <- c()
  
  for (i in vector) {
    for (j in vector) {
      pairwise_vector <- c(pairwise_vector, paste(i, j, sep = "__"))
    }
  }
  
  tmp <- dplyr::setdiff(pairwise_vector, colnames(all_pairs))
  df <- data.frame(matrix(0, nrow = nrow(all_pairs), ncol = length(tmp)))
  colnames(df) <- tmp
  
  all_pairs <- cbind(all_pairs, df)
  
  return(all_pairs)
}