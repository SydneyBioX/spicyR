#' Plots an image with specified from and to cell types.
#' 
#' @param cells A SummarizedExperiment object.
#' @param imageToPlot The ID of the image to be plotted.
#' @param from The "from" cell type.
#' @param to The "to" cell type.
#' @param imageID The name of the imageID column in the SummarizedExperiment object.
#' @param cellType The name of the cellType column in the SummarizedExperiment object. 
#' @param spatialCoords The names of the spatialCoords column if using a SingleCellExperiment.
#' 
#' @return A ggplot object.
#' 
#' @examples
#' data("diabetesData")
#' plotImage(diabetesData, "A09", from = "acinar", to = "alpha")
#' 
#' @export
#' @import ggplot2
#' @importFrom stats density
plotImage = function(cells, 
                     imageToPlot, 
                     from,
                     to,
                     imageID = "imageID", 
                     cellType = "cellType",
                     spatialCoords = c("x", "y")) {
  
  if (!.is_class(cells, "SummarizedExperiment")) {
    stop(paste("Please provide a SummarizedExperiment object as input."))
  }
  cd <- .col_data(cells)
  
  if (!imageID %in% colnames(cd)) {
    stop(paste0(imageID, " not found in colData."))
  }
  
  if (!imageToPlot %in% unique(cd[[imageID]])) {
    stop(paste0("imageToPlot not found in ", imageID, " column."))
  }
  
  if (!cellType %in% colnames(cd)) {
    stop(paste0(cellType, " not found in colData."))
  }
  
  if (length(spatialCoords) != 2) {
    stop(paste("Please provide x and y coordinates columns."))
  }

  # coordinates: spatialCoords() of a SpatialExperiment, else the colData columns
  if (.is_class(cells, "SpatialExperiment")) {
    coords <- as.data.frame(SpatialExperiment::spatialCoords(cells))
    if (!all(spatialCoords %in% colnames(coords))) coords <- coords[, 1:2, drop = FALSE]
    else coords <- coords[, spatialCoords, drop = FALSE]
  } else {
    if (!all(spatialCoords %in% colnames(cd))) {
      stop(paste0(spatialCoords, " not found in colData. "))
    }
    coords <- cd[, spatialCoords, drop = FALSE]
  }
  
  if (!all(c(from, to) %in% unique(cd[[cellType]]))) {
    stop("from and/or to cell types not found in data.")
  }
  
  # filter for specific image
  keep <- cd[[imageID]] == imageToPlot
  cData = data.frame(x = coords[keep, 1],
                     y = coords[keep, 2],
                     cellType = as.character(cd[[cellType]][keep]))
  cData$cellTypeNew <- ifelse(cData$cellType %in% c(from, to), cData$cellType, "Other")
  
  
  pal = setNames(c("#d6b11c", "#850f07"), c(from, to))
  
  ggplot() +
    stat_density_2d(data = cData, aes(x = .data$x, y = .data$y, fill = after_stat(density)), 
                    geom = "raster", 
                    contour = FALSE) +
    geom_point(data = cData[cData$cellTypeNew != "Other", , drop = FALSE],
               aes(x = .data$x, y = .data$y, colour = .data$cellTypeNew), size = 1) +
    scale_color_manual(values = pal) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    theme_classic() +
    labs(title = paste0(imageID, ": ", imageToPlot),
         color = cellType)
}
