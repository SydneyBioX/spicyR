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
#' @import dplyr 
#' @import ggplot2 
#' @import SummarizedExperiment
#' @import SpatialExperiment
plotImage = function(cells, 
                     imageToPlot, 
                     from,
                     to,
                     imageID = "imageID", 
                     cellType = "cellType",
                     spatialCoords = c("x", "y")) {
  
  if (!is(cells, "SummarizedExperiment")) {
    stop(paste("Please provide a SummarizedExperiment object as input."))
  }
  
  if (!imageID %in% colnames(colData(cells))) {
    stop(paste0(imageID, " not found in colData."))
  }
  
  if (!imageToPlot %in% unique(cells[[imageID]])) {
    stop(paste0("imageToPlot not found in ", imageID, " column."))
  }
  
  if (!cellType %in% colnames(colData(cells))) {
    stop(paste0(cellType, " not found in colData."))
  }
  
  if (class(cells) == "SingleCellExperiment") {
    if (!all(spatialCoords %in% colnames(colData(cells)))) {
      stop(paste0(spatialCoords, " not found in colData. "))
    }
  }
  
  if (length(spatialCoords) != 2) {
    stop(paste("Please provide x and y coordinates columns."))
  }
  
  if (!all(c(from, to) %in% unique(colData(cells)[[cellType]]))) {
    stop("from and/or to cell types not found in data.")
  }
  
  if (class(cells) == "SingleCellExperiment") {
    cells = SpatialExperiment(assays = assays(cells), 
                              colData = colData(cells), 
                              rowData = rowData(cells))
    spatialCoords(cells) = as.matrix(colData(cells)[, c(spatialCoords[1], spatialCoords[2])])
  }
  
  # filter for specific image
  subset = cells[, colData(cells)[[imageID]] == imageToPlot]
  
  coords = spatialCoords(subset) |> as.data.frame()
  cData = data.frame(x = coords[["x"]],
                     y = coords[["y"]],
                     cellType = subset[[cellType]])
  
  cData = data.frame(x = coords[["x"]],
                     y = coords[["y"]],
                     cellType = subset[[cellType]])
  
  cData$cellType = as.character(cData$cellType)
  
  cData[[cellType]] = as.character(cData[[cellType]])
  cData = cData |> mutate(cellTypeNew = 
                            ifelse(cellType %in% c(from, to), cellType, "Other"))
  
  
  pal = setNames(c("#d6b11c", "#850f07"), c(from, to))
  
  ggplot() +
    stat_density_2d(data = cData, aes(x = x, y = y, fill = after_stat(density)), 
                    geom = "raster", 
                    contour = FALSE) +
    geom_point(data = cData |> filter(cellTypeNew != "Other"),
               aes(x = x, y = y, colour = cellTypeNew), size = 1) +
    scale_color_manual(values = pal) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    theme_classic() +
    labs(title = paste0(imageID, ": ", imageToPlot),
         color = cellType)
}