#' Plot one image, showing its cell types or the `from` and `to` cells of a pair
#'
#' Without `from` and `to`, every cell is shown, coloured by its cell type. With them, the density of all cells is
#' shown in blue, with the `from` and `to` cells on top. With `r`, a circle of radius `r` is drawn around each `to`
#' cell: the `from` cells inside the circles are those that spicyR counts.
#'
#' @param cells A SummarizedExperiment object.
#' @param imageToPlot The ID of the image to be plotted.
#' @param from The "from" cell type. Optional, with `to`.
#' @param to The "to" cell type. Optional, with `from`.
#' @param imageID The name of the imageID column in the SummarizedExperiment object.
#' @param cellType The name of the cellType column in the SummarizedExperiment object. 
#' @param spatialCoords The names of the spatialCoords column if using a SingleCellExperiment.
#' @param r Optional radius: draw a circle of this radius around each `to` cell.
#' 
#' @return A ggplot object.
#' 
#' @examples
#' data("diabetesData")
#' plotImage(diabetesData, "A09")
#' plotImage(diabetesData, "A09", from = "acinar", to = "alpha")
#' plotImage(diabetesData, "A09", from = "acinar", to = "alpha", r = 50)
#' 
#' @export
#' @import ggplot2
#' @importFrom stats density
#' @concept plot
plotImage = function(cells, 
                     imageToPlot, 
                     from = NULL,
                     to = NULL,
                     imageID = "imageID", 
                     cellType = "cellType",
                     spatialCoords = c("x", "y"),
                     r = NULL) {
  
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
  
  if (is.null(from) != is.null(to)) stop("give both `from` and `to`, or neither.")
  if (is.null(from) && !is.null(r)) stop("`r` draws circles around the `to` cells: give `from` and `to`.")
  if (!all(c(from, to) %in% unique(cd[[cellType]]))) {
    stop("from and/or to cell types not found in data.")
  }
  
  # filter for specific image
  keep <- cd[[imageID]] == imageToPlot
  cData = data.frame(x = coords[keep, 1],
                     y = coords[keep, 2],
                     cellType = as.character(cd[[cellType]][keep]))
  if (is.null(from)) {   # every cell, coloured by its cell type: a type has the same colour in every image
    types <- if (is.factor(cd[[cellType]])) levels(cd[[cellType]]) else sort(unique(as.character(cd[[cellType]])))
    hues <- seq(15, 375, length.out = length(types) + 1)[seq_along(types)]   # ggplot2's default hue palette
    pal <- setNames(grDevices::hcl(h = hues, c = 100, l = 65), types)
    return(ggplot(cData, aes(x = .data$x, y = .data$y, colour = .data$cellType)) +
             geom_point(size = 0.5) +
             scale_colour_manual(values = pal, breaks = intersect(types, cData$cellType)) +
             coord_equal() +
             theme_classic() +
             labs(title = paste0(imageID, ": ", imageToPlot), colour = cellType) +
             guides(colour = guide_legend(override.aes = list(size = 2))))
  }
  cData$cellTypeNew <- ifelse(cData$cellType %in% c(from, to), cData$cellType, "Other")
  
  
  pal = setNames(c("#d6b11c", "#850f07"), c(from, to))
  
  circles <- NULL
  if (!is.null(r)) {
    f <- cData[cData$cellType == to, , drop = FALSE]
    a <- seq(0, 2 * pi, length.out = 41)
    circles <- data.frame(x = rep(f$x, each = 41) + r * cos(a), y = rep(f$y, each = 41) + r * sin(a),
                          id = rep(seq_len(nrow(f)), each = 41))
  }
  p <- ggplot() +
    stat_density_2d(data = cData, aes(x = .data$x, y = .data$y, fill = after_stat(density)), 
                    geom = "raster", 
                    contour = FALSE) +
    geom_point(data = cData[cData$cellTypeNew != "Other", , drop = FALSE],
               aes(x = .data$x, y = .data$y, colour = .data$cellTypeNew), size = 1) +
    scale_color_manual(values = pal) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    coord_equal() +
    theme_classic() +
    labs(title = paste0(imageID, ": ", imageToPlot),
         color = cellType)
  if (!is.null(circles))
    p <- p + geom_path(data = circles, aes(x = .data$x, y = .data$y, group = .data$id), colour = "#850f07",
                       linewidth = 0.3, alpha = 0.7)
  p
}
