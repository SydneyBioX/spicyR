#' A basic plot for SegmentedCells object
#'
#' This function generates a basic x-y plot of the location coordinates and cellType data.
#'
#' @section usage:
#' `plot(x, imageID = NULL)`
#'
#' @param x A SegmentedCells object.
#' @param imageID The image that should be plotted.
#'
#' @return A ggplot object.
#'
#' @examples
#' ### Something that resembles cellProfiler data
#'
#' set.seed(51773)
#'
#' n = 10
#'
#' cells <- data.frame(row.names = seq_len(n))
#' cells$ObjectNumber <- seq_len(n)
#' cells$ImageNumber <- rep(1:2,c(n/2,n/2))
#' cells$AreaShape_Center_X <- runif(n)
#' cells$AreaShape_Center_Y <- runif(n)
#' cells$AreaShape_round <- rexp(n)
#' cells$AreaShape_diameter <- rexp(n, 2)
#' cells$Intensity_Mean_CD8 <- rexp(n, 10)
#' cells$Intensity_Mean_CD4 <- rexp(n, 10)
#'
#' cellExp <- SegmentedCells(cells, cellProfiler = TRUE)
#'
#' ### Cluster cell types
#' markers <- cellMarks(cellExp)
#' kM <- kmeans(markers,2)
#' cellType(cellExp) <- paste('cluster',kM$cluster, sep = '')
#'
#' #plot(cellExp, imageID=1)
#'
#' @name plot-SegmentedCells
#' @rdname plot-SegmentedCells
#' @aliases plot
#' @importFrom ggplot2 ggplot aes geom_point theme_classic labs
#' @importFrom rlang .data
NULL

plot.SegmentedCells <- function(cellData, imageID = NULL) {
    if (is.null(imageID)) {
        imageID <- imageID(cellData)[1]
    }
    
    loc <- as.data.frame(cellSummary(cellData, imageID = imageID))
    if (is.na(loc$cellType[1])) {
        ggplot(loc, aes(x = .data$x, y = .data$y)) + geom_point() + theme_classic() + labs(x = "x", y = "y")
    } else {
        ggplot(loc, aes(
            x = .data$x,
            y = .data$y,
            colour = .data$cellType
        )) + geom_point() +
            theme_classic() + labs(x = "x",
                                   y = "y",
                                   colour = "cell-type")
    }
}


if (!isGeneric("plot"))
    setGeneric("plot", function(x, ...)
        standardGeneric("plot"))

setMethod("plot", signature(x = "SegmentedCells"), function(x, ...) {
    plot.SegmentedCells(x, ...)
})


