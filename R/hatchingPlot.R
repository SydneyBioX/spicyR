#' @export
#' @import ggplot2
hatchingPlot <- function(data, image = NULL, window = "concave", line.spacing = 21, 
    nbp = 250, window.length = 0) {
    if (!is(data, "segmentedCells") | is.null(data$region)) 
        stop("Please provide a segmentedCells object with region information please.")
    
    if (is.null(image)) {
        df <- region(data[1, ], annot = TRUE)
        p <- ggplot(df, aes(x = x, y = y, colour = cellType)) + geom_point() + geom_hatching(aes(region = region), 
            show.legend = TRUE, window = "concave", line.spacing = 21, nbp = 250, 
            window.length = 0)
        q <- p + theme_minimal() + scale_region()
        return(q)
    }
    
    if (!is.null(image)) {
        if (any(!image %in% rownames(data))) 
            stop("Some of the images are not in your segmentedCells object")
        df <- region(cellExp, image = image, annot = TRUE)
        p <- ggplot(df, aes(x = x, y = y, colour = cellType)) + geom_point() + facet_wrap(~imageID) + 
            geom_hatching(aes(region = region), show.legend = TRUE, window = "concave", 
                line.spacing = 21, nbp = 250, window.length = 0)
        q <- p + theme_minimal() + scale_region()
        q
    }
}
