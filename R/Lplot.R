Lplot <- function(x, from, to, imageID=1, rinterval=c(0,100)) {
    
    if (is(x, "SegmentedCellExperiment")) {
        cells <- location(x)
        cells <- cells[cells$imageID==as.character(imageID),]
    } else if (is(x, "data.frame")) {
        cells <- x
        if (is.null(cells$imageID)) {
            cells$imageID <- rep(imageID, nrow(cells))
        }
    } else {
        stop("x is not a data.frame or SegmentedCellExperiment")
    }
    
    pppCell <- pppGenerate(cells)
    
    L <- spatstat::Lcross(pppCell, from=from, to=to)
    
    L.df <- data.frame(x = L$r, y = L$iso)
    
    g <- ggplot2::ggplot(L.df, ggplot2::aes(x=x)) +
             ggplot2::geom_line(ggplot2::aes(y=y), color = "red") +
             ggplot2::theme_classic() +
             ggplot2::xlab("Distance") + 
             ggplot2::ylab("L")
    
    g
}