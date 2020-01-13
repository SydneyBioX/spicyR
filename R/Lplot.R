Lplot <- function(x, from, to, imageID=1, rinterval=c(0,100)) {
    #TODO: Produce error messages
    # e.g. Missing entries in 'cells', 'from' or 'to' not present in cellTypes
    #TODO: make plot with ggplot
    
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
    
    plot(spatstat::Lcross(pppCell, from=from, to=to), xlim=rinterval)
}