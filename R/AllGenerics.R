
### Get morphology information

#' @export
setGeneric("region", function(x, image = NULL, annot = FALSE) standardGeneric("region"))
setMethod("region", "SegmentedCellExperiment", function(x, image = NULL, annot = TRUE) {
    if (!is.null(image)) {
        x <- x[image, ]
    }
    if (is.null(x$region)) 
        stop("There is no region information in your SegmentedCellExperiment yet")
    if (annot) 
        return(data.frame(location(x), region = BiocGenerics::do.call("rbind", x$region)))
    
    BiocGenerics::do.call("rbind", x$region)
})


#' @export
#' @importFrom S4Vectors DataFrame split
setGeneric("region<-", function(x, value, image = NULL) standardGeneric("region<-"))
setMethod("region<-", "SegmentedCellExperiment", function(x, value, image = NULL) {
    if (is.null(image)) 
        image <- rownames(x)
    if (length(value) == length(imageID(x, image))) {
        if (is.null(x$region)) {
            x <- DataFrame(x)
            x$region <- S4Vectors::split(DataFrame(region = value), rep(rownames(x), 
                unlist(lapply(x$location, nrow))))
            x <- new("SegmentedCellExperiment", x)
        }
        if (!is.null(x$region)) 
            x[image, ]@listData$region <- S4Vectors::split(DataFrame(region = value), 
                rep(rownames(x), unlist(lapply(x$location, nrow))))
    }
    x
})


