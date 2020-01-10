#' Generate local indicators of spatial association
#'
#' @param cells A SegmentedCellExperiment or data frame that contains at least the variables x and y, giving the location of each cell, and cellType.
#' @param Rs A vector of the radii that the measures of association should be calculated.
#' @param BPPARAM A BiocParallelParam object.
#' @param window Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length A tuning parameter for controlling the level of concavity when estimating concave windows.
#' 
#' @examples
#' 1+1
#' @export
#' @rdname lisa
#' @importFrom methods is
#' @importFrom BiocParallel bplapply SerialParam
#' @import SegmentedCellExperiment S4Vectors BiocGenerics
lisa <- function(cells, Rs = NULL, BPPARAM = SerialParam(), window = "square", 
    window.length = NULL) {
    
    if (is.data.frame(cells)) {
        if (is.null(cells$cellID)) {
            cat("Creating cellID as it doesn't exist")
            cells$cellID <- paste("cell", seq_len(nrow(cells)), sep = "_")
        }
        
        if (is.null(cells$cellType)) {
            stop("I need celltype")
        }
        
        if (is.null(cells$x) | is.null(cells$y)) {
            stop("I need x and y coordinates")
        }
        
        if (is.null(cells$imageCellID)) {
            cells$imageCellID <- paste("cell", seq_len(nrow(cells)), sep = "_")
        }
        if (length(unique(cells$imageCellID)) != nrow(cells)) 
            stop("The number of rows in cells does not equal the number of uniqueCellIDs")
        
        if (is.null(cells$imageID)) {
            cat("There is no imageID. I'll assume this is only one image and create an arbitrary imageID")
            cells$imageID <- "image1"
        }
        
        location <- split(S4Vectors::DataFrame(cells[, c("imageID", "imageCellID", "cellID", "x", "y", 
            "cellType")]), cells$imageID)
    }
    
    if (is(cells, "SegmentedCellExperiment")) {
        location <- location(cells, bind = FALSE)
    }
  
   if (is.null(Rs)){
     loc = do.call('rbind', location)
     range <- max(loc$x) - min(loc$x)
     maxR <- range/5
     Rs = seq(from = maxR / 20, maxR, length.out = 20)
   } 
    
    curveList <- BiocParallel::bplapply(location, generateCurves, Rs = Rs, BPPARAM = BPPARAM, 
        window = window, window.length = window.length)
    
    curves <- do.call("rbind", curveList)
    curves <- curves[location(cells)$cellID, ]
    return(curves)
}


#' @importFrom grDevices chull
#' @importFrom spatstat owin
#' @importFrom concaveman concaveman
makeWindow <- function(data, window = "square", window.length = 21) {
  data = data.frame(data)
    ow <- spatstat::owin(xrange = range(data$x), yrange = range(data$y))
    
    if (window == "convex") {
        ch <- grDevices::chull(data[, c("x", "y")])
        poly <- data[, c("x", "y")][rev(ch), ]
        colnames(poly) <- c("x", "y")
        ow <- spatstat::owin(xrange = range(data$x), yrange = range(data$y), poly = poly)
    }
    if (window == "concave") {
        ch <- concaveman::concaveman(as.matrix(data[, c("x", "y")]))
        poly <- as.data.frame(ch[nrow(ch):1, ])
        range1 <- max(poly[, 1]) - min(poly[, 1])
        range2 <- max(poly[, 2]) - min(poly[, 2])
        ch <- concaveman::concaveman(do.call("rbind", lapply(as.list(as.data.frame(t(poly))), 
            function(x) cbind(rnorm(1000, x[1], range1/1000), rnorm(1000, x[2], range2/1000)))), 
            length_threshold = window.length)
        poly <- as.data.frame(ch[nrow(ch):1, ])
        colnames(poly) <- c("x", "y")
        ow <- spatstat::owin(xrange = range(poly$x), yrange = range(poly$y), poly = poly)
    }
    ow
}


#' @importFrom spatstat ppp
#' @importFrom spatstat localLcross
generateCurves <- function(data, Rs, window, window.length, ...) {
    
    ow <- makeWindow(data, window, window.length)
    p1 <- spatstat::ppp(x = data$x, y = data$y, window = ow, marks = as.factor(data$cellType))
    
    
    locIJ <- lapply(as.list(levels(p1$marks)), function(j) {
        
        locI <- lapply(as.list(levels(p1$marks)), function(i) {
            
            iID <- data$cellID[p1$marks == i]
            jID <- data$cellID[p1$marks == j]
            locR <- matrix(NA, length(iID), length(Rs))
            rownames(locR) <- iID
            
            if (length(jID) > 1 & length(iID) > 1) {
                
                localL <- spatstat::localLcross(p1, from = i, to = j, verbose = FALSE)
                ur <- sapply(Rs, function(x) which.min(abs(localL$r - x))[1])
                locR <- t(apply(as.matrix(localL)[, grep("iso", colnames(localL))], 
                  2, function(x) (x - localL$theo)[ur]))
                rownames(locR) <- iID
            }
            colnames(locR) <- paste(j, round(Rs, 2), sep = "_")
            
            locR
        })
        do.call("rbind", locI)
    })
    do.call("cbind", locIJ)
}




