
lisa <- function(cells, Rs = c(2, 20, 50, 100), BPPARAM = SerialParam(), window = "square", window.length = 20) {
    
 # if(is.data.frame(cells)){   
 #    if (is.null(cells$cellID)) {
 #        cat("Creating cellID as it doesn't exist")
 #        cells$cellID <- paste("cell", seq_len(nrow(cells)), sep = "_")
 #    }
 #   
 #   if (is.null(cells$cellType)) {
 #     stop("I need celltype")
 #    }
 #    
 #   if(is.null(cells$uniqueCellID)){
 #     cells$uniqueCellID <- paste("cell", seq_len(nrow(cells)), sep = "_")
 #   }
 #    if (length(cells$uniqueCellID) != nrow(cells) )
 #        stop("The number of rows in cells does not equal the number of uniqueCellIDs")
 #    
 #  if (is.null(cells$imageID)) {
 #    cat("There is no imageID. I'll assume this is only one image and create an arbitrary imageID")
 #    cells$imageID <- "image1"
 #  }
 #  
 #  location <- cells[,c('uniqueCellID','cellID','imageID','x','y', 'cellType')]
 # }
  

    
    
    if(is.null(uniqueCellID(cells))) cells <- makeUniqueCellID(cells)
    location <- location(cells)
    
    curveList <- BiocParallel::bplapply(location, generateCurves, Rs = Rs, BPPARAM = BPPARAM, window = window, window.length = window.length, ...)
                                                
                                                #Rs = Rs, BPPARAM = BPPARAM, window = window, window.length = window.length)
    
    curves <- do.call("rbind", curveList)
    curves <- curves[cells$cellID, ]
    return(curves)
}


makeWindow <- function(data, window = FALSE, window.length = 21){
  ow <- spatstat::owin(xrange = range(data$x), yrange = range(data$y))
  
  if (window == "convex") {
    ch <- grDevices::chull(data[, c("x", "y")])
    poly <- data[, c("x", "y")][rev(ch), ]
    colnames(poly) <- c("x", "y")
    ow <- spatstat::owin(xrange = range(data$x), yrange = range(data$y), 
                         poly = poly)
  }
  if (window == "concave") {
    ch <- concaveman::concaveman(as.matrix(data[, c("x", "y")]), length = 20)
    poly <- as.data.frame(ch[nrow(ch):1, ])
    ch <- concaveman::concaveman(do.call("rbind", lapply(as.list(as.data.frame(t(poly))), 
                                                         function(x) cbind(rnorm(1000, x[1]), rnorm(1000, x[2])))), length = window.length)
    poly <- as.data.frame(ch[nrow(ch):1, ])
    colnames(poly) <- c("x", "y")
    ow <- spatstat::owin(xrange = range(poly$x), yrange = range(poly$y), 
                         poly = poly)
  }
  
}



generateCurves <- function(data, Rs, window, window.length, ...) {
  
  ow <- makeWindow
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
      
      return(locR)
    })
    return(do.call("rbind", locI))
  })
  return(do.call("cbind", locIJ))
}



