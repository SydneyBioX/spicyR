#' Generate local indicators of spatial association
#'
#' @param cells A SegmentedCells or data frame that contains at least the 
#' variables x and y, giving the  coordinates of each cell, and cellType.
#' @param Rs A vector of the radii that the measures of association should be calculated.
#' @param BPPARAM A BiocParallelParam object.
#' @param window Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length A tuning parameter for controlling the level of concavity 
#' when estimating concave windows.
#' @param whichParallel Should the function use parallization on the imageID or 
#' the cellType.
#' @param sigma A numeric variable used for scaling when filting inhomogeneous L-curves.
#' @param fast A logical describing whether to use a fast approximation of the 
#' inhomogeneous local L-curves.
#'
#' @return A matrix of LISA curves
#'
#' @examples
#' # Read in data as a SegmentedCells objects
#' isletFile <- system.file("extdata","isletCells.txt.gz", package = "spicyR")
#' cells <- read.table(isletFile, header=TRUE)
#' cellExp <- SegmentedCells(cells, cellProfiler = TRUE)
#'
#' # Cluster cell types
#' markers <- cellMarks(cellExp)
#' kM <- kmeans(markers,8)
#' cellType(cellExp) <- paste('cluster',kM$cluster, sep = '')
#'
#' # Generate LISA
#' lisaCurves <- lisa(cellExp)
#'
#' # Cluster the LISA curves
#' kM <- kmeans(lisaCurves,2)
#' region(cellExp) <- paste('region',kM$cluster,sep = '_')
#'
#' @export
#' @rdname lisa
#' @importFrom methods is
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocGenerics do.call rbind
lisa <-
  function(cells,
           Rs = NULL,
           BPPARAM = BiocParallel::SerialParam(),
           window = "convex",
           window.length = NULL,
           whichParallel = 'imageID',
           sigma = NULL,
           fast = TRUE) {
    if (is.data.frame(cells)) {
      if (is.null(cells$cellID)) {
        message("Creating cellID as it doesn't exist")
        cells$cellID <-
          paste("cell", seq_len(nrow(cells)), sep = "_")
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
        message(
          "There is no imageID. I'll assume this is only one image and create an arbitrary imageID"
        )
        cells$imageID <- "image1"
      }
      
      cellSummary <-
        split(S4Vectors::DataFrame(cells[, c("imageID",
                                             "imageCellID",
                                             "cellID",
                                             "x",
                                             "y",
                                             "cellType")]), cells$imageID)
    }
    
    if (is(cells, "SegmentedCells")) {
      cellSummary <- cellSummary(cells, bind = FALSE)
    }
    
    if (is.null(Rs)) {
      # loc = do.call('rbind', cellSummary)
      # range <- max(loc$x) - min(loc$x)
      # maxR <- range / 5
      # Rs = seq(from = maxR / 20, maxR, length.out = 20)
      Rs = c(20, 50, 100, 200)
    }
    
    BPimage = BPcellType = BiocParallel::SerialParam()
    if (whichParallel == 'imageID')
      BPimage <- BPPARAM
    if (whichParallel == 'cellType')
      BPcellType <- BPPARAM
    
    if(!fast){
      message("Generating local L-curves. ")
      if(identical(BPimage, BPcellType)) 
        message("You might like to consider setting BPPARAM to run the calculations in parallel.")
    curveList <-
      BiocParallel::bplapply(
        cellSummary,
        generateCurves,
        Rs = Rs,
        window = window,
        window.length = window.length,
        BPcellType = BPcellType,
        BPPARAM = BPimage,
        sigma = sigma
      )
    }
    
    if(fast){
      
      message("Generating local L-curves. If you run out of memory, try 'fast = FALSE'.")
      
      curveList <-
        BiocParallel::bplapply(
          cellSummary,
          inhomLocalL,
          Rs = Rs,
          window = window,
          window.length = window.length,
          BPPARAM = BPimage,
          sigma = sigma
        )
    }
    
    curves <- do.call("rbind", curveList)
    curves <- curves[as.character(cellSummary(cells)$cellID), ]
    return(curves)
  }



#' @importFrom spatstat owin convexhull
#' @importFrom concaveman concaveman
makeWindow <-
  function(data,
           window = "square",
           window.length = NULL) {
    data = data.frame(data)
    ow <-
      spatstat::owin(xrange = range(data$x), yrange = range(data$y))
    
    if (window == "convex") {
      p <- ppp(data$x, data$y, ow)
      ow <- spatstat::convexhull(p)
      
    }
    if (window == "concave") {
      message("Concave windows are temperamental. Try choosing values of window.length > and < 1 if you have problems.")
      if(is.null(window.length)){
        window.length <- (max(data$x) - min(data$x))/20
      }else{
        window.length <- (max(data$x) - min(data$x))/20 * window.length
      }
      dist <- (max(data$x) - min(data$x)) / (length(data$x))
      bigDat <-
        do.call("rbind", lapply(as.list(as.data.frame(t(data[, c("x", "y")]))), function(x)
          cbind(
            x[1] + c(0, 1, 0,-1,-1, 0, 1,-1, 1) * dist,
            x[2] + c(0, 1, 1, 1,-1,-1,-1, 0, 0) * dist
          )))
      ch <-
        concaveman::concaveman(bigDat,
                               length_threshold = window.length,
                               concavity = 1)
      poly <- as.data.frame(ch[nrow(ch):1, ])
      colnames(poly) <- c("x", "y")
      ow <-
        spatstat::owin(
          xrange = range(poly$x),
          yrange = range(poly$y),
          poly = poly
        )
      
    }
    ow
  }


#' @importFrom spatstat ppp localLcross
#' @importFrom BiocParallel bplapply
generateCurves <-
  function(data,
           Rs,
           window,
           window.length,
           BPcellType = BPcellType,
           sigma = sigma,
           ...) {
    ow <- makeWindow(data, window, window.length)
    p1 <-
      spatstat::ppp(
        x = data$x,
        y = data$y,
        window = ow,
        marks = data$cellType
      )
    
    if (!is.null(sigma)) {
      d <- spatstat::density.ppp(p1, sigma = sigma)
      d <- d / mean(d)
    }
    
    
    locIJ <-
      BiocParallel::bplapply(as.list(levels(p1$marks)), function(j) {
        locI <- lapply(as.list(levels(p1$marks)), function(i) {
          iID <- data$cellID[p1$marks == i]
          jID <- data$cellID[p1$marks == j]
          locR <- matrix(NA, length(iID), length(Rs))
          rownames(locR) <- iID
          
          if (length(jID) > 1 & length(iID) > 1) {
            if (!is.null(sigma)) {
              dFrom <- d * (sum(p1$marks == i) - 1) / spatstat::area(ow)
              dTo <-
                d * (sum(p1$marks == j) - 1) / spatstat::area(ow)
              localL <-
                spatstat::localLcross.inhom(
                  p1,
                  from = i,
                  to = j,
                  verbose = FALSE,
                  lambdaFrom = dFrom,
                  lambdaTo = dTo
                )
            } else{
              localL <-
                spatstat::localLcross(
                  p1,
                  from = i,
                  to = j,
                  verbose = FALSE
                )
            }
            ur <-
              vapply(Rs, function(x)
                which.min(abs(localL$r - x))[1], numeric(1))
            locR <-
              t(apply(as.matrix(localL)[, grep("iso", colnames(localL))],
                      2, function(x)
                        ((x - localL$theo)/localL$theo)[ur]))
            rownames(locR) <- iID
          }
          colnames(locR) <-
            paste(j, round(Rs, 2), sep = "_")
          
          locR
        })
        do.call("rbind", locI)
      }, BPPARAM = BPcellType)
    do.call("cbind", locIJ)
  }

#' @importFrom spatstat edge.Ripley nearest.valid.pixel area marks
weightCounts <- function(dt, X, maxD) {
  maxD <- as.numeric(as.character(maxD))
  
  # edge correction
  e <- spatstat::edge.Ripley(X, rep(maxD, length(X$x)))

  n <- apply(table(spatstat::marks(X)),1,function(x)x)
  A <- spatstat::area(X)
  lam <- n/A
  #V <- maxD^2/lam/pi
  #V <- (1/4 - 1/64/(lam*pi*maxD^2))maxD^2/lam/pi
  
  
  len = 300
  lambda <- (seq(1,300,length.out = len)/100)^4
  mL <- max(lambda)
  
  
  V <- NULL
  for(i in 1:len){
    V[i] <- var(sqrt(rpois(10000,lambda[i]))  ) 
  }
  
  lambda <- lambda^(1/4)
  V <- V
  
  f <- loess(V~lambda,span = 0.1)
  # plot(f,log = "xy")
  # points(lambda, fitted(f), type = "l")
  # lambda4 <- lambda^4
  # points(lambda, lambda4 - 4*lambda4^2, type = "l", col = 2)
  
  lambda <- maxD^2*lam*pi
  lambda <- lambda[marks(X)]/e
  pred <- predict(f,lambda^(1/4))
  pred[lambda < 0.001] = (lambda - 4*lambda^2)[lambda < 0.001]
  pred[lambda > mL] = 0.25
  
  
  
  # 
  # V <- lambda - c(2*lambda - 5/4*lambda^2)^2
  # 
  # 
  # V <- (maxD^2*lam*pi - (maxD^2*lam*pi)^2)
  # V[maxD^2*lam*pi>0.5] <- 0.25
  # V <- V/lam/pi
  
  V <- pred/lam[marks(X)]/pi*e
  
  
  # count and scale
  mat <- as.matrix(dt)
  #rownames(mat) <- dt$i
  # mat <- sweep(mat, 1, e, "*")
  # mat <- sqrt(mat / pi)
  # mat[is.na(mat)] <- 0
  # mat <- sweep(mat-maxD, 1, sqrt(V), "/") 
  
  mat <- sweep(mat, 1, e, "*")
  mat <- sweep(mat-pi*maxD^2, 1, sqrt(pi*maxD^2/lam[marks(X)]), "/")
  colnames(mat) <- paste(maxD, colnames(mat), sep = "_")
  mat
}

#' @importFrom spatstat ppp closepairs density.ppp marks
#' @import data.table
inhomLocalL <-
  function (data,
            Rs = c(20, 50, 100, 200),
            sigma = 10000,
            window = "convex",
            window.length = NULL,
            minLambda = 0.05) {
    ow <- makeWindow(data, window, window.length)
    X <-
      spatstat::ppp(
        x = data$x,
        y = data$y,
        window = ow,
        marks = data$cellType
      )
    
    if (is.null(Rs))
      Rs = c(20, 50, 100, 200)
    if (is.null(sigma))
      sigma = 100000
    
    Rs <- unique(c(0, sort(Rs)))
    
    den <- spatstat::density.ppp(X, sigma = sigma)
    den <- den / mean(den)
    
    p <- spatstat::closepairs(X, max(Rs), what = "ijd")
    n <- X$n
    p$j <- data$cellID[p$j]
    p$i <- data$cellID[p$i]
    
    cT <- data$cellType
    names(cT) <- data$cellID
    
    p$d <- cut(p$d, Rs, labels = Rs[-1], include.lowest = TRUE)
    
    # inhom density
    cat(length(unique(p$i)),"\n",n,"\n", length(levels(p$i)),"\n")
    np <- spatstat::nearest.valid.pixel(X$x, X$y, den)
    w <- den$v[cbind(np$row, np$col)]
    names(w) <- data$cellID
    lam <- table(marks(X))/spatstat::area(X)
    w <- w*as.numeric(lam[spatstat::marks(X)])
    D <- tapply(1/w,marks(X),sum)/spatstat::area(X)
    p$wt <- 1/w[p$j]#/D[cT[p$j]]#/w[p$i]*lam[marks(X)[p$i]]/D[marks(X)[p$i]]
    rm(np)
  

    p$j <- cT[p$j]
    p$i <- factor(p$i, levels = data$cellID)
    
    #p$wt <- 1/lam[marks(X)[p$j]]
    #p <- data.table::setDT(p)
    
    p <- as.data.frame(p)
    xt <- xtabs(wt ~ i+j+d, p)
    r <- NULL
    for(i in dimnames(xt)$d){
      r <- cbind(r, weightCounts(xt[,,i], X, i))
    }
    
    
    # p <- setkey(p, "i", "j", "d")
    # r <- p[CJ(i, j, d, unique = TRUE)][, lapply(.SD, sum), by = .(i, j, d)][is.na(wt), wt := 0]
    # message(r[1,])
    # r <- data.table::dcast(r, i + d ~ j, value.var = "wt", fill = 0)
    # r <- split(r, r$d)
    # r <- lapply(r, weightCounts, X)
    # r <- do.call("cbind", r)

    
    # 
    # r <- p[, N := sum(wt), by = .(i, j, d), drop = FALSE]
    # r <- r[,wt:=NULL]
    # r <- unique(r)
    # r <- data.table::dcast(r, i + d ~ j, value.var = "N", fill = 0)
    # r <-
    #   data.table::melt(r,
    #                    c("i", "d"),
    #                    variable.name = "cellType",
    #                    value.name = "N")
    # r <- data.table::dcast(r, i + cellType ~ d, value.var = "N", fill = 0)
    # r <-
    #   data.table::melt(r,
    #                    c("i", "cellType"),
    #                    variable.name = "d",
    #                    value.name = "N")
    # r <- r[, N := cumsum(N), by = list(i, cellType)]
    # r <-
    #   data.table::dcast(r, i + d ~ cellType, value.var = "N", fill = 0)
    # r <- split(r, r$d)
    # r <- lapply(r, weightCounts, X)
    # r <- do.call("cbind", r)
     # message("finish")
     r[data$cellID,]
    
  }

