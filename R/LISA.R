
localLCurves <- function(cells,Rs = c(1,20,50,100), window = 'square', window.length = 20){


  if(is.null(cells$cellID)){
    cat("Creating cellID as it doesn't exist")
    cells$cellID <- paste('cell',seq_len(nrow(cells)),sep = '_')
  }

  if(length(cells$cellID) != length(unique(cells$cellID))) stop("The number of rows in cells does not equal the number of unique cellIDs")


  generateCurves <- BiocParallel::bplapply(as.list(unique(cells$SampleID)), function(k, cells, Rs, cellTypes, window, window.length){
    use <- cells$SampleID == k

    data <- cells[use,]
    ow <- spatstat::owin(xrange = range(data$x), yrange = range(data$y))

    if(window=='convex'){
      ch = grDevices::chull(data[,c('x','y')])
      poly <- data[,c('x','y')][rev(ch),]
      colnames(poly) = c('x','y')
      ow <- spatstat::owin(xrange = range(data$x), yrange = range(data$y), poly = poly)
    }
    if(window=='concave'){
      ch = concaveman::concaveman(as.matrix(data[,c('x','y')]),length = 20)
      poly <- as.data.frame(ch[nrow(ch):1,])
      ch = concaveman::concaveman(do.call('rbind',lapply(as.list(as.data.frame(t(poly))),function(x)cbind(rnorm(1000,x[1]),rnorm(1000,x[2])))), length = window.length)
      poly <- as.data.frame(ch[nrow(ch):1,])
      colnames(poly) = c('x','y')
      ow <- spatstat::owin(xrange = range(poly$x), yrange = range(poly$y), poly = poly)
    }


    p1 <- spatstat::ppp(x = data$x, y = data$y, window = ow, marks = as.factor(data$cellType))


    locIJ <- lapply(as.list(cellTypes), function(j){

      locI <- lapply(as.list(cellTypes),function(i){

        iID <- cells$cellID[use][p1$marks == i]
        jID <- cells$cellID[use][p1$marks == j]
        locR <- matrix(NA, length(iID), length(Rs))
        rownames(locR) = iID

        if(length(jID)>1&length(iID)>1){

          localL <- spatstat::localLcross(p1,from = i, to = j, verbose = FALSE)
          ur <- sapply(Rs,function(x)which.min(abs(localL$r-x))[1])
          locR <- t(apply(as.matrix(localL)[,grep('iso',colnames(localL))],2,function(x)(x-localL$theo)[ur]))
          rownames(locR) <- iID
        }
        colnames(locR) <- paste(j,round(Rs,2),sep = '_')

        return(locR)
      })
      return(do.call('rbind',locI))
    })
    return(do.call('cbind',locIJ))
  }, cells = cells, Rs = Rs, cellTypes = unique(cells$cellType), window = window, window.length = window.length)

  curves <- do.call('rbind',generateCurves)
  curves <- curves[cells$cellID,]
  return(curves)
}





localLCurvesInhom <- function(cells,Rs = c(1,20,50,100), window = 'square', window.length = 20){


  if(is.null(cells$cellID)){
    cat("Creating cellID as it doesn't exist")
    cells$cellID <- paste('cell',seq_len(nrow(cells)),sep = '_')
  }

  if(length(cells$cellID) != length(unique(cells$cellID))) stop("The number of rows in cells does not equal the number of unique cellIDs")


  generateCurves <- BiocParallel::bplapply(as.list(unique(cells$SampleID)), function(k, cells, Rs, cellTypes, window, window.length){
    use <- cells$SampleID == k

    data <- cells[use,]
    ow <- spatstat::owin(xrange = range(data$x), yrange = range(data$y))

    if(window=='convex'){
      ch = grDevices::chull(data[,c('x','y')])
      poly <- data[,c('x','y')][rev(ch),]
      colnames(poly) = c('x','y')
      ow <- spatstat::owin(xrange = range(data$x), yrange = range(data$y), poly = poly)
    }
    if(window=='concave'){
      ch = concaveman::concaveman(as.matrix(data[,c('x','y')]),length = 20)
      poly <- as.data.frame(ch[nrow(ch):1,])
      ch = concaveman::concaveman(do.call('rbind',lapply(as.list(as.data.frame(t(poly))),function(x)cbind(rnorm(1000,x[1]),rnorm(1000,x[2])))), length = window.length)
      poly <- as.data.frame(ch[nrow(ch):1,])
      colnames(poly) = c('x','y')
      ow <- spatstat::owin(xrange = range(poly$x), yrange = range(poly$y), poly = poly)
    }


    p1 <- spatstat::ppp(x = data$x, y = data$y, window = ow, marks = as.factor(data$cellType))


    locIJ <- lapply(as.list(cellTypes), function(j){

      locI <- lapply(as.list(cellTypes),function(i){

        iID <- cells$cellID[use][p1$marks == i]
        jID <- cells$cellID[use][p1$marks == j]
        locR <- matrix(NA, length(iID), length(Rs))
        rownames(locR) = iID

        if(length(jID)>1&length(iID)>1){

          localL <- spatstat::localLcross.inhom(p1,from = i, to = j, verbose = FALSE)
          ur <- sapply(Rs,function(x)which.min(abs(localL$r-x))[1])
          locR <- t(apply(as.matrix(localL)[,grep('iso',colnames(localL))],2,function(x)(x-localL$theo)[ur]))
          rownames(locR) <- iID
        }
        colnames(locR) <- paste(j,round(Rs,2),sep = '_')

        return(locR)
      })
      return(do.call('rbind',locI))
    })
    return(do.call('cbind',locIJ))
  }, cells = cells, Rs = Rs, cellTypes = unique(cells$cellType), window = window, window.length = window.length)

  curves <- do.call('rbind',generateCurves)
  curves <- curves[cells$cellID,]
  return(curves)
}





localLCurves2 <- function(cells,Rs = c(1,20,50,100), window = 'square', window.length = 20){


  if(is.null(cells$cellID)){
    cat("Creating cellID as it doesn't exist")
    cells$cellID <- paste('cell',seq_len(nrow(cells)),sep = '_')
  }

  if(length(cells$cellID) != length(unique(cells$cellID))) stop("The number of rows in cells does not equal the number of unique cellIDs")


  generateCurves <- lapply(as.list(unique(cells$SampleID)), function(k, cells, Rs, cellTypes, window, window.length){
    use <- cells$SampleID == k

    data <- cells[use,]
    ow <- spatstat::owin(xrange = range(data$x), yrange = range(data$y))

    if(window=='convex'){
      ch = grDevices::chull(data[,c('x','y')])
      poly <- data[,c('x','y')][rev(ch),]
      colnames(poly) = c('x','y')
      ow <- spatstat::owin(xrange = range(data$x), yrange = range(data$y), poly = poly)
    }
    if(window=='concave'){
      ch = concaveman::concaveman(as.matrix(data[,c('x','y')]),length = 20)
      poly <- as.data.frame(ch[nrow(ch):1,])
      ch = concaveman::concaveman(do.call('rbind',lapply(as.list(as.data.frame(t(poly))),function(x)cbind(rnorm(1000,x[1]),rnorm(1000,x[2])))), length = window.length)
      poly <- as.data.frame(ch[nrow(ch):1,])
      colnames(poly) = c('x','y')
      ow <- spatstat::owin(xrange = range(poly$x), yrange = range(poly$y), poly = poly)
    }


    p1 <- spatstat::ppp(x = data$x, y = data$y, window = ow, marks = as.factor(data$cellType))


    locIJ <- BiocParallel::bplapply(as.list(cellTypes), function(j,cells, Rs){

      locI <- lapply(as.list(cellTypes),function(i){

        iID <- cells$cellID[use][p1$marks == i]
        jID <- cells$cellID[use][p1$marks == j]
        locR <- matrix(NA, length(iID), length(Rs))
        rownames(locR) = iID

        if(length(jID)>1&length(iID)>1){

          localL <- spatstat::localLcross(p1,from = i, to = j, verbose = FALSE)
          ur <- sapply(Rs,function(x)which.min(abs(localL$r-x))[1])
          locR <- t(apply(as.matrix(localL)[,grep('iso',colnames(localL))],2,function(x)(x-localL$theo)[ur]))
          rownames(locR) <- iID
        }
        colnames(locR) <- paste(j,round(Rs,2),sep = '_')

        return(locR)
      })
      return(do.call('rbind',locI))
    }, cells = cells, Rs = Rs)
    return(do.call('cbind',locIJ))
  }, cells = cells, Rs = Rs, cellTypes = unique(cells$cellType), window = window, window.length = window.length)

  curves <- do.call('rbind',generateCurves)
  curves <- curves[cells$cellID,]
  return(curves)
}





localLCurves2Inhom <- function(cells,Rs = c(1,20,50,100), window = 'square', window.length = 20){


  if(is.null(cells$cellID)){
    cat("Creating cellID as it doesn't exist")
    cells$cellID <- paste('cell',seq_len(nrow(cells)),sep = '_')
  }

  if(length(cells$cellID) != length(unique(cells$cellID))) stop("The number of rows in cells does not equal the number of unique cellIDs")


  generateCurves <- lapply(as.list(unique(cells$SampleID)), function(k, cells, Rs, cellTypes, window, window.length){
    use <- cells$SampleID == k

    data <- cells[use,]
    ow <- spatstat::owin(xrange = range(data$x), yrange = range(data$y))

    if(window=='convex'){
      ch = grDevices::chull(data[,c('x','y')])
      poly <- data[,c('x','y')][rev(ch),]
      colnames(poly) = c('x','y')
      ow <- spatstat::owin(xrange = range(data$x), yrange = range(data$y), poly = poly)
    }
    if(window=='concave'){
      ch = concaveman::concaveman(as.matrix(data[,c('x','y')]),length = 20)
      poly <- as.data.frame(ch[nrow(ch):1,])
      ch = concaveman::concaveman(do.call('rbind',lapply(as.list(as.data.frame(t(poly))),function(x)cbind(rnorm(1000,x[1]),rnorm(1000,x[2])))), length = window.length)
      poly <- as.data.frame(ch[nrow(ch):1,])
      colnames(poly) = c('x','y')
      ow <- spatstat::owin(xrange = range(poly$x), yrange = range(poly$y), poly = poly)
    }


    p1 <- spatstat::ppp(x = data$x, y = data$y, window = ow, marks = as.factor(data$cellType))


    locIJ <- BiocParallel::bplapply(as.list(cellTypes), function(j,cells, Rs){

      locI <- lapply(as.list(cellTypes),function(i){

        iID <- cells$cellID[use][p1$marks == i]
        jID <- cells$cellID[use][p1$marks == j]
        locR <- matrix(NA, length(iID), length(Rs))
        rownames(locR) = iID

        if(length(jID)>1&length(iID)>1){

          localL <- spatstat::localLcross.inhom(p1,from = i, to = j, verbose = FALSE)
          ur <- sapply(Rs,function(x)which.min(abs(localL$r-x))[1])
          locR <- t(apply(as.matrix(localL)[,grep('iso',colnames(localL))],2,function(x)(x-localL$theo)[ur]))
          rownames(locR) <- iID
        }
        colnames(locR) <- paste(j,round(Rs,2),sep = '_')

        return(locR)
      })
      return(do.call('rbind',locI))
    }, cells = cells, Rs = Rs)
    return(do.call('cbind',locIJ))
  }, cells = cells, Rs = Rs, cellTypes = unique(cells$cellType), window = window, window.length = window.length)

  curves <- do.call('rbind',generateCurves)
  curves <- curves[cells$cellID,]
  return(curves)
}
