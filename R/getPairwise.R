getStat <- function(cells, from, to, d) {
    pppCell <- pppGenerate(cells)
    
    L <- spatstat::Lcross(pppCell, from=from, to=to, correction="best")
    
    r <- L$r[L$r <= d]
    iso = L$iso[L$r <= d]
    factor = min(r[r > 0])
    
    sum(iso)*factor
    
    r = L$r[L$r<=d]
    iso = L$iso[L$r<=d]
    
    max(iso)
}

getStat2 <- function(cells, from, to, d) {
    pppCell <- pppGenerate(cells)
    
    L <- spatstat::Lcross(pppCell, from=from, to=to, correction="best")
    
    r <- L$r[L$r <= d]
    iso = L$iso[L$r <= d]
    
    max(iso)
}

getPairwise <- function(x, from, to, dist=50, integrate=TRUE) {
    if (is(x, "SegmentedCellExperiment")) {
      cells <- location(cellExp, bind = FALSE)
    } else if (is(x, "data.frame")) {
      cells <- split(x, x$imageID)
    } else {
      stop("x is not a data.frame or SegmentedCellExperiment")
    }
  
    pairwiseVals <- lapply(cells, 
                           getStat, 
                           from = "CD8-PD1-PDL1-",
                           to = "CD8-PD1-PDL1-",
                           d = 50)
    
    if (integrate) {
        pairwiseVals <- lapply(cells, 
                               getStat, 
                               from = "CD8-PD1-PDL1-",
                               to = "CD8-PD1-PDL1-",
                               d = 50)
    } else {
        pairwiseVals <- lapply(cells, 
                               getStat, 
                               from = "CD8-PD1-PDL1-",
                               to = "CD8-PD1-PDL1-",
                               d = 50)
    }
    
    unlist(pairwiseVals)
}
