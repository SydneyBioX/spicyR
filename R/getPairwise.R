getStat <- function(pppCell, from, to, dist) {
  L <- Lcross(cell, from=from, to=to, correction="best")
  
  r <- L$r[L$r <= d]
  iso = L$iso[L$r <= d]
  factor = min(r[r > 0])
  
  sum(iso)*factor
  
  r = LVal$r[LVal$r<=d]
  iso = LVal$iso[LVal$r<=d]
  statmat[i,j] = max(iso)
  
}

getStat2 <- function(pppCell, from, to, dist) {
    L <- Lcross(cell, from=from, to=to, correction="best")
    
    r <- L$r[L$r <= d]
    iso = L$iso[L$r <= d]
    
    max(iso)
}

getPairwise <- function(cells, from, to, dist=50, integrate=TRUE){
    cellSplit <- split(cells,cells$patient)
    pppCell <- lapply(cellSplit,pppGenerate)
    
    if (integrate) lapply(pppCell, getStat, from, to, dist)
    else lapply(pppCell, getStat2, from, to, dist)
}

