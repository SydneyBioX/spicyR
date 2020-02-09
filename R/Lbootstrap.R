Lbootstrap <- function(x, 
                       from, 
                       to, 
                       imageID=1,
                       rinterval=c(0,100),
                       nsim=19,
                       alternative=c("two.sided", "less", "greater")) {
  
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
    
    E <- spatstat::envelope(pppCell,
                            spatstat::Lcross,
                            nsim = nsim,
                            from = from,
                            to = to,
                            simulate = expression(spatstat::rlabel(pppCell)),
                            correction = "best",
                            savepatterns = TRUE)
    
    test <- spatstat::dclf.test(E,
                                rinterval = rinterval,
                                alternative = alternative)
    
    test
}

LBootstrapMulti <- function(x, 
                            pVal = 0.05,
                            rinterval=c(0,100), 
                            nsim=19) {
  
  if (is(x, "SegmentedCellExperiment")) {
    cells <- location(x, bind = FALSE)
  } else if (is(x, "data.frame")  | is(x, "DataFrame")) {
    cells <- x
    cells <- split(cells, cells$imageID)
  } else {
    stop("x is not a data.frame, DataFrame or SegmentedCellExperiment")
  }
  
  m <- as.character(unique(location(x)$cellType))
  m1 <- rep(m, times = length(m))
  m2 <- rep(m, each = length(m))
  
  labels <- paste(m1, m2, sep="_")
  
  rankList <- lapply(cells, 
                     LGetRankMulti,
                     allLabels = labels,
                     seed = seed, 
                     rinterval = rinterval, 
                     nsim = nsim)
  
  rankMat <- do.call(rbind, rankList)
  
  pValLocalisation <- rankMat/nsim*2
  pValAntiLocalisation <- (nsim+1-rankMat)/nsim*2
  
  percentSignificantLocalisation <- 
    colSums(pValLocalisation < pVal, na.rm = TRUE) /
    nrow(pValLocalisation)*100
  percentSignificantLocalisation <- 
    matrix(percentSignificantLocalisation,
           nrow = length(m),
           ncol = length(m),
           dimnames = list(m,m))
  
  percentSignificantAntiLocalisation <- 
    colSums(pValAntiLocalisation < pVal, na.rm = TRUE) / 
    nrow(pValAntiLocalisation)*100
  percentSignificantAntiLocalisation <- 
    matrix(percentSignificantAntiLocalisation,
           nrow = length(m),
           ncol = length(m),
           dimnames = list(m,m))
  
  percentSignificant <- list(percentSignificantLocalisation,
                             percentSignificantAntiLocalisation)
  names(percentSignificant) <- c("Localisation", "Anti-Localisation")
}

LGetRankMulti <- function(cells, 
                          allLabels,
                          rinterval=c(0,100), 
                          nsim=19) {
  
  pppCell <- pppGenerate(cells)
  
  m <- as.character(unique(pppCell$marks))
  
  m1 <- rep(m, times = length(m))
  m2 <- rep(m, each = length(m))
  
  labels <- paste(m1, m2, sep="_")
  
  MoreArgs = list(pppCell=pppCell, rinterval = rinterval, nsim = nsim)
  
  ranks <- mapply(LGetRank,
                  from = m1,
                  to = m2,
                  MoreArgs = MoreArgs)
  names(ranks) <- labels
  
  
  allRanks <- rep(NA, length(allLabels))
  names(allRanks) <- allLabels
  allRanks[names(ranks)] <- ranks
  
  allRanks
}

LGetRank <- function(pppCell,
                     from=from, 
                     to=to, 
                     rinterval=c(0,50), 
                     nsim=19) {
  
  E <- spatstat::envelope(pppCell,
                          spatstat::Lcross,
                          nsim = nsim,
                          from = from,
                          to = to,
                          simulate = expression(spatstat::rlabel(pppCell)),
                          correction = "best",
                          savepatterns = TRUE)
  
  test <- spatstat::dclf.test(E, 
                              rinterval = rinterval, 
                              alternative = "greater")
  
  test$statistic$rank
}