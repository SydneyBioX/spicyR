getPairwise <- function(x, from, to, dist=50, integrate=TRUE) {
  if (is(x, "SegmentedCellExperiment")) {
    cells <- location(x, bind = FALSE)
  } else if (is(x, "data.frame")) {
    cells <- split(x, x$imageID)
  } else {
    stop("x is not a data.frame or SegmentedCellExperiment")
  }
  
  if (integrate) {
    pairwiseVals <- lapply(cells, 
                           getStat, 
                           from = from,
                           to = to,
                           dist = dist)
  } else {
    pairwiseVals <- lapply(cells, 
                           getStat2, 
                           from = from,
                           to = to,
                           dist = dist)
  }
  
  unlist(pairwiseVals)
}

getStat <- function(cells, from, to, dist) {
  pppCell <- pppGenerate(cells)
  
  L <- tryCatch({
    spatstat::Lcross(pppCell, from=from, to=to, correction="best")
  }, error = function(e) {
  })
  
  if (length(class(L)) == 1) {
    return(NA)
  }
  
  r <- L$r[L$r <= dist]
  iso = L$iso[L$r <= dist]
  factor = min(r[r > 0])
  
  sum(iso)*factor
  
  r = L$r[L$r <= dist]
  iso = L$iso[L$r <= dist]
  
  max(iso)
}

getStat2 <- function(cells, from, to, dist) {
  pppCell <- pppGenerate(cells)
  
  L <- try(spatstat::Lcross(pppCell, from=from, to=to, correction="best"))
  
  if (class(L) == "try-error") {
    return(NA)
  }
  
  r <- L$r[L$r <= dist]
  iso = L$iso[L$r <= dist]
  
  max(iso)
}

spatialREM <- function(x, donor, condition, count1, count2) {
  
  spatialData <- data.frame(Pair = x,
                            Donor = donor,
                            Condition = condition)
  
  count1 <- sqrt(count1)
  count2 <- sqrt(count2)
  
  pairwiseVals <- x
  
  z <- (pairwiseVals-mean(pairwiseVals))^2
  z1 <- mgcv::gam(z~ti(x,y))
  w <- 1/sqrt(z1$fitted.values-min(z1$fitted.values)+1)
  w <- w/sum(w)
  
  mixed.lmer <- lme4::lmer(Pair ~ Group + (1|Donor), 
                           data = spatialData, 
                           weights = w)
  
  mixed.lmer
}

spatialREMBootstrap <- function(mixed.lmer, nsim=19) {
  m <- lme4::fixef(mixed.lmer)[2]
  
  bootCoef <- lme4::bootMer(mixed.lmer, 
                            lme4::fixef, 
                            nsim=nsim, 
                            re.form = NA)
  
  bootCoef <- bootCoef$t[,2]
  
  allCoefs <- c(m, bootCoef)
  
  p.val1 <- 1-sum(allCoefs>0)/length(allCoefs)
  p.val2 <- 1-sum(allCoefs<0)/length(allCoefs)
  
  p.val <- c(p.val1, p.val2)
  names(p.val) <- c("Greater", "Lesser")
  
  p.val
}

