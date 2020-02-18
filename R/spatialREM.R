#' Get statistic from pairwise L curve of a single image.
#'
#' @param x A SegmentedCellExperiment or data frame that contains at least the variables x and y, giving the location of each cell, and cellType.
#' @param from The 'from' cellType for generating the L curve.
#' @param to The 'to' cellType for generating the L curve.
#' @param dist The distance at which the statistic is obtained.
#' @param integrate Should the statistic be the integral from 0 to dist, or the value of the L curve at dist.
#'
#' @return Statistic from pairwise L curve of a single image.
#' @export
#' 
#' @import SegmentedCellExperiment spatstat
#'
#' @examples
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

#' Title
#'
#' @param cells 
#' @param from 
#' @param to 
#' @param dist 
#'
#' @return
#' 
#' @import spatstat
#'
#' @examples
getStat <- function(cells, from, to, dist) {
    pppCell <- pppGenerate(cells)
    
    L <- tryCatch({
        Lcross(pppCell, from=from, to=to, correction="best")
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

#' Title
#'
#' @param cells 
#' @param from 
#' @param to 
#' @param dist 
#'
#' @return
#' @import spatstat
#'
#' @examples
getStat2 <- function(cells, from, to, dist) {
    pppCell <- pppGenerate(cells)
    
    L <- tryCatch({
        Lcross(pppCell, from=from, to=to, correction="best")
    }, error = function(e) {
    })
    
    if (length(class(L)) == 1) {
      return(NA)
    }
    
    if (class(L) == "try-error") {
      return(NA)
    }
    
    r <- L$r[L$r <= dist]
    iso = L$iso[L$r <= dist]
    
    max(iso)
}

#' Perform mixed-effects modelling on spatial statistics.
#'
#' @param x A SegmentedCellExperiment or data frame that contains at least the variables x and y, giving the location of each cell, and cellType.
#' @param pairwise Statistic from pairwise L curve of a single image
#' @param condition Vector of conditions corresponding to each image.
#' @param from The 'from' cellType for generating the L curve.
#' @param to The 'to' cellType for generating the L curve.
#'
#' @return SegmentedCellExperiment lmerModLmerTest object.
#' @export
#' 
#' @import lme4 lmerTest mgcv
#'
#' @examples
spatialREM <- function(x, pairwise, from, to) {
  
    cells <- x
    subject <- phenotype(cells)$subject
    condition <- phenotype(cells)$condition
    
    cellSplit <- location(cells, bind = FALSE)
    count1 <- unlist(lapply(cellSplit, 
                           function(x) sum(x$cellType == from)))
    count2 <- unlist(lapply(cellSplit, 
                            function(x) sum(x$cellType == to)))
    
    spatialData <- data.frame(Pair = pairwise,
                              Subject = subject,
                              Condition = condition)
    
    count1 <- sqrt(count1)
    count2 <- sqrt(count2)
    
    filter <- !is.na(spatialData$Pair)
    
    spatialData <- spatialData[filter,]
    count1 <- count1[filter]
    count2 <- count2[filter]
    
    pairwise <- pairwise[filter]
    
    z <- (pairwise-mean(pairwise))^2
    z1 <- gam(z~ti(count1,count2))
    w <- 1/sqrt(z1$fitted.values-min(z1$fitted.values)+1)
    w <- w/sum(w)
    
    mixed.lmer <- lmer(Pair ~ Condition + (1|Subject), 
                       data = spatialData, 
                       weights = w)
    
    mixed.lmer
}

#' Performs bootstrapping to estimate p-value.
#'
#' @param mixed.lmer lmerModLmerTest object.
#' @param nsim Number of simulations.
#'
#' @return Vector of two p-values.
#' @export
#' 
#' @import lme4 lmerTest
#'
#' @examples
spatialREMBootstrap <- function(mixed.lmer, nsim=19) {
    m <- fixef(mixed.lmer)[2]
    
    bootCoef <- bootMer(mixed.lmer, 
                        fixef, 
                        nsim=nsim, 
                        re.form=NA)
    
    bootCoef <- bootCoef$t[,2]
    
    allCoefs <- c(m, bootCoef)
    
    p.val1 <- 1-sum(allCoefs>0)/length(allCoefs)
    p.val2 <- 1-sum(allCoefs<0)/length(allCoefs)
    
    p.val <- c(p.val1, p.val2)
    names(p.val) <- c("Greater", "Lesser")
    
    p.val
}

#' Performs multiple mixed-effects modelling on spatial statistics.
#'
#' @param x A SegmentedCellExperiment or data frame that contains at least the variables x and y, giving the location of each cell, and cellType.
#' @param whichCondition Vector containing the two conditions to be analysed.
#' @param dist The distance at which the statistic is obtained.
#' @param integrate Should the statistic be the integral from 0 to dist, or the value of the L curve at dist.
#' @param subject Vector of subject IDs corresponding to each image if x is a data frame.
#' @param condition Vector of conditions corresponding to each image if x is a data frame.
#' @param nsim Number of simulations to perform. If empty, the p-value from lmerTest is used.
#'
#' @return Data frame of p-values.
#' @export
#'
#' @import SegmentedCellExperiment lme4 lmerTest
#'
#' @examples
spatialREMMulti <- function(x,
                            whichCondition=NULL,
                            dist=50,
                            integrate=TRUE,
                            subject=NULL,
                            condition=NULL,
                            nsim=NULL) {

    if (is(x, "SegmentedCellExperiment")) {
        cells <- location(x, bind = FALSE)
        phenotype <- phenotype(x)
        condition <- phenotype$condition
  
        if (length(whichCondition) == 0) {
            imagesToGet <- which(condition %in% unique(condition)[c(1,2)])
            cells <- cells[imagesToGet]
            phenotype <- phenotype[imagesToGet,]
        } else if (length(whichCondition) == 2) {
            imagesToGet <- which(condition %in% whichCondition)
            cells <- cells[imagesToGet]
            phenotype <- phenotype[imagesToGet,]
            condition <- condition[imagesToget]
        } else (
            stop("Two conditions are required")
        )
    
        m <- as.character(unique(location(x)$cellType))
  
    } else if (is(x, "data.frame")  | is(x, "DataFrame")) {
  
        if (length(condition==0)) {
            stop("Specify the conditions")
        } else if (length(subject==0)) {
            stop("Specify the subjects")
        }
    
        cells <- x
        cells <- split(cells, cells$imageID)
        phenotype <- data.frame(imageID = unique(cells$imageID),
                                condition = condition,
                                subject = subject)
    
        m <- as.character(unique(x$cellType))
      } else {
          stop("x is not a data.frame, DataFrame or SegmentedCellExperiment")
      }
    
      if (length(unique(condition)) != 2) {
          stop("There must be two unique conditions")
      }
    
      m1 <- rep(m, times = length(m))
      m2 <- rep(m, each = length(m))
    
      labels <- paste(m1, m2, sep="_")
    
      MoreArgs1 <- list(x = x, dist = dist, integrate = integrate)
    
      pairwise <- mapply(getPairwise,
                         from = m1,
                         to = m2,
                         MoreArgs = MoreArgs1)
    
      pairwise[is.na(pairwise)] <- NA
    
      pairwise <- as.list(data.frame(pairwise))
      names(pairwise) <- labels
    
      MoreArgs2 <- list(cells = cells, phenotype = phenotype, nsim = nsim)
    
      mixed.lmer <- mapply(spatialREMForMulti,
                           pairwise = pairwise,
                           from = m1,
                           to = m2,
                           MoreArgs = MoreArgs2,
                           SIMPLIFY = FALSE)
    
      if (length(nsim) > 0) {
          p <- lapply(mixed.lmer, spatialREMBootstrap, nsim = nsim)
      
          p <- do.call(rbind, p)
      
          df <- data.frame(From = m1,
                           To = m2,
                           Greater = unname(p[,1]),
                           Lesser = unname(p[,2]))
      
          df
      } else {
          p <- do.call(rbind, mixed.lmer)
      
          df <- data.frame(From = m1,
                           To = m2,
                           Greater = rep(1, nrow(p)),
                           Lesser = rep(1, nrow(p)))
      
          df[p[,1]==1,3] <- unname(p[p[,1]==1,2])
          df[p[,1]==0,4] <- unname(p[p[,1]==0,2])
      
          df
      }
}

#' Title
#'
#' @param pairwise 
#' @param from 
#' @param to 
#' @param cells 
#' @param phenotype 
#' @param nsim 
#'
#' @return
#' @import lme4 lmerTest mgcv
#'
#' @examples
spatialREMForMulti <- function(pairwise, from, to, cells, phenotype, nsim) {
    spatialData <- data.frame(Pair = pairwise,
                              subject = factor(phenotype$subject),
                              Condition = factor(phenotype$condition))
  
    count1 <- unlist(lapply(cells, function(x) sum(x$cellType == from)))
    count1 <- as.numeric(count1)
    count2 <- unlist(lapply(cells, function(x) sum(x$cellType == to)))
    count2 <- as.numeric(count2)
    filter <- !is.na(spatialData$Pair)
  
    spatialData <- spatialData[filter,]
    count1 <- count1[filter]
    count2 <- count2[filter]
  
    z <- (spatialData$Pair-mean(spatialData$Pair))^2
  
    z1 <- gam(z~ti(count1,count2))
    w <- 1/sqrt(z1$fitted.values-min(z1$fitted.values)+1)
    w <- w/sum(w)
  
    mixed.lmer <- lmer(Pair ~ Condition + (1|subject),
                             data = spatialData,
                             weights = w)
  
    if (length(nsim > 0)) {
        return(mixed.lmer)
    } else {
        coefficient <- c(summary(mixed.lmer)$coefficients[2,1]>0,
                         summary(mixed.lmer)$coefficients[2,5])
        names(coefficient) <- c("Sign", "P-Value")
        return(coefficient)
    }
}


#' Plots result of spatialREMMulti.
#'
#' @param df Data frame obtained from spatialREMMulti.
#' @param fdr TRUE if FDR correction is used.
#' @param breaks Vector of 3 numbers giving breaks used in pheatmap. The first number is the minimum, the second is the maximum, the third is the number of breaks.
#' @param col Vector of colours to use in pheatmap.
#'
#' @return
#' @export
#' 
#' @import pheatmap
#'
#' @examples
spatialREMMultiPlot <- function(df,
                                fdr=FALSE,
                                breaks=c(-5, 5, 0.5),
                                col=c("blue", "white", "red")) {
    pVal <- pmin(df[,3], df[,4])
  
    marks <- unique(df[,2])
  
    if (min(pVal) == 0) {
        pVal <- pVal + 10^floor(log10(min(pVal[pVal>0])))
    }
  
    if (fdr) {
        pVal <- p.adjust(pVal, method = "fdr")
    }
  
    isGreater <- df[,3] > df[,4]
  
    pVal <- log10(pVal)
  
    pVal[isGreater] <- abs(pVal[isGreater])
  
    pVal <- matrix(pVal, nrow = length(marks))
    colnames(pVal) <- marks
    rownames(pVal) <- marks
  
    breaks <- seq(from = breaks[1], to = breaks[2], by = breaks[3])
    pal <- colorRampPalette(col)(length(breaks))
  
    heatmap <- pheatmap(pVal,
                        col = pal,
                        breaks = breaks,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE)
  
    heatmap
}

#' Title
#'
#' @param cells 
#'
#' @return
#' @import spatstat
#'
#' @examples
pppGenerate <- function(cells) {
    pppCell <- ppp(cells$x,
                   cells$y,
                   xrange = c(0,max(cells$x)),
                   yrange = c(0,max(cells$y)),
                   marks = cells$cellType)
    
    pppCell
}
