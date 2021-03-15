#' Performs spatial tests on spatial cytometry data.
#'
#' @param cells A SegmentedCells or data frame that contains at least the 
#' variables x and y, giving the location coordinates of each cell, and cellType.
#' @param condition Vector of conditions to be tested corresponding to each image
#'  if cells is a data frame.
#' @param subject Vector of subject IDs corresponding to each image if cells is
#'  a data frame.
#' @param covariates Vector of covariate names that should be included in the 
#' mixed effects model as fixed effects.
#' @param from vector of cell types which you would like to compare to the to vector
#' @param to vector of cell types which you would like to compare to the from vector
#' @param dist The distance at which the statistic is obtained.
#' @param integrate Should the statistic be the integral from 0 to dist, or the 
#' value of the L curve at dist.
#' @param nsim Number of simulations to perform. If empty, the p-value from lmerTest is used.
#' @param verbose logical indicating whether to output messages.
#' @param weights logical indicating whether to include weights based on cell counts.
#' @param window Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length A tuning parameter for controlling the level of concavity 
#' when estimating concave windows.
#' @param BPPARAM A BiocParallelParam object.
#' @param sigma A numeric variable used for scaling when fitting inhomogeneous L-curves.
#' @param minLambda Minimum value for density for scaling when fitting inhomogeneous L-curves.
#' @param Rs A vector of the radii that the measures of association should be calculated.
#' @param fast A logical describing whether to use a fast approximation of the 
#' inhomogeneous L-curves.
#' @param ... Other options to pass to bootstrap.
#'
#' @return Data frame of p-values.
#' @export
#'
#' @examples
#' data("diabetesData")
#'
#' # Test with random effect for patient on only one pairwise combination of cell types.
#' spicy(diabetesData, condition = "stage", subject = "case", 
#'       from = "Tc", to = "Th")
#' 
#' # Test all pairwise combination of cell types without random effect of patient.
#' #spicyTest <- spicy(diabetesData, condition = "stage", subject = "case")
#'
#' # Test all pairwise combination of cell types with random effect of patient.
#' #spicy(diabetesData, condition = "condition", subject = "subject")
#'
#' # Test all pairwise combination of cell types with random effect of patient using 
#' # a bootstrap to calculate significance.
#' #spicy(diabetesData, condition = "stage", subject = "case", nsim = 10000)
#' 
#' @aliases
#' spicy
#' spicy,spicy-method
#' @importFrom mgcv gam ti
#' @importFrom BiocParallel SerialParam
spicy <- function(cells,
                  condition = NULL,
                  subject = NULL,
                  covariates = NULL,
                  from = NULL,
                  to = NULL,
                  dist = NULL,
                  integrate = TRUE,
                  nsim = NULL,
                  verbose = TRUE,
                  weights = TRUE,
                  window = "convex",
                  window.length = NULL,
                  BPPARAM = BiocParallel::SerialParam(),
                  sigma = NULL,
                  Rs = NULL,
                  minLambda = 0.05,
                  fast = TRUE,
                  ...) {
    if (!is(cells, "SegmentedCells")) {
        stop('cells needs to be a SegmentedCells object')
    }
    
    if (is.null(from))
        from <- as.character(unique(cellType(cells)))
    if (is.null(to))
        to <- as.character(unique(cellType(cells)))
    
    from <- as.character(unique(from))
    to <- as.character(unique(to))
    
    if (any((!to %in% cellType(cells)) |
            (!from %in% cellType(cells))))
        stop("to and from need to be cell type in your SegmentedCells")
    
    nCells <- table(imageID(cells), cellType(cells))
    
    
    m1 <- rep(from, times = length(to))
    m2 <- rep(to, each = length(from))
    labels <- paste(m1, m2, sep = "_")
    
    ## Find pairwise associations
    
    if(fast){
        
        pairwiseAssoc <- getPairwise(cells,
                                Rs = Rs,
                                sigma = sigma,
                                window = window,
                                window.length = window.length,
                                minLambda = minLambda,
                                from = from,
                                to = to,
                                fast = fast)
        
    }else{
    
    

    
    MoreArgs1 <- list(cells = cells, dist = dist, window = window, window.length = window.length, fast = FALSE)
    
    if (verbose)
        message("Calculating pairwise spatial associations")
    
    pairwiseAssoc <- mapply(getPairwise,
                            from = m1,
                            to = m2,
                            MoreArgs = MoreArgs1)
    colnames(pairwiseAssoc) <- labels
    
    }
    
    count1 <- as.vector(nCells[, m1])
    count2 <- as.vector(nCells[, m2])
    
    resSq <-
        as.vector(apply(pairwiseAssoc, 2, function(x)
            (x - mean(x, na.rm = TRUE))^2))
    
    toWeight <- !is.na(as.vector(pairwiseAssoc))
    resSqToWeight <- resSq[toWeight]
    count1ToWeight <- count1[toWeight]
    count2ToWeight <- count2[toWeight]
    
    if (weights) {
        weightFunction <- mgcv::gam(resSqToWeight ~ ti(count1ToWeight, count2ToWeight))
    } else {
        weightFunction <- NULL
    }
    
    pairwiseAssoc <- as.list(data.frame(pairwiseAssoc))
    names(pairwiseAssoc) <- labels
    
    ## Linear model
    if (is.null(subject) & !is.null(condition)) {
        if (verbose)
            message("Testing for spatial differences across conditions")
        
        MoreArgs2 <-
            list(
                cells = cells,
                condition = condition,
                covariates = covariates,
                weightFunction = weightFunction,
                cellCounts <- table(imageID(cells), cellType(cells)),
                pheno <- as.data.frame(imagePheno(cells))
            )
        
        
        linearModels <- mapply(
            spatialLM,
            spatAssoc = pairwiseAssoc,
            from = m1,
            to = m2,
            MoreArgs = MoreArgs2,
            SIMPLIFY = FALSE
        )
        
        df <- cleanLM(linearModels, nsim, BPPARAM = BPPARAM)
    }
    
    ## Mixed effects model
    if ((!is.null(subject)) & !is.null(condition)) {
        if (verbose)
            message(
                "Testing for spatial differences across conditions accounting for multiple images per subject"
            )
        
        MoreArgs2 <-
            list(
                cells = cells,
                subject = subject,
                condition = condition,
                covariates = covariates,
                weightFunction = weightFunction,
                cellCounts <- table(imageID(cells), cellType(cells)),
                pheno <- as.data.frame(imagePheno(cells))
            )
        
        mixed.lmer <- mapply(
            spatialMEM,
            spatAssoc = pairwiseAssoc,
            from = m1,
            to = m2,
            MoreArgs = MoreArgs2,
            SIMPLIFY = FALSE
        )

        mixed.lmer <- lapply(mixed.lmer, function(x){
            if(x@devcomp$cmp["REML"]==-Inf)return(NA)
            x
        })

        df <- cleanMEM(mixed.lmer, nsim, BPPARAM = BPPARAM)

    }
    
    
    
    df$pairwiseAssoc <- pairwiseAssoc
    df$comparisons <- data.frame(from = m1,
                                 to = m2,
                                 labels = labels)
    
    df <- new('SpicyResults', df)
    df
}

#' @importFrom dplyr bind_rows
cleanLM <- function(linearModels, nsim,  BPPARAM) {
    if (length(nsim) > 0) {
        boot <- bplapply(linearModels, spatialLMBootstrap, nsim = nsim, BPPARAM = BPPARAM)
        #p <- do.call(rbind, p)
        tBoot <- lapply(boot, function(coef) {
            if(length(coef)>1){
            coef <- as.data.frame(t(coef))
            coef <-
                split(coef, c("coefficient", "se", "statistic", "p.value"))
            }else{
                coef <- list(coefficient = NA, p.value = NA, se = NA, statistics = NA)
            }
            coef
            
        })
        
        df <- apply(do.call(rbind, tBoot), 2, function(x)
            do.call(rbind, x))
        df <- lapply(df, function(x) {
            rownames(x) <- names(linearModels)
            x
        })
    } else {
        tLm <- lapply(linearModels, function(LM) {
            if(is(LM,"lm")){
            coef <- as.data.frame(t(summary(LM)$coef))
            coef <-
                split(coef, c("coefficient", "se", "statistic", "p.value"))
            }else{
                coef <- list(coefficient = NA, p.value = NA, se = NA, statistics = NA)
            }
            coef
        })
        
        df <- apply(do.call(rbind, tLm), 2, function(x){
            x[is.na(x)] <- list(data.frame(`(Intercept)`=NA, check.names=FALSE))
            dplyr::bind_rows(x)})
        df <- lapply(df, function(x) {
            x <- as.data.frame(x)
            rownames(x) <- names(linearModels)
            x
        })
    }
    df
}

cleanMEM <- function(mixed.lmer, nsim, BPPARAM) {
    if (length(nsim) > 0) {

        boot <- BiocParallel::bplapply(mixed.lmer, spatialMEMBootstrap, nsim = nsim, BPPARAM = BPPARAM)
        #p <- do.call(rbind, p)
        tBoot <- lapply(boot, function(coef) {
            if(length(coef)>1){
                coef <- as.data.frame(t(coef))
                coef <-
                    split(coef, c("coefficient", "se", "statistic", "p.value"))
            }else{
                coef <- list(coefficient = NA, p.value = NA, se = NA, statistics = NA)
            }
            coef
        })
        
        
        df <- apply(do.call(rbind, tBoot), 2, function(x)
            do.call(rbind, x))
        df <- lapply(df, function(x) {
            rownames(x) <- names(mixed.lmer)
            x
        })
        
    } else {
        tLmer <- lapply(mixed.lmer, function(lmer) {
            if(is(lmer,"lmerMod")){
            coef <- as.data.frame(t(summary(lmer)$coef))
            coef <-
                split(coef,
                      c("coefficient", "se", "df", "statistic", "p.value"))
        }else{
            coef <- list(coefficient = NA, p.value = NA, se = NA, statistics = NA)
        }
            coef
            })
        
        df <- apply(do.call(rbind, tLmer), 2, function(x)
            do.call(rbind, x))
        df <- lapply(df, function(x) {
            rownames(x) <- names(mixed.lmer)
            x
        })
        
    }
    df
}

#' Get statistic from pairwise L curve of a single image.
#'
#' @param cells A SegmentedCells or data frame that contains at least the 
#' variables x and y, giving the location coordinates of each cell, and cellType.
#' @param from The 'from' cellType for generating the L curve.
#' @param to The 'to' cellType for generating the L curve.
#' @param dist The distance at which the statistic is obtained.
#' @param window Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length A tuning parameter for controlling the level of concavity 
#' when estimating concave windows.
#' @param Rs A vector of the radii that the measures of association should be calculated.
#' @param sigma A numeric variable used for scaling when fitting inhomogeneous L-curves.
#' @param minLambda Minimum value for density for scaling when fitting inhomogeneous L-curves.
#' @param fast A logical describing whether to use a fast approximation of the 
#' inhomogeneous L-curves.
#' @param BPPARAM A BiocParallelParam object.
#'
#' @return Statistic from pairwise L curve of a single image.
#'
#'
#' @examples
#' data("diabetesData")
#' pairAssoc <- getPairwise(diabetesData)
#' @export
#' @importFrom BiocParallel bplapply
getPairwise <- function(cells, from = unique(cellType(cells)), to = unique(cellType(cells)), dist = NULL, window = "convex", window.length, Rs = NULL, sigma = NULL, minLambda = 0.05, fast = TRUE, BPPARAM=BiocParallel::SerialParam()) {
    cells2 <- cellSummary(cells, bind = FALSE)
    
    if(fast){
        pairwiseVals <- BiocParallel::bplapply(cells2,
                               inhomLPair,
                               Rs = c(20, 50),
                               sigma = sigma,
                               window = window,
                               window.length = window.length,
                               minLambda = minLambda,
                               from = from,
                               to = to, BPPARAM=BPPARAM)
        return(do.call("rbind",pairwiseVals ))
    }else{

    pairwiseVals <- lapply(cells2,
                           getStat,
                           from = from,
                           to = to,
                           dist = dist, window, window.length)
    
    return(unlist(pairwiseVals))
    }
    

}


#' @importFrom spatstat.core Lcross
getStat <- function(cells, from, to, dist, window, window.length) {
    pppCell <- pppGenerate(cells, window, window.length)
    
    L <- tryCatch({
        spatstat.core::Lcross(pppCell,
               from = from,
               to = to,
               correction = "best")
    }, error = function(e) {
        
    })
    
    if (!is(L,"fv")) {
        return(NA)
    }
    
    if (is.null(dist))
        dist <- max(L$r)
    
    theo = L$theo[L$r <= dist]
    iso = L$iso[L$r <= dist]
    mean(iso - theo)
}

# Performs bootstrapping to estimate p-value.

#' @importFrom lme4 fixef
#' @importFrom lmerTest lmer
#' @importFrom stats formula weights
spatialMEMBootstrap <- function(mixed.lmer, nsim = 19) {

    functionToReplicate <- function(x) {
        toGet <- sample(nrow(x@frame), replace = TRUE)
        
        spatialData <- x@frame[toGet,]
        spatialData$weights <- spatialData$`(weights)`
        
        mixed.lmer1 <- tryCatch({lmerTest::lmer(formula(x),
                                     data = spatialData,
                                     weights = weights)},
                               error = function(e) {
                                   
                               })
        if (!is(mixed.lmer1,"lmerMod")) {
            return(c(NA,NA))
        }
        
        summary(mixed.lmer1)$coef[, "t value"]
    }
    if(!is(mixed.lmer, "lmerMod"))return(NA)
    stats <- replicate(nsim, functionToReplicate(x = mixed.lmer))
    stats <- t(stats)
    fe <- lme4::fixef(mixed.lmer)
    s1 <- sweep(stats,2,sign(fe),"*")
    pval <- colMeans(s1 < 0, na.rm = TRUE)*2 #+ colMeans(sweep(s1,2,2*abs(summary(mixed.lmer)$coef[, "t value"]),">"))
    
    df <-
        data.frame(
            coefficient = fe,
            se = apply(stats, 2, stats::sd, na.rm = TRUE),
            statistic = fe/apply(stats, 2, stats::sd, na.rm = TRUE),
            p.value = pval
        )
    df
}

#' @importFrom stats p.adjust
.show_SpicyResults <- function(df) {
    pval <- as.data.frame(df$p.value)
    cond <- colnames(pval)[grep('condition', colnames(pval))]
    message(df$test)
    message("Number of cell type pairs: ", nrow(pval), "\n")
    message("Number of differentially localised cell type pairs: \n")
    if (nrow(pval) == 1)
        print(sum(pval[cond] < 0.05, na.rm = TRUE))
    if (nrow(pval) > 1)
        print(colSums(apply(pval[cond], 2, p.adjust, 'fdr') < 0.05, na.rm = TRUE))
    
}
setMethod("show", signature(object = "SpicyResults"), function(object) {
    .show_SpicyResults(object)
})

#' @importFrom lmerTest lmer
#' @importFrom stats predict weights
#' @importFrom methods is
spatialMEM <-
    function(spatAssoc,
             from,
             to,
             cells,
             subject,
             condition,
             covariates,
             weightFunction,
             cellCounts,
             pheno) {
        
        spatAssoc[is.na(spatAssoc)] <- 0

        count1 <- cellCounts[, from]
        count2 <- cellCounts[, to]
        # filter <- !is.na(spatAssoc)
        # 
        # if (sum(filter) < 3)
        #     return(NA)
        
        spatialData <-
            data.frame(spatAssoc,
                       condition = pheno[, condition],
                       subject = pheno[, subject],
                       pheno[covariates])
        
        if (is.null(weightFunction)) {
            w <- rep(1, length(count1))
        } else{
            z1 <- predict(weightFunction, data.frame(count1ToWeight = as.numeric(count1), 
                                                     count2ToWeight = as.numeric(count2)))
            w <- 1 / sqrt(z1 - min(z1) + 1)
            w <- w / sum(w)
        }
        
        formula <- 'spatAssoc ~ condition + (1|subject)'
        
        if (!is.null(covariates))
            formula <-
            paste('spatAssoc ~ condition + (1|subject)',
                  paste(covariates, collapse = '+'),
                  sep = "+")
        spatialData$weights = w
        
        mixed.lmer <- tryCatch({lmerTest::lmer(formula(formula),
                           data = spatialData,
                           weights = weights)},
                           error = function(e) {
                               
                           })
        if (!is(mixed.lmer,"lmerMod")) {
            return(NA)
        }
        
        mixed.lmer
    }

#' @importFrom stats predict lm
spatialLM <-
    function(spatAssoc,
             from,
             to,
             cells,
             condition,
             covariates,
             weightFunction,
             cellCounts,
             pheno) {
        
        count1 <- cellCounts[, from]
        count2 <- cellCounts[, to]
        
        #print(pheno)
        spatialData <-
            data.frame(spatAssoc, condition = pheno[, condition], pheno[,covariates])
        # spatialData <- spatialData[filter, ]
        # count1 <- count1[filter]
        # count2 <- count2[filter]
        #print(spatialData)
        
        if (is.null(weightFunction)) {
            w <- rep(1, length(count1))
        } else {
            z1 <- predict(weightFunction, data.frame(count1ToWeight = as.numeric(count1), 
                                                     count2ToWeight = as.numeric(count2)))
            w <- 1 / sqrt(z1 - min(z1) + 1)
            w <- w / sum(w)
        }
        
        formula <- 'spatAssoc ~ condition'
        
        if (!is.null(covariates))
            formula <-
            paste('spatAssoc ~ condition',
                  paste(covariates, collapse = '+'),
                  sep = "+")
        lm1 <- tryCatch({lm(formula(formula),
                  data = spatialData,
                  weights = w)},
                   error = function(e) {
                      
                  })
                  
                  if (!is(lm1,"lm")) {
                      return(NA)
                  }
                  
        
        lm1
    }

#' @importFrom stats sd coef
spatialLMBootstrap <- function(linearModels, nsim=19) {
    functionToReplicate <- function(x) {
        
        toGet <- sample(nrow(x$model), replace = TRUE)
        
        spatAssocBoot <- x$model$spatAssoc[toGet]
        conditionBoot <- x$model$condition[toGet]
        weightsBoot <- x$weights[toGet]
        
        spatialDataBoot <- data.frame(spatAssoc = spatAssocBoot,
                                      condition = conditionBoot)
        
        lm1 <- tryCatch({lm(spatAssoc ~ condition,
                  data = spatialDataBoot,
                  weights = weightsBoot)},
                  error = function(e) {
                      
                  })
        if (!is(lm1,"lm")) {
            return(NA)
        }
        
        lm1$coefficients[2]
    }
    if(!is(linearModels, "lm"))return(NA)
    stats <- replicate(nsim, functionToReplicate(x = linearModels))
    
    df <- coef(summary(linearModels))
    
    
    pval <- pmin(mean(stats < 0,na.rm = TRUE), mean(stats > 0, na.rm = TRUE), na.rm = TRUE) * 2
    df[2,4] <- pval
    df
}

#' Plots result of signifPlot.
#'
#' @param results Data frame obtained from spicy.
#' @param fdr TRUE if FDR correction is used.
#' @param breaks Vector of 3 numbers giving breaks used in pheatmap. The first 
#' number is the minimum, the second is the maximum, the third is the number of breaks.
#' @param colors Vector of colours to use in pheatmap.
#' @param marksToPlot Vector of marks to include in pheatmap.
#'
#' @return a pheatmap object
#'
#' @examples
#' data(spicyTest)
#' signifPlot(spicyTest, breaks=c(-3, 3, 0.5))
#' 
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom stats p.adjust
signifPlot <- function(results,
                       fdr = FALSE,
                       breaks = c(-5, 5, 0.5),
                       colors = c("blue", "white", "red"),
                       marksToPlot = NULL) {
    pVal <- results$p.value[,2]
    marks <- unique(results$comparisons$from)
    
    if (is.null(marksToPlot)) marksToPlot <- marks
    
    if (min(pVal,na.rm=TRUE) == 0) {
        pVal[pVal == 0] <-
            pVal[pVal == 0] + 10 ^ floor(log10(min(pVal[pVal > 0],na.rm = TRUE)))
    }
    
    if (fdr) {
        pVal <- p.adjust(pVal, method = "fdr")
    }
    
    isGreater <- results$coefficient[,2] > 0
    
    pVal <- log10(pVal)
    
    pVal[isGreater] <- abs(pVal[isGreater])
    
    pVal <- matrix(pVal, nrow = length(marks))
    colnames(pVal) <- marks
    rownames(pVal) <- marks
    
    
    breaks <- seq(from = breaks[1], to = breaks[2], by = breaks[3])
    pal <- colorRampPalette(colors)(length(breaks))
    
    heatmap <- pheatmap(
        pVal[marksToPlot, marksToPlot],
        color = pal,
        breaks = breaks,
        cluster_rows = FALSE,
        cluster_cols = FALSE
    )
    
    heatmap
}

###########################
#
#  Generate distances
#
###########################


#' @importFrom spatstat.geom owin convexhull ppp
#' @importFrom concaveman concaveman
makeWindow <-
    function(data,
             window = "square",
             window.length = NULL) {
        data = data.frame(data)
        ow <-
            spatstat.geom::owin(xrange = range(data$x), yrange = range(data$y))
        
        if (window == "convex") {
            p <- spatstat.geom::ppp(data$x, data$y, ow)
            ow <- spatstat.geom::convexhull(p)
            
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
                spatstat.geom::owin(
                    xrange = range(poly$x),
                    yrange = range(poly$y),
                    poly = poly
                )
            
        }
        ow
    }




#' @importFrom spatstat.geom ppp
pppGenerate <- function(cells, window, window.length) {
    ow <- makeWindow(cells, window, window.length)
    pppCell <- spatstat.geom::ppp(
        cells$x,
        cells$y,
        window = ow,
        marks = cells$cellType
    )
    
    pppCell
}



#' @importFrom spatstat.core density.ppp
#' @importFrom spatstat.geom closepairs nearest.valid.pixel area ppp
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join
inhomLPair <- 
    function (data,
              Rs = c(20, 50, 100),
              sigma = 10000,
              window = "convex",
              window.length = NULL,
              minLambda = 0.05, from = NULL, to = NULL) {
        ow <- makeWindow(data, window, window.length)
        X <-
            spatstat.geom::ppp(
                x = data$x,
                y = data$y,
                window = ow,
                marks = data$cellType
            )
        
        if (is.null(Rs))
            Rs = c(20, 50, 100)
        if (is.null(sigma))
            sigma = 100000
        
        maxR <- min(ow$xrange[2]- ow$xrange[1], ow$yrange[2]- ow$yrange[1])/2.01
        Rs <- unique(pmin(c(0, sort(Rs)),maxR))
        
        den <- spatstat.core::density.ppp(X, sigma = sigma)
        den <- den / mean(den)
        den$v <- pmax(den$v, minLambda)
        
        if(is.null(from)) from <- levels(data$cellType)
        if(is.null(to)) to <- levels(data$cellType)
        
        use <- data$cellType %in% c(from, to)
        data <- data[use,]
        X <- X[use,]
        
        
        p <- spatstat.geom::closepairs(X, max(Rs), what = "ijd")
        n <- X$n
        p$j <- data$cellID[p$j]
        p$i <- data$cellID[p$i]
        
        cT <- data$cellType
        names(cT) <- data$cellID
        
        p$d <- cut(p$d, Rs, labels = Rs[-1], include.lowest = TRUE)
        
        # inhom density
        np <- spatstat.geom::nearest.valid.pixel(X$x, X$y, den)
        w <- den$v[cbind(np$row, np$col)]
        names(w) <- data$cellID
        p$wt <- 1/w[p$j]*mean(w)
        rm(np)
        
        
        lam <- table(data$cellType)/spatstat.geom::area(X)
       # p$wt <- as.numeric(p$wt/lam[cT[p$j]])
        
   
        p$cellTypeJ <- cT[p$j]
        p$cellTypeI <- cT[p$i]
        p$i <- factor(p$i, levels = data$cellID)
        
        
        edge <- sapply(Rs[-1],function(x)borderEdge(X,x), simplify = FALSE)
        edge <- do.call("cbind",edge)
        edge <- as.data.frame(edge)
        colnames(edge) <- Rs[-1]
        edge$i <- data$cellID
        edge <- tidyr::pivot_longer(edge,-i,"d")
        
        p <- dplyr::left_join(as.data.frame(p), edge, c("i", "d"))
        p$d <- factor(p$d, levels = Rs[-1])
        
        use <- p$cellTypeI %in% from & p$cellTypeJ %in% to
        p <- p[use,]
    
        r <- inhomL(p, lam, X, Rs)
        wt <- r$wt
        names(wt) <- paste(r$cellTypeI, r$cellTypeJ, sep = "_")
        
        m1 <- rep(from, times = length(to))
        m2 <- rep(to, each = length(from))
        labels <- paste(m1, m2, sep = "_")
        
        assoc <- rep(NA, length(labels))
        names(assoc) <- labels
        assoc <- wt[labels]
        names(assoc) <- labels
        
        assoc
    }


#' @importFrom data.table as.data.table setkey CJ .SD ":="
#' @importFrom spatstat.geom area
inhomL <-
    function (p, lam, X, Rs) {
        r <- data.table::as.data.table(p)
        r$wt <- r$wt/r$value/as.numeric(lam[r$cellTypeJ])/as.numeric(lam[r$cellTypeI])/spatstat.geom::area(X)
        r <- r[,j:=NULL]
        r <- r[,value:=NULL]
        r <- r[,i:=NULL]
        data.table::setkey(r, d, cellTypeI, cellTypeJ)
        r <- r[data.table::CJ(d, cellTypeI, cellTypeJ, unique = TRUE)
        ][, lapply(.SD, sum), by = .(d, cellTypeI, cellTypeJ)
        ][is.na(wt), wt := 0]
        r <- r[, wt := cumsum(wt), by = list(cellTypeI, cellTypeJ)]
        r <- r[, list(wt=sum(sqrt(wt/pi))), by=.(cellTypeI, cellTypeJ)]
        r$wt <- r$wt - sum(Rs)
        
        r <- as.data.frame(r)
        
        r
        
    }




#' @importFrom spatstat.geom union.owin border inside.owin solapply intersect.owin area discs
borderEdge <- function(X, maxD){
    W <-X$window
    bW <- spatstat.geom::union.owin(spatstat.geom::border(W,maxD, outside = FALSE),
                               spatstat.geom::border(W,2, outside = TRUE))
    inB <- spatstat.geom::inside.owin(X$x, X$y, bW)
    e <- rep(1, X$n)
    if(any(inB)){
    circs <-spatstat.geom:: discs(X[inB], maxD, separate = TRUE)
    circs <- spatstat.geom::solapply(circs, spatstat.geom::intersect.owin, X$window)
    areas <- unlist(lapply(circs, spatstat.geom::area))/(pi*maxD^2)
    e[inB] <- areas
    }
    
    e
}
