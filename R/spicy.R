#' Performs spatial tests on spatial cytometry data.
#'
#' @param cells A segmentedCells or data frame that contains at least the variables x and y, giving the location of each cell, and cellType.
#' @param condition Vector of conditions to be tested corresponding to each image if cells is a data frame.
#' @param subject Vector of subject IDs corresponding to each image if cells is a data frame.
#' @param covariates Vector of covariate names that should be included in the mixed effects model as fixed effects.
#' @param from vector of cell types which you would like to compare to the to vector
#' @param to vector of cell types which you would like to compare to the from vector
#' @param dist The distance at which the statistic is obtained.
#' @param integrate Should the statistic be the integral from 0 to dist, or the value of the L curve at dist.
#' @param nsim Number of simulations to perform. If empty, the p-value from lmerTest is used.
#' @param verbose logical indicating whether to output messages.
#' @param ... Other options to pass to bootstrap.
#'
#' @return Data frame of p-values.
#' @export
#'
#' @examples
#' data("melanomaResponders")
#'
#' # Test with random effect for patient on only one pairwise combination of cell types.
#' spicy(melanomaResponders, condition = "condition", subject = "subject", from = "CD8+PD1+PDL1-", to = "CD8-PD1+PDL1+")
#' \dontrun{
#' # Test all pairwise combination of cell types without random effect of patient.
#' spicyTest <- spicy(melanomaResponders, condition = "condition")
#'
#'  Test all pairwise combination of cell types with random effect of patient.
#' #spicy(melanomaResponders, condition = "condition", subject = "subject")
#'
#' # Test all pairwise combination of cell types with random effect of patient using a bootstrap to calculate significance.
#' spicy(melanomaResponders, condition = "condition", subject = "subject", nsim = 199)
#' }
#' @aliases
#' spicy
#' spicy,spicy-method
#' @importFrom mgcv gam ti
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
                  ...) {
  if (!is(cells, "segmentedCells")) {
    stop('cells needs to be a segmentedCells object')
  }
  
  if (is.null(from))
    from <- as.character(unique(cellType(cells)))
  if (is.null(to))
    to <- as.character(unique(cellType(cells)))
  
  from <- as.character(unique(from))
  to <- as.character(unique(to))
  
  if (any((!to %in% cellType(cells)) |
          (!from %in% cellType(cells))))
    stop("to and from need to be cell type in your segmentedCells")
  
  nCells <- table(imageID(cells), cellType(cells))
  
  ## Find pairwise associations
  
  m1 <- rep(from, times = length(to))
  m2 <- rep(to, each = length(from))
  labels <- paste(m1, m2, sep = "_")
  
  MoreArgs1 <- list(cells = cells, dist = dist)
  
  if (verbose)
    message("Calculating pairwise spatial associations")
  
  pairwiseAssoc <- mapply(getPairwise,
                          from = m1,
                          to = m2,
                          MoreArgs = MoreArgs1)
  
  
  
  count1 <- as.vector(nCells[, m1])
  count2 <- as.vector(nCells[, m2])
  
  resSq <-
    as.vector(apply(pairwiseAssoc, 2, function(x)
      (x - mean(x, na.rm = TRUE))^2))
  
  weightFunction <- gam(resSq ~ ti(count1, count2))
  
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
        weightFunction = weightFunction
      )
    
    linearModels <- mapply(
      spatialLM,
      spatAssoc = pairwiseAssoc,
      from = m1,
      to = m2,
      MoreArgs = MoreArgs2,
      SIMPLIFY = FALSE
    )
    
    df <- cleanLM(linearModels)
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
        weightFunction = weightFunction
      )
    
    mixed.lmer <- mapply(
      spatialMEM,
      spatAssoc = pairwiseAssoc,
      from = m1,
      to = m2,
      MoreArgs = MoreArgs2,
      SIMPLIFY = FALSE
    )
    df <- cleanMEM(mixed.lmer, nsim)
    
  }
  
  
  
  df$pairwiseAssoc <- pairwiseAssoc
  df$comparisons <- data.frame(from = m1,
                               to = m2,
                               labels = labels)
  
  df <- new('spicy', df)
  df
}




cleanLM <- function(linearModels) {
  tLm <- lapply(linearModels, function(LM) {
    coef <- as.data.frame(t(summary(LM)$coef))
    coef <-
      split(coef, c("coefficient", "se", "statistic", "p.value"))
  })
  
  df <- apply(do.call(rbind, tLm), 2, function(x)
    do.call(rbind, x))
  df <- lapply(df, function(x) {
    rownames(x) <- names(linearModels)
    x
  })
  df
}









cleanMEM <- function(mixed.lmer, nsim) {
  if (length(nsim) > 0) {
    boot <- lapply(mixed.lmer, spatialMEMBootstrap, nsim = nsim)
    #p <- do.call(rbind, p)
    tBoot <- lapply(boot, function(coef) {
      coef <- as.data.frame(t(coef))
      coef <- split(coef, c("coefficient", "se", "p.value"))
    })
    
    
    df <- apply(do.call(rbind, tBoot), 2, function(x)
      do.call(rbind, x))
    df <- lapply(df, function(x) {
      rownames(x) <- names(mixed.lmer)
      x
    })
    
  } else {
    tLmer <- lapply(mixed.lmer, function(lmer) {
      coef <- as.data.frame(t(summary(lmer)$coef))
      coef <-
        split(coef,
              c("coefficient", "se", "df", "statistic", "p.value"))
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
#' @param cells A segmentedCells or data frame that contains at least the variables x and y, giving the location of each cell, and cellType.
#' @param from The 'from' cellType for generating the L curve.
#' @param to The 'to' cellType for generating the L curve.
#' @param dist The distance at which the statistic is obtained.
#'
#' @return Statistic from pairwise L curve of a single image.
#'
#'
#' @examples
#' data("melanomaResponders")
#' pairAssoc <- getPairwise(melanomaResponders)
#' @export
getPairwise <- function(cells, from, to, dist = NULL) {
  cells <- location(cells, bind = FALSE)
  
  pairwiseVals <- lapply(cells,
                         getStat,
                         from = from,
                         to = to,
                         dist = dist)
  
  unlist(pairwiseVals)
}


#' @importFrom spatstat Lcross
getStat <- function(cells, from, to, dist) {
  pppCell <- pppGenerate(cells)
  
  L <- tryCatch({
    Lcross(pppCell,
           from = from,
           to = to,
           correction = "best")
  }, error = function(e) {
    
  })
  
  if (length(class(L)) == 1) {
    return(NA)
  }
  
  if (is.null(dist))
    dist <- max(L$r)
  
  theo = L$theo[L$r <= dist]
  iso = L$iso[L$r <= dist]
  mean(iso - theo)
}



# Performs bootstrapping to estimate p-value.

#' @importFrom lme4 fixef bootMer
spatialMEMBootstrap <- function(mixed.lmer, nsim = 19) {
  bootCoef <- bootMer(mixed.lmer,
                      fixef,
                      nsim = nsim,
                      re.form = NA)
  
  stats <- bootCoef$t
  fe <- fixef(mixed.lmer)
  pval = pmin(colMeans(stats < 0), colMeans(stats > 0)) * 2
  df <-
    data.frame(
      coefficient = fe,
      se = apply(stats, 2, stats::sd),
      p.value = pval
    )
  df
}







.show_spicy <- function(df) {
  pval <- as.data.frame(df$p.value)
  cond <- colnames(pval)[grep('condition', colnames(pval))]
  cat(df$test)
  cat("Number of cell type pairs: ", nrow(pval), "\n")
  cat("Number of differentially localised cell type pairs: \n")
  if (nrow(pval) == 1)
    print(sum(pval[cond] < 0.05))
  if (nrow(pval) > 1)
    print(colSums(apply(pval[cond], 2, p.adjust, 'fdr') < 0.05))
  
}
setMethod("show", signature(object = "spicy"), function(object) {
  .show_spicy(object)
})





#' @importFrom lmerTest lmer
spatialMEM <-
  function(spatAssoc,
           from,
           to,
           cells,
           subject,
           condition,
           covariates,
           weightFunction) {
    cellCounts <- table(imageID(cells), cellType(cells))
    
    count1 <- cellCounts[, from]
    count2 <- cellCounts[, to]
    filter <- !is.na(spatAssoc)
    
    if (sum(filter) < 3)
      return(NA)
    
    pheno <- as.data.frame(phenotype(cells))
    spatialData <-
      data.frame(spatAssoc,
                 condition = pheno[, condition],
                 subject = pheno[, subject],
                 pheno[covariates])
    
    spatialData <- spatialData[filter, ]
    count1 <- count1[filter]
    count2 <- count2[filter]
    
    if (is.null(weightFunction)) {
      w <- rep(1, length(count1))
    } else{
      z1 <- predict(weightFunction, data.frame(count1, count2))
      w <- 1 / sqrt(z1 - min(z1) + 1)
      w <- w / sum(w)
    }
    
    formula <- 'spatAssoc ~ condition + (1|subject)'
    if (!is.null(covariates))
      formula <-
      paste('spatAssoc ~ condition + (1|subject)',
            paste(covariates, collapse = '+'),
            sep = "+")
    mixed.lmer <- lmer(formula(formula),
                       data = spatialData,
                       weights = w)
    mixed.lmer
  }




spatialLM <-
  function(spatAssoc,
           from,
           to,
           cells,
           condition,
           covariates,
           weightFunction) {
    cellCounts <- table(imageID(cells), cellType(cells))
    
    count1 <- cellCounts[, from]
    count2 <- cellCounts[, to]
    filter <- !is.na(spatAssoc)
    
    if (sum(filter) < 3)
      return(NA)
    
    pheno <- as.data.frame(phenotype(cells))
    spatialData <-
      data.frame(spatAssoc, condition = pheno[, condition], pheno[covariates])
    
    spatialData <- spatialData[filter, ]
    count1 <- count1[filter]
    count2 <- count2[filter]
    
    if (is.null(weightFunction)) {
      w <- rep(1, length(count1))
    } else{
      z1 <- predict(weightFunction, data.frame(count1, count2))
      w <- 1 / sqrt(z1 - min(z1) + 1)
      w <- w / sum(w)
    }
    
    formula <- 'spatAssoc ~ condition'
    if (!is.null(covariates))
      formula <-
      paste('spatAssoc ~ condition',
            paste(covariates, collapse = '+'),
            sep = "+")
    lm1 <- lm(formula(formula),
              data = spatialData,
              weights = w)
    lm1
  }


#' Plots result of signifPlot.
#'
#' @param results Data frame obtained from spicy.
#' @param fdr TRUE if FDR correction is used.
#' @param breaks Vector of 3 numbers giving breaks used in pheatmap. The first number is the minimum, the second is the maximum, the third is the number of breaks.
#' @param col Vector of colours to use in pheatmap.
#'
#' @return a pheatmap object
#'
#' @examples
#' \dontrun{
#' example(spicy)
#' signifPlot(spicyTest, breaks=c(-3, 3, 0.5))
#' }
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
signifPlot <- function(results,
                       fdr = FALSE,
                       breaks = c(-5, 5, 0.5),
                       col = c("blue", "white", "red")) {
  pVal <- results$p.value$conditionResponders
  marks <- unique(results$comparisons$from)
  
  if (min(pVal) == 0) {
    pVal[pVal == 0] <-
      pVal[pVal == 0] + 10 ^ floor(log10(min(pVal[pVal > 0])))
  }
  
  if (fdr) {
    pVal <- p.adjust(pVal, method = "fdr")
  }
  
  isGreater <- results$coefficient$conditionResponders > 0
  
  pVal <- log10(pVal)
  
  pVal[isGreater] <- abs(pVal[isGreater])
  
  pVal <- matrix(pVal, nrow = length(marks))
  colnames(pVal) <- marks
  rownames(pVal) <- marks
  
  breaks <- seq(from = breaks[1], to = breaks[2], by = breaks[3])
  pal <- grDevices::colorRampPalette(col)(length(breaks))
  
  heatmap <- pheatmap(
    pVal,
    col = pal,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
  
  heatmap
}



#' @importFrom spatstat ppp
pppGenerate <- function(cells) {
  pppCell <- ppp(
    cells$x,
    cells$y,
    xrange = c(0, max(cells$x)),
    yrange = c(0, max(cells$y)),
    marks = cells$cellType
  )
  
  pppCell
}
