## The original spicyR test (method = "image"): a per-image L-function summary per pair, compared between
## conditions by a weighted linear or mixed model (Canete et al. 2022). Called by spicy(); not exported.
#' @noRd
.spicyImage <- function(cells,
                  condition,
                  subject = NULL,
                  covariates = NULL,
                  imageID = "imageID",
                  cellType = "cellType",
                  spatialCoords = c("x", "y"),
                  r = NULL,
                  sigma = NULL,
                  from = NULL,
                  to = NULL,
                  alternateResult = NULL,
                  cores = 1,
                  minLambda = 0.05,
                  weights = TRUE,
                  weightsByPair = FALSE,
                  weightFactor = 1,
                  weightZThreshold = 0.1,
                  window = "convex",
                  window.length = NULL,
                  edgeCorrect = TRUE,
                  includeZeroCells = FALSE,
                  verbose = FALSE,
                  BPPARAM = NULL,
                  imageIDCol = imageID,
                  cellTypeCol = cellType,
                  spatialCoordCols = spatialCoords,
                  nCores = cores,
                  Rs = r,
                  ...) {
  
  user_args = as.list(match.call())[-1]
  user_vals = lapply(user_args, eval, envir = parent.frame())
  # spicy() always passes imageID, cellType, spatialCoords, r and cores on; a deprecated argument given
  # with its replacement at the default takes over.
  defaults <- list(imageID = "imageID", cellType = "cellType", spatialCoords = c("x", "y"), r = NULL, cores = 1)
  deprecated <- list(imageID = "imageIDCol", cellType = "cellTypeCol", spatialCoords = "spatialCoordCols",
                     r = "Rs", cores = c("nCores", "BPPARAM"))
  for (arg in names(deprecated)) {
    if (any(deprecated[[arg]] %in% names(user_vals)) && arg %in% names(user_vals) &&
        identical(user_vals[[arg]], defaults[[arg]])) {
      user_vals[arg] <- NULL
    }
  }
  argumentChecks("spicy", user_vals)

  # Workers for the per-pair weights and survival models (BiocParallel params are accepted for
  # backward compatibility and converted to a worker count).
  cores <- .n_workers(cores, BPPARAM)

  if (.is_class(cells, "SummarizedExperiment") || is.data.frame(cells)) {
    cells <- .format_data(
      cells, imageID, cellType, spatialCoords, verbose
    )
  }

  if (is.null(from) || is.null(to)) {
    if (is.null(from)) {
      from <- as.character(unique(getCellType(cells)))
    }
    if (is.null(to)) {
      to <- as.character(unique(getCellType(cells)))
    }

    m1 <- rep(from, times = length(to))
    m2 <- rep(to, each = length(from))
    labels <- paste(m1, m2, sep = "__")
  } else {
    # from and to are paired element by element; a length-one side is recycled (spicyR 1.x failed in the
    # weight fit when, e.g., one `from` met several `to`)
    n_pair <- max(length(from), length(to))
    m1 <- rep_len(from, n_pair)
    m2 <- rep_len(to, n_pair)
    labels <- paste(m1, m2, sep = "__")
    if (any(duplicated(labels))) stop("There are duplicated from-to pairs")
  }


  if (any((!to %in% getCellType(cells)) | (!from %in% getCellType(cells)))) {
    stop("to and from need to be cell type in your data")
  }

  nCells <- table(getImageID(cells), getCellType(cells))
  
  
  if (!is.null(condition)) {
    conditionVector <- as.data.frame(getImagePheno(cells))[condition][, 1]
    
    if (!inherits(conditionVector, "Surv")) {
      wasFactor <- is.factor(conditionVector)
      
      if (!wasFactor) {
        conditionVector <- as.factor(conditionVector)
      }
      
      conditionVector <- droplevels(conditionVector)
      conditionVector <- stats::relevel(conditionVector, ref = levels(conditionVector)[1])
      
      if (!wasFactor || TRUE) {  
        message(
          paste0(
            if (!wasFactor) "Coercing condition into factor. " else "",
            "Dropping unused levels. Using ",
            condition, " = ", levels(conditionVector)[1],
            " as base comparison group. If this is not the desired base group,",
            " please convert cells$", condition, " into a factor and change the order of levels(cells$",
            condition, ") so that the base group is at index 1."
          )
        )
      }
    }
  }
  
  

  ## Check whether the subject parameter has a one-to-one mapping with image
  if (!is.null(subject)) {
    if (nrow(as.data.frame(unique(cells[, subject]))) == nrow(as.data.frame(unique(cells[, imageID])))) {
      subject <- NULL
      
      if(inherits(conditionVector, "Surv")) {
        warning("Your specified subject parameter has a one-to-one mapping with imageID. Converting to a coxph model instead of cox mixed effects model.")
      } else{
        warning("Your specified subject parameter has a one-to-one mapping with imageID. Converting to a linear model instead of mixed model.")
      }
      
    }  
  }

  
  spicyResult = list()

  
  comparisons <- data.frame(from = m1,
                            to = m2,
                            labels = labels)

  ## Find pairwise associations

  if (is.null(alternateResult)) {
    pairwiseAssoc <- getPairwise(cells,
        r = r,
        sigma = sigma,
        from = from,
        to = to,
        minLambda = minLambda,
        window = window,
        window.length = window.length,
        edgeCorrect = edgeCorrect,
        includeZeroCells = includeZeroCells,
        cores = cores
    )
    pairwiseAssoc <- as.data.frame(pairwiseAssoc)
    pairwiseAssoc <- pairwiseAssoc[labels]
    
  }
  
  
  if (!is.null(alternateResult)) {
    pairwiseAssoc <- alternateResult
    
    # Checking if Kontextual result.
    if (isTRUE(attr(alternateResult, "kontextualResult"))) {
      pairwiseAssoc <- alternateResult
      
      labels <- names(pairwiseAssoc)
      
      comparisons <- .split_labels(labels, c("from", "to", "parent"))
      comparisons$labels <- paste(comparisons$from, comparisons$to, comparisons$parent, sep = "__")
      
      m1 <- comparisons$from
      m2 <- comparisons$to
      
      spicyResult$isKontextual = TRUE
    
    }
  }
  
  
  weightFunction <- getWeightFunction(
    pairwiseAssoc, nCells, m1, m2, cores, weights, weightsByPair, weightFactor,
    weightZThreshold
  )
  
  # Matrix needed for survival analysis
  pairwiseAssocMatrix = pairwiseAssoc
  
  pairwiseAssoc <- as.list(data.frame(pairwiseAssoc))
  names(pairwiseAssoc) <- labels
  

  
  ## Survival model
  if(inherits(conditionVector, "Surv")){
    
    survivalResult = spatialSurv(
      measurementMat = pairwiseAssocMatrix,
      condition = condition,
      pheno = as.data.frame(getImagePheno(cells)),
      covariates = covariates,
      subject = subject,
      weights = weightFunction,
      BPPARAM = cores
    )
    
    spicyResult$survivalOutcome = conditionVector
    spicyResult = append(spicyResult, list("survivalResults" = survivalResult))
    
  }


  ## Linear model
  if (!inherits(conditionVector, "Surv") && is.null(subject) && !is.null(condition)) {
    if (verbose) {
      message("Testing for spatial differences across conditions")
    }

    MoreArgs2 <-
      list(
        cells = cells,
        condition = condition,
        covariates = covariates,
        cellCounts = table(getImageID(cells), getCellType(cells)),
        pheno = as.data.frame(getImagePheno(cells))
      )


    linearModels <- mapply(
      spatialLM,
      spatAssoc = pairwiseAssoc,
      from = m1,
      to = m2,
      weightFunction = weightFunction,
      MoreArgs = MoreArgs2,
      SIMPLIFY = FALSE
    )

    lmResult <- cleanLM(linearModels)
    spicyResult = append(spicyResult, lmResult)
  }


  ## Mixed effects model
  if (!inherits(conditionVector, "Surv") && (!is.null(subject)) && !is.null(condition)) {
    if (verbose) {
      message(
        "Testing for spatial differences across conditions accounting for multiple images per subject" # nolint
      )
    }

    MoreArgs2 <-
      list(
        cells = cells,
        subject = subject,
        condition = condition,
        covariates = covariates,
        cellCounts = table(getImageID(cells), getCellType(cells)),
        pheno = as.data.frame(getImagePheno(cells))
      )


    mixed.lmer <- mapply(
      spatialMEM,
      spatAssoc = pairwiseAssoc,
      from = m1,
      to = m2,
      weightFunction = weightFunction,
      MoreArgs = MoreArgs2,
      SIMPLIFY = FALSE
    )


    melmResult <- cleanMEM(mixed.lmer)
    spicyResult = append(spicyResult, melmResult)
  }

  
  if(!is.null(condition)){
    spicyResult$condition <- conditionVector
  }

  if (!is.null(subject)) {
    spicyResult$subject <- as.data.frame(getImagePheno(cells))[subject][, 1]
  }
  
  spicyResult$pairwiseAssoc <- pairwiseAssoc
  spicyResult$comparisons <- comparisons
  
  spicyResult$weights <- weightFunction
  spicyResult$nCells <- nCells
  
  spicyResult$imageIDs <- as.data.frame(getImagePheno(cells))["imageID"][, 1]
  spicyResult$alternateResult <- ifelse(is.null(alternateResult), FALSE, TRUE)
  
  spicyResult <- methods::new("SpicyResults", spicyResult)
  spicyResult
}




cleanLM <- function(linearModels, BPPARAM = NULL) {
  tLm <- lapply(linearModels, function(LM) {
    if (is(LM, "lm")) {
      coef <- as.data.frame(t(summary(LM)$coef))
      coef <-
        split(coef, c("coefficient", "se", "statistic", "p.value"))
    } else {
      n <- data.frame(NA)
      colnames(n) <- "(Intercept)"
      coef <- list(coefficient = n, se = n, statistic = n, p.value = n)
    }
    coef
  })


  df <- do.call("rbind", tLm)

  df <- suppressWarnings(apply(df, 2, function(x) {
    .bind_rows(x)
  }))

  df <- lapply(df, function(x) {
    rownames(x) <- names(linearModels)
    x
  })
  df
}


cleanMEM <- function(mixed.lmer, BPPARAM = NULL) {
  tLmer <- lapply(mixed.lmer, function(lmer) {
    if (is.matrix(lmer)) {
      coef <- as.data.frame(t(lmer))
      coef <-
        split(
          coef,
          c("coefficient", "se", "df", "statistic", "p.value")
        )
    } else {
      n <- data.frame(NA)
      colnames(n) <- "(Intercept)"
      coef <- list(
        coefficient = n, df = n, p.value = n, se = n, statistic = n
      )
    }
    coef
  })

  df <- do.call("rbind", tLmer)

  df <- suppressWarnings(apply(df, 2, function(x) {
    .bind_rows(x)
  }))

  df <- lapply(df, function(x) {
    rownames(x) <- names(mixed.lmer)
    x
  })
  df
}

#' Get statistic from pairwise L curve of a single image.
#'
#' @param cells A SummarizedExperiment that contains at least the
#'     variables x and y, giving the location coordinates of each cell, and
#'     cellType.
#' @param imageID
#'     The name of the imageID column if using a SingleCellExperiment or SpatialExperiment.
#' @param cellType The name of the cellType column if using a SingleCellExperiment or SpatialExperiment.
#' @param spatialCoords The names of the spatialCoords column if using a SingleCellExperiment.
#' @param r A vector of the radii that the measures of association should be calculated over.
#' @param sigma A numeric variable used for scaling when fitting inhomogenous L-curves.
#' @param from The 'from' cellType for generating the L curve.
#' @param to The 'to' cellType for generating the L curve.
#' @param cores Number of cores to use for parallel processing or a BiocParallel MulticoreParam or SerialParam object.
#' @param minLambda Minimum value density for scaling when fitting inhomogeneous L-curves.
#' @param window Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length A tuning parameter for controlling the level of concavity when estimating concave windows.
#' @param edgeCorrect A logical indicating whether to perform edge correction.
#' @param includeZeroCells A logical indicating whether to include cells with zero counts in the pairwise association.
#' @param BPPARAM \{DEPRECATED\} A BiocParallel MulticoreParam or SerialParam object. 
#' @param imageIDCol \{DEPRECATED\} The name of the imageID column if using a SingleCellExperiment or SpatialExperiment.
#' @param cellTypeCol \{DEPRECATED\} The name of the cellType column if using a SingleCellExperiment or SpatialExperiment.
#' @param spatialCoordCols \{DEPRECATED\} The names of the spatialCoords column if using a SingleCellExperiment.
#' @param nCores \{DEPRECATED\} Number of cores to use for parallel processing or a BiocParallel MulticoreParam or SerialParam object.
#' @param Rs \{DEPRECATED\} A vector of the radii that the measures of association should be calculated over.
#' calculation.
#' @return Statistic from pairwise L-curve of a single image.
#' @examples
#' data("diabetesData")
#' # Subset by imageID for fast example
#' selected_cells <- diabetesData[
#'   , SummarizedExperiment::colData(diabetesData)$imageID == "A09"
#' ]
#' pairAssoc <- getPairwise(selected_cells)
#' @export
getPairwise <- function(
    cells,
    imageID = "imageID",
    cellType = "cellType",
    spatialCoords = c("x", "y"),
    r = NULL,
    sigma = NULL,
    from = NULL,
    to = NULL,
    cores = 1, 
    minLambda = 0.05,
    window = "convex",
    window.length = NULL,
    edgeCorrect = TRUE,
    includeZeroCells = FALSE,
    BPPARAM = NULL,
    imageIDCol = imageID,
    cellTypeCol = cellType,
    spatialCoordCols = spatialCoords,
    nCores = cores,
    Rs = r
    ) {
    
    
    user_args = as.list(match.call())[-1]
      
    tryCatch({
      user_vals = lapply(user_args, eval, envir = parent.frame())
      argumentChecks("getPairwise", user_vals)
    }, error = function(e) {
    if (grepl("object 'cells' not found", e$message)) {
        message("Skipping argument checks as `getPairwise()` is being called within `spicy()`")
    } else {
      stop(e)
    }
  })

  
  if (.is_class(cells, "SummarizedExperiment")) {
    cells <- .format_data(
      cells, imageID, cellType, spatialCoords, FALSE
    )
  }

  nThreads <- .n_workers(cores, BPPARAM)

  # Square and convex windows without density weighting: one threaded C++
  # call for all images, with `cores` threads.
  lev <- levels(cells$cellType)
  if (is.null(sigma) && window %in% c("square", "convex") && !is.null(lev) &&
      all(c(from, to) %in% lev)) {
    return(getPairwiseThreaded(
      cells, Rs, from, to, window, edgeCorrect, includeZeroCells, nThreads
    ))
  }

  # Inhomogeneous (sigma) or concave windows: per image in R, with spatstat.
  .need("spatstat.geom", "for `sigma` (inhomogeneous L) or concave windows")
  if (!is.null(sigma)) .need("spatstat.explore", "for `sigma` (inhomogeneous L)")
  if (window == "concave") .need("concaveman", "for concave windows")

  cells2 <- getCellSummary(cells, bind = FALSE)


  if (is.null(from)) from <- levels(cells2$cellType)
  if (is.null(to)) to <- levels(cells2$cellType)

  pairwiseVals <- .par_lapply(cells2,
    inhomLPair,
    Rs = Rs,
    sigma = sigma,
    window = window,
    window.length = window.length,
    minLambda = minLambda,
    from = from,
    to = to,
    edgeCorrect = edgeCorrect,
    includeZeroCells = includeZeroCells,
    cores = nThreads
  )
  return(do.call("rbind", pairwiseVals))

  unlist(pairwiseVals)
}




# getPairwise() for square and convex windows without density weighting: the
# inhomLPair() computation for every image in one C++ call, with the images
# spread over nThreads threads (src/pairwise.cpp).
getPairwiseThreaded <- function(cells, Rs, from, to, window, edgeCorrect,
                                includeZeroCells, nThreads) {
  img <- droplevels(as.factor(cells$imageID))
  o <- order(as.integer(img))
  lev <- levels(cells$cellType)
  if (is.null(Rs)) Rs <- c(20, 50, 100)
  if (is.null(from)) from <- lev
  if (is.null(to)) to <- lev
  m1 <- rep(from, times = length(to))
  m2 <- rep(to, each = length(from))
  res <- getPairwiseCpp(
    as.numeric(cells$x[o]), as.numeric(cells$y[o]), as.integer(cells$cellType)[o],
    c(0L, cumsum(tabulate(as.integer(img), nlevels(img)))), length(lev), as.numeric(Rs),
    window == "square", lev %in% from, lev %in% to, match(m1, lev), match(m2, lev),
    edgeCorrect, includeZeroCells, as.integer(nThreads)
  )
  dimnames(res) <- list(levels(img), paste(m1, m2, sep = "__"))
  res
}


#' Get proportions from a SummarizedExperiment.
#'
#' @param cells A SingleCellExperiment, SpatialExperiment or data.frame.
#' @param feature The feature of interest
#' @param imageID The imageID's
#'
#' @return Proportions
#'
#'
#' @examples
#' data("diabetesData")
#' prop <- getProp(diabetesData)
#' @export
getProp <- function(cells, feature = "cellType", imageID = "imageID") {
  if (is.data.frame(cells)) {
    df <- cells[, c(imageID, feature)]
  }

  if (.is_class(cells, "SummarizedExperiment")) {
    .need("SummarizedExperiment", "to use SummarizedExperiment-based inputs")
    df <- as.data.frame(
      SummarizedExperiment::colData(cells)
    )[, c(imageID, feature)]
  }


  tab <- table(df[, imageID], df[, feature])
  tab <- sweep(tab, 1, rowSums(tab), "/")
  as.data.frame.matrix(tab)
}


#' @importFrom stats p.adjust
.show_SpicyResults <- function(df) {
  if (identical(df$method, "cell")) {
    tab <- df$cellResults
    what <- if (!is.null(df$survivalOutcome)) "association with survival" else
      paste0(levels(df$condition)[-1], " vs ", levels(df$condition)[1], collapse = "; ")
    scale <- if (!is.null(df$k)) paste0("k = ", df$k, " nearest neighbours") else paste0("r = ", paste(df$r, collapse = ", "))
    cat("spicyR (cell-level test): ", length(unique(paste(tab$from, tab$to))), " pairs, ", what, ", ", scale, "\n", sep = "")
    if (is.null(df$subject)) cat("Units: ", length(df$imageID), " images (no subject given: each image is a patient)\n", sep = "")
    else cat("Units: ", length(unique(df$subject)), " patients with ", length(df$imageID), " images\n", sep = "")
    cat("BH-adjusted p < 0.05: ", sum(tab$p_adj < 0.05, na.rm = TRUE), " pairs", sep = "")
    if (!is.null(tab$adjusted_p_adj)) cat(" (", sum(tab$adjusted_p_adj < 0.05, na.rm = TRUE), " after adjusting for abundance)", sep = "")
    cat("\nSee topPairs() and $cellResults.\n")
    return(invisible(df))
  }
  pval <- as.data.frame(df$p.value)
  cond <- colnames(pval)[grep("condition", colnames(pval))]
  cat("spicyR (image-level test): ", nrow(pval), " pairs\n", sep = "")
  cat("Pairs with BH-adjusted p < 0.05:\n")
  if (nrow(pval) == 1) print(sum(pval[cond] < 0.05, na.rm = TRUE))
  if (nrow(pval) > 1) print(colSums(apply(pval[cond], 2, p.adjust, "fdr") < 0.05, na.rm = TRUE))
  invisible(df)
}
setMethod(
  "show", methods::signature(object = "SpicyResults"), function(object) {
    .show_SpicyResults(object)
  }
)

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

    # count1 <- cellCounts[, from]
    # count2 <- cellCounts[, to]

    spatialData <-
      data.frame(
        "spatAssoc" = spatAssoc,
        "condition" = pheno[, condition],
        "subject" = pheno[, subject],
        "covariates" = pheno[covariates]
      )

    names(spatialData)[grep("covariates", names(spatialData))] <- covariates

    formula <- "spatAssoc ~ condition + (1|subject)"

    if (!is.null(covariates)) {
      formula <-
        paste("spatAssoc ~ condition + (1|subject)",
          paste(covariates, collapse = "+"),
          sep = "+"
        )
    }

    spatialData$weights <- weightFunction

    # Fit in C++ (lmerCoefTable()), on the rows and design lme4 would use. Cases
    # it does not cover fall through to lmerTest below.
    tab <- tryCatch({
      fd <- droplevels(spatialData[stats::complete.cases(spatialData), , drop = FALSE])
      X <- stats::model.matrix(stats::formula(sub(" \\+ \\(1\\|subject\\)", "", formula)), fd)
      nSubject <- length(unique(fd$subject))
      if (nSubject >= 2 && nSubject < nrow(fd) && all(fd$weights > 0) &&
          qr(X)$rank == ncol(X)) {
        lmerCoefTable(X, fd$spatAssoc, fd$weights, fd$subject)
      }
    }, error = function(e) NULL)
    if (!is.null(tab)) return(tab)

    # Fallback for the cases the C++ fit does not cover.
    .need("lmerTest", "for this mixed model (the built-in fit could not be used)")
    mixed.lmer <- suppressWarnings(suppressMessages(tryCatch(
      {
        lmerTest::lmer(stats::formula(formula),
          data = spatialData,
          weights = spatialData$weights # TODO: weights does not exist or is a funciton.
        )
      },
      error = function(e) {

      }
    )))
    if (!is(mixed.lmer, "lmerMod")) {
      return(NA)
    }
    if (is(mixed.lmer, "lmerModLmerTest") && mixed.lmer@devcomp$cmp["REML"] == -Inf) {
      return(NA)
    }


    summary(mixed.lmer)$coef
  }

# Coefficient table of lmerTest::lmer(y ~ X + (1 | subject), weights = w), as
# summary() gives it with Satterthwaite degrees of freedom. The REML fit and
# the exact derivatives lmerTest approximates with numDeriv come from C++
# (src/lmerRI.cpp); the rest follows lmerTest's as_lmerModLT() and contest1D().
# Returns NULL if the C++ fit is not possible.
#' @importFrom stats pt
lmerCoefTable <- function(X, y, w, subject) {
  g <- as.integer(factor(subject))
  fit <- lmerRandomIntercept(X, y, w, g, max(g))
  if (!fit$ok) return(NULL)
  eh <- eigen(fit$hessian, symmetric = TRUE)
  pos <- eh$values > 1e-8
  vcovVarpar <- 2 * eh$vectors[, pos, drop = FALSE] %*%
    diag(1 / eh$values[pos], nrow = sum(pos)) %*% t(eh$vectors[, pos, drop = FALSE])
  varCon <- diag(fit$vcov)
  grad <- cbind(diag(fit$dvcov_theta), diag(fit$dvcov_sigma))
  df <- 2 * varCon^2 / rowSums((grad %*% vcovVarpar) * grad)
  se <- sqrt(varCon)
  tval <- fit$beta / se
  tab <- cbind(fit$beta, se, df, tval, 2 * stats::pt(abs(tval), df, lower.tail = FALSE))
  dimnames(tab) <- list(colnames(X), c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)"))
  tab
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


    spatialData <-
      data.frame(
        "spatAssoc" = spatAssoc,
        "condition" = pheno[, condition],
        "covariates" = pheno[, covariates]
      )

    names(spatialData)[grep("covariates", names(spatialData))] <- covariates

    formula <- "spatAssoc ~ condition"

    if (!is.null(covariates)) {
      formula <-
        paste("spatAssoc ~ condition",
          paste(covariates, collapse = "+"),
          sep = "+"
        )
    }

    lm1 <- tryCatch(
      {
        stats::lm(stats::formula(formula),
          data = spatialData,
          weights = weightFunction # TODO: check that this works correctly
        )
      },
      error = function(e) {

      }
    )

    if (!is(lm1, "lm")) {
      return(NA)
    }


    lm1
  }

#' @importFrom survival coxph
spatialSurv <- function(measurementMat,
                        condition,
                        pheno,
                        covariates = NULL,
                        subject = NULL,
                        weights = NULL,
                        remove = NULL,
                        BPPARAM = NULL) {
  
  if (!is.null(subject)) .need("coxme", "for survival models with a `subject`")
  result <- .par_lapply(colnames(measurementMat), function(test) {
    measurementCol <- measurementMat[, test]
    ind <- !(measurementCol %in% remove)
    
    spatialData = data.frame("measurementCol" = measurementCol, "surv" = pheno[, condition])
    
    formula = "surv ~ measurementCol"
    
    if (!is.null(subject)) {
      spatialData = data.frame(spatialData, "subject" = pheno[, subject])
      
      formula <- paste(formula, "(1 | subject)", sep = "+")
    }
    
    if (!is.null(covariates)) {
      spatialData = data.frame(spatialData, "covariates" = pheno[, covariates])
      names(spatialData)[grep("covariates", names(spatialData))] <- covariates
      
      formula <- paste(formula, paste(covariates, collapse = "+"), sep = "+")
    }
    
    formula = stats::formula(formula)
    
    if(is.null(weights)) {
      spatialData$weights = 1
    } else {
      spatialData$weights = weights[[test]]
    }
    
    if(is.null(subject)) {
    # If no random effects must use coxph
    fit <- survival::coxph(formula, data = spatialData, weights = spatialData$weights, subset = ind)
    result <- summary(fit)$coefficients["measurementCol", c("coef", "se(coef)", "Pr(>|z|)")]
  } else{
    # Drop factors for any factor columns in spatialData
    isFactor <- vapply(spatialData, is.factor, logical(1))
    spatialData[isFactor] <- lapply(spatialData[isFactor], droplevels)
    
    # Mixed effects survival
    fit <- coxme::coxme(formula, data = spatialData, weights = spatialData$weights, subset = ind)
    
    # from print.coxme() function
    beta <- unname(fit$coefficients)
    nvar <- length(beta)
    nfrail <- nrow(fit$var) - nvar
    
    se <- sqrt(diag(as.matrix(fit$var))[nfrail + 1:nvar])
    p_val <- 1 - stats::pchisq((beta / se) ^ 2, 1)
    
    # Need to select first indexes of each value incase of multiple betas
    result <- c("coef" = beta[1], "se(coef)" = se[1], "Pr(>|z|)" = p_val[1])
  }
  
  return(result)
  }, cores = .n_workers(BPPARAM))

result <- .bind_rows(result)
result <- data.frame(
  test = colnames(measurementMat),
  coef = result[["coef"]],
  se.coef = result[["se(coef)"]],
  p.value = result[["Pr(>|z|)"]]
)
result <- result[order(result$p.value), , drop = FALSE]
rownames(result) <- NULL

return(result)
}


###########################
#
#  Generate distances
#
###########################


makeWindow <-
  function(data,
           window = "square",
           window.length = NULL) {
    data <- data.frame(data)
    ow <-
      spatstat.geom::owin(xrange = range(data$x), yrange = range(data$y))

    if (window == "convex") {
      p <- spatstat.geom::ppp(data$x, data$y, ow)
      ow <- spatstat.geom::convexhull(p)
    }
    if (window == "concave") {
      message("Concave windows are temperamental. Try choosing values of window.length > and < 1 if you have problems.") # nolint
      if (is.null(window.length)) {
        window.length <- (max(data$x) - min(data$x)) / 20
      } else {
        window.length <- (max(data$x) - min(data$x)) / 20 * window.length
      }
      dist <- (max(data$x) - min(data$x)) / (length(data$x))
      bigDat <-
        do.call(
          "rbind",
          lapply(as.list(as.data.frame(t(data[, c("x", "y")]))), function(x) {
            cbind(
              x[1] + c(0, 1, 0, -1, -1, 0, 1, -1, 1) * dist,
              x[2] + c(0, 1, 1, 1, -1, -1, -1, 0, 0) * dist
            )
          })
        )
      ch <-
        concaveman::concaveman(bigDat,
          length_threshold = window.length,
          concavity = 1
        )
      poly <- as.data.frame(ch[nrow(ch):1, ]) # nolint
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




inhomLPair <- function(data,
                       Rs = c(20, 50, 100),
                       sigma = NULL,
                       window = "convex",
                       window.length = NULL,
                       minLambda = 0.05,
                       from = NULL,
                       to = NULL,
                       edgeCorrect = TRUE,
                       includeZeroCells = TRUE) {
  ow <- makeWindow(data, window, window.length)


  X <-
    spatstat.geom::ppp(
      x = data$x,
      y = data$y,
      window = ow,
      marks = data$cellType,
      check = FALSE
    )

  if (is.null(Rs)) {
    Rs <- c(20, 50, 100)
  }

  maxR <- min(ow$xrange[2] - ow$xrange[1], ow$yrange[2] - ow$yrange[1]) / 2.01
  Rs <- unique(pmin(c(0, sort(Rs)), maxR))

  if (!is.null(sigma)) {
    den <- spatstat.explore::density.ppp(X, sigma = sigma)
    den <- den / mean(den)
    den$v <- pmax(den$v, minLambda)
  }

  if (is.null(from)) from <- levels(data$cellType)
  if (is.null(to)) to <- levels(data$cellType)

  use <- data$cellType %in% c(from, to)
  if (all(!use)) {
    return(NA)
  }
  data <- data[use, ]
  X <- X[use, ]

  # inhom density
  wt <- rep(1, X$n)
  if (!is.null(sigma)) {
    np <- spatstat.geom::nearest.valid.pixel(X$x, X$y, den)
    w <- den$v[cbind(np$row, np$col)]
    wt <- 1 / w * mean(w)
    rm(np)
  }

  lam <- table(data$cellType) / spatstat.geom::area(X)
  lev <- levels(data$cellType)

  edge <- matrix(1, X$n, length(Rs) - 1)
  if (edgeCorrect) {
    for (k in seq_len(length(Rs) - 1)) edge[, k] <- borderEdge(X, Rs[k + 1])
  }

  # Pair counting, the L-function and the averaging over radii run in C++
  # (src/inhomL.cpp). The bin labels are passed as the R code compared them.
  L <- inhomLCpp(
    X$x, X$y, as.integer(data$cellType), length(lev), Rs,
    as.numeric(as.character(Rs[-1])), lev %in% from, lev %in% to,
    wt, as.numeric(lam), spatstat.geom::area(X), edge, edgeCorrect
  )
  dimnames(L) <- list(lev, lev)
  L <- as.data.frame(as.table(L), stringsAsFactors = FALSE)
  L <- L[!is.na(L$Freq), ]
  wt <- L$Freq
  names(wt) <- paste(L$Var1, L$Var2, sep = "__")

  m1 <- rep(from, times = length(to))
  m2 <- rep(to, each = length(from))
  labels <- paste(m1, m2, sep = "__")

  assoc <- rep(-sum(Rs), length(labels))
  names(assoc) <- labels
  if (!includeZeroCells) assoc[!(m1 %in% X$marks & m2 %in% X$marks)] <- NA
  assoc[names(wt)] <- wt
  names(assoc) <- labels

  assoc
}




#' @useDynLib spicyR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
borderEdge <- function(X, maxD) {
  W <- X$window
  bW <- spatstat.geom::union.owin(
    spatstat.geom::border(W, maxD, outside = FALSE),
    spatstat.geom::border(W, 2, outside = TRUE)
  )
  inB <- spatstat.geom::inside.owin(X$x, X$y, bW)
  e <- rep(1, X$n)
  if (any(inB)) {
    # Same as area(intersect.owin(discs(X[inB], maxD), W)): each disc is the
    # 128-gon spatstat.geom::disc() builds, clipped to the window in C++.
    rings <- if (W$type == "rectangle") {
      list(list(x = W$xrange[c(1, 2, 2, 1)], y = W$yrange[c(1, 1, 2, 2)]))
    } else {
      W$bdry
    }
    areas <- discWindowArea(X$x[inB], X$y[inB], maxD, 128L, rings)
    e[inB] <- areas / (pi * maxD^2)
  }

  e
}

#' @importFrom stats quantile
calcWeights <- function(rS, M1, M2, nCells, weightFactor, weightZThreshold = 0.1) {
  count1 <- as.vector(nCells[, M1])
  count2 <- as.vector(nCells[, M2])
  rS <- as.vector(rS)
  toWeight <- !is.na(rS)
  resSqToWeight <- rS[toWeight]
  count1ToWeight <- count1[toWeight]
  count2ToWeight <- count2[toWeight]
  if (length(count1ToWeight) <= 20) {
    warning("A cell type pair is seen less than 20 times, not using weights for this pair.") # nolint
    weightFunction <- rep(1, length(count1))
    return(weightFunction)
  }
  z1 <- mpdWeightFit(
    log10(resSqToWeight + 1), log10(count1ToWeight + 1), log10(count2ToWeight + 1),
    log10(as.numeric(count1) + 1), log10(as.numeric(count2) + 1)
  )
  if (is.null(z1)) {
    # Fallback when the C++ fit fails.
    .need("scam", "for the weight model (the built-in fit failed)")
    weightFunction <- scam::scam(
      log10(resSqToWeight + 1) ~ s(log10(count1ToWeight + 1), bs = "mpd") + s(log10(count2ToWeight + 1), bs = "mpd") # nolint
    ) # , optimizer = "nlm.fd")

    z1 <- suppressWarnings(stats::predict(weightFunction, data.frame(
      count1ToWeight = as.numeric(count1),
      count2ToWeight = as.numeric(count2)
    )))
  }

  # Floor for 1/z1 (caps the maximum weight): the 1st percentile of the weight-model
  # predictions above `weightZThreshold`. The default cutoff is calibrated to the
  # L-function's numeric scale; for a statistic whose residual variance is small in
  # absolute terms (e.g. getPairwiseProp()'s obs/exp ratio) every prediction can sit
  # below it, making the quantile NA and every weight NA -- pass weightZThreshold = 0
  # for such statistics.
  zFloor <- stats::quantile(z1[z1 > weightZThreshold], 0.01, na.rm = TRUE)
  if (is.na(zFloor)) {
    warning("Weight model predictions are all <= weightZThreshold; returning unweighted (weights = 1).") # nolint
    return(rep(1, length(count1)))
  }
  w <- 1 / pmax(z1, zFloor)
  w <- w / mean(w, na.rm = TRUE)
  w^weightFactor
}




# The basis of scam's monotone-decreasing P-spline smooth, s(x, bs = "mpd")
# with its defaults (k = 10, m = 2), as smooth.construct.mpd.smooth.spec()
# builds it: B-splines on evenly spaced knots, reparameterised so that
# exp(coefficients) give a decreasing curve, then centred.
mpdBasis <- function(x, q = 10, m = 2) {
  xk <- rep(0, q + m + 2)
  xk[(m + 2):(q + 1)] <- seq(min(x), max(x), length = q - m)
  for (i in 1:(m + 1)) xk[i] <- xk[m + 2] - (m + 2 - i) * (xk[m + 3] - xk[m + 2])
  for (i in (q + 2):(q + m + 2)) xk[i] <- xk[q + 1] + (i - q - 1) * (xk[m + 3] - xk[m + 2])
  Sig <- matrix(-1, q, q)
  Sig[upper.tri(Sig)] <- 0
  Sig[, 1] <- -Sig[, 1]
  X <- splines::splineDesign(xk, x, ord = m + 2)[, -1] %*% Sig[-1, -1]
  cmX <- colMeans(X)
  list(X = sweep(X, 2, cmX), knots = xk, cmX = cmX, Sigma = Sig, ord = m + 2)
}

# The prediction matrix of an mpdBasis() smooth at new x, as
# Predict.matrix.mpd.smooth() builds it: linear beyond the fitted range.
mpdPredict <- function(b, x) {
  ll <- b$knots[b$ord]
  ul <- b$knots[length(b$knots) - b$ord + 1]
  X <- matrix(0, length(x), ncol(b$Sigma))
  inside <- x >= ll & x <= ul
  if (any(inside)) X[inside, ] <- splines::splineDesign(b$knots, x[inside], b$ord)
  D <- splines::splineDesign(b$knots, c(ll, ll, ul, ul), b$ord, c(0, 1, 0, 1))
  below <- x < ll
  above <- x > ul
  if (any(below)) X[below, ] <- cbind(1, x[below] - ll) %*% D[1:2, ]
  if (any(above)) X[above, ] <- cbind(1, x[above] - ul) %*% D[3:4, ]
  sweep((X %*% b$Sigma)[, -1, drop = FALSE], 2, b$cmX)
}

# Minimise b'Qb - 2 c'b subject to b >= lower (primal active-set method). This
# is the problem pcls() solves for scam.fit()'s starting coefficients.
boxQP <- function(Q, c, lower) {
  b <- ifelse(is.finite(lower), pmax(lower, 0.1), 0)
  active <- rep(FALSE, length(c))
  for (it in seq_len(10 * length(c))) {
    f <- !active
    target <- b
    target[f] <- solve(Q[f, f, drop = FALSE], c[f] - Q[f, active, drop = FALSE] %*% b[active])
    blocked <- f & target < lower
    if (!any(blocked)) {
      b <- target
      lambda <- drop(Q %*% b - c)
      if (!any(active & lambda < 0)) return(b)
      active[which(active)[which.min(lambda[active])]] <- FALSE
    } else {
      t <- (b[blocked] - lower[blocked]) / (b[blocked] - target[blocked])
      j <- which(blocked)[which.min(t)]
      b <- b + min(t) * (target - b)
      b[j] <- lower[j]
      active[j] <- TRUE
    }
  }
  b
}

# calcWeights()'s model scam(y ~ s(x1, bs = "mpd") + s(x2, bs = "mpd")) with
# GCV smoothing parameters, fitted in C++ (src/scamMono.cpp) along the same
# path scam() takes: penalties scaled as mgcv's smoothCon() scales them, the
# pcls() start at sp = 0.05, then scam's BFGS search. Returns the predictions
# at (x1new, x2new), or NULL if the C++ fit fails.
#' @importFrom splines splineDesign
mpdWeightFit <- function(y, x1, x2, x1new, x2new) {
  b1 <- mpdBasis(x1)
  b2 <- mpdBasis(x2)
  X <- cbind(1, b1$X, b2$X)
  k <- ncol(b1$X)
  D <- crossprod(diff(diag(k)))
  S1 <- S2 <- matrix(0, ncol(X), ncol(X))
  S1[1 + seq_len(k), 1 + seq_len(k)] <- D / (norm(D) / norm(b1$X, type = "I")^2)
  S2[1 + k + seq_len(k), 1 + k + seq_len(k)] <- D / (norm(D) / norm(b2$X, type = "I")^2)
  iv <- c(FALSE, rep(TRUE, 2 * k))
  XtX <- crossprod(X)
  Xty <- drop(crossprod(X, y))
  start <- boxQP(XtX + 0.05 * (S1 + S2), Xty, ifelse(iv, 1e-12, -Inf))
  fit <- scamMonoFit(XtX, Xty, sum(y^2), length(y), list(S1, S2), iv, start,
                     c(1L, 1L + k), c(k, k))
  if (!fit$ok) return(NULL)
  z <- drop(cbind(1, mpdPredict(b1, x1new), mpdPredict(b2, x2new)) %*% fit$coef)
  names(z) <- seq_along(z) # as predict.scam() names them
  z
}


getWeightFunction <- function(
    pairwiseAssoc,
    nCells,
    m1,
    m2,
    BPPARAM,
    weights,
    weightsByPair,
    weightFactor,
    weightZThreshold = 0.1) {
  if (!weights) {
    weightFunction <- rep(1, nrow(pairwiseAssoc) * ncol(pairwiseAssoc))
    pair <- rep(colnames(pairwiseAssoc), each = nrow(pairwiseAssoc))
    weightFunction <- split(weightFunction, pair)
    return(weightFunction)
  }


  resSq <- apply(pairwiseAssoc, 2, function(x) {
    if (sum(!is.na(x)) > 1) {
      if (stats::sd(x, na.rm = TRUE) > 0) {
        return((x - mean(x, na.rm = TRUE))^2)
      } else {
        return(rep(NA, length(x)))
      }
    } else {
      return(rep(NA, length(x)))
    }
  })

  if (weightsByPair) {
    weightFunction <- .par_mapply(
      calcWeights,
      rS = as.list(as.data.frame(resSq)), M1 = m1, M2 = m2,
      MoreArgs = list(nCells = nCells, weightFactor = weightFactor, weightZThreshold = weightZThreshold),
      SIMPLIFY = FALSE, cores = .n_workers(BPPARAM)
    )
  } else {
    weightFunction <- calcWeights(m1, m2, rS = resSq, nCells, weightFactor, weightZThreshold)
    pair <- rep(colnames(pairwiseAssoc), each = nrow(pairwiseAssoc))
    weightFunction <- split(weightFunction, pair)
  }
  weightFunction[colnames(pairwiseAssoc)]
}





#' @importFrom methods is
prepCellSummary <- function(
    cells, spatialCoords, cellType, imageID, bind = FALSE) {
  getCellSummary(cells, bind = bind)
}




#' Perform a simple wilcoxon-rank-sum test or t-test on the columns of a data
#' frame
#'
#' @param df A data.frame or SingleCellExperiment, SpatialExperiment
#' @param condition The condition of interest
#' @param type The type of test, "wilcox", "ttest" or "survival".
#' @param feature
#'     Can be used to calculate the proportions of this feature for each image
#' @param imageID The imageID's if presenting a SingleCellExperiment
#'
#' @return Proportions
#'
#'
#' @examples
#'
#' # Test for an association with long-duration diabetes
#' # This is clearly ignoring the repeated measures...
#' data("diabetesData")
#' diabetesData <- spicyR:::.format_data(
#'   diabetesData, "imageID", "cellType", c("x", "y"), FALSE
#' )
#' props <- getProp(diabetesData)
#' condition <- spicyR:::getImagePheno(diabetesData)$stage
#' names(condition) <- spicyR:::getImagePheno(diabetesData)$imageID
#' condition <- condition[condition %in% c("Long-duration", "Onset")]
#' test <- colTest(props[names(condition), ], condition)
#' @export
#' @importFrom stats wilcox.test t.test
#' @importFrom S4Vectors as.data.frame
colTest <- function(
    df, 
    condition, 
    type = NULL, 
    feature = NULL, 
    imageID = "imageID") {
  if (.is_class(df, "SingleCellExperiment") || .is_class(df, "SpatialExperiment")) {
    if (is.null(feature)) stop("'feature' is still null")

    if (is.null(type) && length(condition) == 1) {
      type <- "ttest"
    } else if (is.null(type) && length(condition) == 2) {
      type <- "survival"
    } else if (is.null(type)) {
      stop("Invalid nuber of columns in condition. Must 1 or 2 (survival).")
    }

    x <- .col_data(df)
    x <- x[, c(imageID, condition)]
    x <- unique(x)
    condition <- x[[condition]]
    names(condition) <- x[[imageID]]

    df <- getProp(df, imageID = imageID, feature = feature)
    condition <- condition[rownames(df)]
  } else {
    stopifnot(
      "Number of oberservation does not match between df and outcome" = {
        dim(df)[1] == dim(condition)[1]
      }
    )

    if (is.null(type) && (is.atomic(condition) || dim(condition)[2] == 1)) {
      type <- "ttest"
    } else if (is.null(type) && dim(condition)[2] == 2) {
      type <- "survival"
    } else if (is.null(type)) {
      stop("Invalid nuber of columns in condition. Must 1 or 2 (survival).")
    }
  }

  if (type == "survival") {
    test <- colCoxTests(df, condition)
    test <- signif(test, 2)
    names(test)[names(test) == "p.value"] <- "pval"
  } else {
    test <- apply(df, 2, function(x) {
      if (type == "wilcox") {
        test <- stats::wilcox.test(x ~ condition)
      } else if (type == "ttest") {
        test <- stats::t.test(x ~ condition)
      }
      signif(c(test$estimate, tval = test$statistic, pval = test$p.value), 2)
    })
  }
  if (type != "survival") test <- as.data.frame(t(test))
  test$adjPval <- signif(stats::p.adjust(test$pval, "fdr"), 2)
  test$cluster <- rownames(test)
  test <- test[order(test$pval), ]
  test
}

# A Cox model (Efron ties, Wald test) of `outcome` (Surv or a time/event matrix) on each column of
# `measurements`, as ClassifyR::colCoxTests() gave it: a data.frame of coef, se.coef and p.value with
# one row per column. Uses the package's C++ Cox fit, with survival::coxph() if that fails.
colCoxTests <- function(measurements, outcome) {
  measurements <- as.matrix(measurements)
  outcome <- as.matrix(outcome)
  time <- as.numeric(outcome[, 1])
  event <- as.integer(outcome[, 2])
  out <- t(vapply(seq_len(ncol(measurements)), function(j) {
    x <- as.numeric(measurements[, j])
    ok <- stats::complete.cases(time, event, x)
    fit <- tryCatch(stats_cox_fit(time[ok], event[ok], matrix(x[ok], ncol = 1)), error = function(e) NULL)
    if (!is.null(fit) && isTRUE(fit$ok) && is.finite(fit$se[1])) {
      return(c(fit$beta[1], fit$se[1], fit$p[1]))
    }
    cf <- tryCatch(
      summary(survival::coxph(survival::Surv(time[ok], event[ok]) ~ x[ok]))$coefficients[1, c(1, 3, 5)],
      error = function(e) rep(NA_real_, 3)
    )
    as.numeric(cf)
  }, numeric(3)))
  output <- data.frame(coef = out[, 1], se.coef = out[, 2], p.value = out[, 3])
  rownames(output) <- colnames(measurements)
  output
}

#' Produces a dataframe showing L-function metric for each imageID entry.
#'
#' @param results
#'  Spicy test result obtained from spicy.
#' @param pairName
#'  A string specifying the pairwise interaction of interest. If NULL, all
#'  pairwise interactions are shown.
#'
#' @return A data.frame containing the colData related to the results.
#' @export
#'
#' @examples
#'
#' data(spicyTest)
#' df <- bind(spicyTest)
#'
#' @export
bind <- function(results,
                 pairName = NULL) {
  df <- data.frame(
    imageID = results$imageID,
    condition = results$condition
  )

  if (!is.null(results$subject)) {
    df <- cbind(df, subject = results$subject)
  }

  if (is.null(pairName)) {
    df <- cbind(df, do.call(cbind, results$pairwiseAssoc))
  } else {
    df[[pairName]] <- results$pairwiseAssoc[[pairName]]
  }

  return(df)
}
