#' Performs spatial tests on spatial cytometry data.
#'
#' @param cells A SummarizedExperiment or data frame that contains at least the  variables
#'   x and y, giving the location coordinates of each cell, and cellType.
#' @param condition A character specifying which column which contains the condition or `Surv` objects.
#' @param subject Vector of subject IDs corresponding to each image if cells is
#'   a data frame.
#' @param covariates Vector of covariate names that should be included in the
#'   mixed effects model as fixed effects.
#' @param imageID The name of the imageID column if using a SingleCellExperiment or SpatialExperiment.
#' @param cellType The name of the cellType column if using a SingleCellExperiment or SpatialExperiment.
#' @param spatialCoords The names of the spatialCoords column if using a SingleCellExperiment.
#' @param r A vector of the radii that the measures of association should be calculated over.
#' @param sigma A numeric variable used for scaling when fitting inhomogenous L-curves.
#' @param from vector of cell types which you would like to compare to the to vector.
#' @param to vector of cell types which you would like to compare to the from vector.
#' @param alternateResult A pairwise association statistic between each combination of celltypes in
#'   each image.
#' @param cores Number of cores to use for parallel processing or a BiocParallel MulticoreParam or SerialParam object.
#' @param minLambda Minimum value density for scaling when fitting inhomogeneous L-curves.
#' @param weights logical indicating whether to include weights based on cell counts.
#' @param weightsByPair logical indicating whether weights should be calculated for each cell type
#'   pair.
#' @param weightFactor numeric that controls the convexity of the weight function.
#' @param weightZThreshold numeric; the minimum weight-model prediction
#'   (\code{log10(resSq + 1)}) counted when choosing the \code{1/z} weight-cap
#'   floor. The default (\code{0.1}) suits the L-function's numeric scale. Pass
#'   \code{0} for a statistic whose residual variance is small in absolute terms
#'   (e.g. the observed/expected ratio from \code{\link{getPairwiseProp}}), where
#'   every prediction can otherwise sit below the default and collapse the
#'   weights to \code{NA}.
#' @param window 	Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length A tuning parameter for controlling the level of concavity when estimating concave windows.
#' @param edgeCorrect A logical indicating whether to perform edge correction.
#' @param includeZeroCells 	A logical indicating whether to include cells with zero counts in the pairwise association calculation.
#' @param verbose logical indicating whether to output messages.
#' @param BPPARAM \{DEPRECATED\} A BiocParallel MulticoreParam or SerialParam object. 
#' @param imageIDCol \{DEPRECATED\} The name of the imageID column if using a SingleCellExperiment or SpatialExperiment.
#' @param cellTypeCol \{DEPRECATED\} The name of the cellType column if using a SingleCellExperiment or SpatialExperiment.
#' @param spatialCoordCols \{DEPRECATED\} The names of the spatialCoords column if using a SingleCellExperiment.
#' @param nCores \{DEPRECATED\} Number of cores to use for parallel processing or a BiocParallel MulticoreParam or SerialParam object.
#' @param Rs \{DEPRECATED\} A vector of the radii that the measures of association should be calculated over.
#' @param ... Other options
#' @return Data frame of p-values.
#' @export
#'
#' @examples
#' data("diabetesData")
#'
#' # Test with random effect for patient on a pairwise combination of cell
#' # types.
#' spicy(diabetesData,
#'   condition = "stage", subject = "case",
#'   from = "Tc", to = "Th"
#' )
#'
#' # Test all pairwise combinations of cell types without random effect of
#' # patient.
#' \dontrun{
#' spicyTest <- spicy(diabetesData, condition = "stage", subject = "case")
#' }
#'
#' # Test all pairwise combination of cell types with random effect of patient.
#' \dontrun{
#' spicy(diabetesData, condition = "condition", subject = "subject")
#' }
#'
#' @aliases
#' spicy
#' spicy,spicy-method
#' @importFrom scam scam
#' @importFrom rlang .data
#' @importFrom tibble column_to_rownames
#' @importFrom lifecycle deprecate_soft
spicy <- function(cells,
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
  argumentChecks("spicy", user_vals)
  
  if (is.null(BPPARAM)) {
    # Built on first use: making a MulticoreParam takes ~0.2 s, and only
    # per-pair weights and survival models use it.
    delayedAssign("BPPARAM", {
      if (cores > 1 && .Platform$OS.type != "windows") {
        BiocParallel::MulticoreParam(workers = cores)
      } else if (cores > 1) {
        BiocParallel::SnowParam(workers = cores)
      } else {
        BiocParallel::SerialParam()
      }
    })
  }
  
  if (is(cells, "SummarizedExperiment") || is(cells, "data.frame")) {
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
    m1 <- from
    m2 <- to
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
      conditionVector <- relevel(conditionVector, ref = levels(conditionVector)[1])
      
      if (!wasFactor || TRUE) {  
        cli::cli_inform(
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
      
      comparisons <- data.frame(labels) |>
        tidyr::separate(
          col = labels,
          into = c("from", "to", "parent"),
          sep = "__"
        ) |>
        dplyr::mutate(labels = paste(.data$from, .data$to, .data$parent, sep = "__"))
      
      m1 <- comparisons$from
      m2 <- comparisons$to
      
      spicyResult$isKontextual = TRUE
    
    }
  }
  
  
  weightFunction <- getWeightFunction(
    pairwiseAssoc, nCells, m1, m2, BPPARAM, weights, weightsByPair, weightFactor,
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
      BPPARAM = BPPARAM
    )
    
    spicyResult$survivalOutcome = conditionVector
    spicyResult = append(spicyResult, list("survivalResults" = survivalResult))
    
  }


  ## Linear model
  if (!inherits(conditionVector, "Surv") && is.null(subject) && !is.null(condition)) {
    if (verbose) {
      cli::cli_inform("Testing for spatial differences across conditions")
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

    lmResult <- cleanLM(linearModels, BPPARAM = BPPARAM)
    spicyResult = append(spicyResult, lmResult)
  }


  ## Mixed effects model
  if (!inherits(conditionVector, "Surv") && (!is.null(subject)) && !is.null(condition)) {
    if (verbose) {
      cli::cli_inform(
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


    melmResult <- cleanMEM(mixed.lmer, BPPARAM = BPPARAM)
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




#' @importFrom dplyr bind_rows
cleanLM <- function(linearModels, BPPARAM) {
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
    dplyr::bind_rows(x)
  }))

  df <- lapply(df, function(x) {
    rownames(x) <- names(linearModels)
    x
  })
  df
}


#' @importFrom dplyr bind_rows
cleanMEM <- function(mixed.lmer, BPPARAM) {
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
    dplyr::bind_rows(x)
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
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel MulticoreParam
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

  
  if (is(cells, "SummarizedExperiment")) {
    cells <- .format_data(
      cells, imageID, cellType, spatialCoords, FALSE
    )
  }
    
  # Square and convex windows without density weighting: one threaded C++
  # call for all images, with `cores` threads.
  lev <- levels(cells$cellType)
  if (is.null(sigma) && window %in% c("square", "convex") && !is.null(lev) &&
      all(c(from, to) %in% lev)) {
    nThreads <- if (!is.null(BPPARAM)) {
      BiocParallel::bpnworkers(BPPARAM)
    } else if (is.numeric(cores)) {
      cores
    } else {
      BiocParallel::bpnworkers(cores)
    }
    return(getPairwiseThreaded(
      cells, Rs, from, to, window, edgeCorrect, includeZeroCells, nThreads
    ))
  }

  if (is.null(BPPARAM)) {
    if (cores > 1 && .Platform$OS.type != "windows") {
      BPPARAM = BiocParallel::MulticoreParam(workers = cores)
    } else if (cores > 1) {
      BPPARAM = BiocParallel::SnowParam(workers = cores)
    } else {
      BPPARAM = BiocParallel::SerialParam()
    } 
  }  

  cells2 <- getCellSummary(cells, bind = FALSE)


  if (is.null(from)) from <- levels(cells2$cellType)
  if (is.null(to)) to <- levels(cells2$cellType)

  pairwiseVals <- BiocParallel::bplapply(cells2,
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
    BPPARAM = BPPARAM
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
#' @importFrom SummarizedExperiment colData
getProp <- function(cells, feature = "cellType", imageID = "imageID") {
  if (is.data.frame(cells)) {
    df <- cells[, c(imageID, feature)]
  }

  if (is(cells, "SingleCellExperiment") || is(cells, "SpatialExperiment")) {
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
  pval <- as.data.frame(df$p.value)
  cond <- colnames(pval)[grep("condition", colnames(pval))]
  message(df$test)
  message("Number of cell type pairs: ", nrow(pval), "\n")
  message("Number of differentially localised cell type pairs: \n")
  if (nrow(pval) == 1) {
    print(sum(pval[cond] < 0.05, na.rm = TRUE))
  }
  if (nrow(pval) > 1) {
    print(colSums(apply(pval[cond], 2, p.adjust, "fdr") < 0.05, na.rm = TRUE))
  }
}
setMethod(
  "show", methods::signature(object = "SpicyResults"), function(object) {
    .show_SpicyResults(object)
  }
)

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
#' @importFrom dplyr bind_rows mutate rename select arrange across where 
#' @importFrom coxme coxme
#' @importFrom BiocParallel bplapply SerialParam
spatialSurv <- function(measurementMat,
                        condition,
                        pheno,
                        covariates = NULL,
                        subject = NULL,
                        weights = NULL,
                        remove = NULL,
                        BPPARAM = NULL) {
  
  result <- bplapply(colnames(measurementMat), function(test) {
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
    spatialData = spatialData %>%
      mutate(across(where(is.factor), ~ droplevels(.)))
    
    # Mixed effects survival
    fit <- coxme::coxme(formula, data = spatialData, weights = spatialData$weights, subset = ind)
    
    # from print.coxme() function
    beta <- unname(fit$coefficients)
    nvar <- length(beta)
    nfrail <- nrow(fit$var) - nvar
    
    se <- sqrt(diag(as.matrix(fit$var))[nfrail + 1:nvar])
    p_val <- 1 - pchisq((beta / se) ^ 2, 1)
    
    # Need to select first indexes of each value incase of multiple betas
    result <- c("coef" = beta[1], "se(coef)" = se[1], "Pr(>|z|)" = p_val[1])
  }
  
  return(result)
  }, BPPARAM = BPPARAM)

result <- result |>
  dplyr::bind_rows() |>
  dplyr::mutate(test = colnames(measurementMat)) |>
  dplyr::rename("se.coef" = "se(coef)", "p.value" = "Pr(>|z|)") |>
  dplyr::select(test, coef, se.coef, p.value) |>
  dplyr::arrange(p.value)

return(result)
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
    data <- data.frame(data)
    ow <-
      spatstat.geom::owin(xrange = range(data$x), yrange = range(data$y))

    if (window == "convex") {
      p <- spatstat.geom::ppp(data$x, data$y, ow)
      ow <- spatstat.geom::convexhull(p)
    }
    if (window == "concave") {
      cli::cli_inform("Concave windows are temperamental. Try choosing values of window.length > and < 1 if you have problems.") # nolint
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




#' @importFrom spatstat.explore density.ppp
#' @importFrom spatstat.geom nearest.valid.pixel area ppp
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




#' @importFrom spatstat.geom union.owin border inside.owin
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

#' @importFrom scam scam
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


#' @importFrom BiocParallel bpmapply
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
    weightFunction <- BiocParallel::bpmapply(
      calcWeights,
      rS = as.list(as.data.frame(resSq)), M1 = m1, M2 = m2, BPPARAM = BPPARAM,
      MoreArgs = list(nCells = nCells, weightFactor, weightZThreshold), SIMPLIFY = FALSE
    )
  } else {
    weightFunction <- calcWeights(m1, m2, rS = resSq, nCells, weightFactor, weightZThreshold)
    pair <- rep(colnames(pairwiseAssoc), each = nrow(pairwiseAssoc))
    weightFunction <- split(weightFunction, pair)
  }
  weightFunction[colnames(pairwiseAssoc)]
}





#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
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
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom stats wilcox.test t.test
#' @importFrom S4Vectors as.data.frame
#' @importFrom ClassifyR colCoxTests
colTest <- function(
    df, 
    condition, 
    type = NULL, 
    feature = NULL, 
    imageID = "imageID") {
  if (is(df, "SingleCellExperiment") || is(df, "SpatialExperiment")) {
    if (is.null(feature)) stop("'feature' is still null")

    if (is.null(type) && length(condition) == 1) {
      type <- "ttest"
    } else if (is.null(type) && length(condition) == 2) {
      type <- "survival"
    } else if (is.null(type)) {
      stop("Invalid nuber of columns in condition. Must 1 or 2 (survival).")
    }

    x <- df@colData
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
    test <- ClassifyR::colCoxTests(df, condition)
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
    df <- cbind(df, results$pairwiseAssoc[[pairName]])
  }

  return(df)
}
