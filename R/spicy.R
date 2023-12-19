#' Performs spatial tests on spatial cytometry data.
#'
#' @param cells
#'   A SegmentedCells or data frame that contains at least the  variables x and
#'   y, giving the location coordinates of each cell, and cellType.
#' @param condition
#'   Vector of conditions to be tested corresponding to each image if cells is
#'   a data frame.
#' @param subject Vector of subject IDs corresponding to each image if cells is
#'   a data frame.
#' @param covariates Vector of covariate names that should be included in the
#'   mixed effects model as fixed effects.
#' @param from
#'   vector of cell types which you would like to compare to the to vector
#' @param to
#'   vector of cell types which you would like to compare to the from vector
#' @param alternateResult
#'   An pairwise association statistic between each combination of celltypes in each image.
#' @param verbose logical indicating whether to output messages.
#' @param weights
#'   logical indicating whether to include weights based on cell counts.
#' @param weightsByPair
#'   logical indicating whether weights should be calculated for each cell type
#'   pair.
#' @param weightFactor
#'   numeric that controls the convexity of the weight function.
#' @param window
#'   Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length
#'   A tuning parameter for controlling the level of concavity when estimating
#'   concave windows.
#' @param BPPARAM A BiocParallelParam object.
#' @param sigma
#'   A numeric variable used for scaling when fitting inhomogeneous L-curves.
#' @param minLambda
#'   Minimum value for density for scaling when fitting inhomogeneous L-curves.
#' @param Rs
#'   A vector of radii that the measures of association should be calculated.
#' @param edgeCorrect A logical indicating whether to perform edge correction.
#' @param includeZeroCells
#'   A logical indicating whether to include cells with zero counts in the
#'   pairwise association calculation.
#' @param imageID The image ID if using SingleCellExperiment.
#' @param cellType The cell type if using SingleCellExperiment.
#' @param spatialCoords The spatial coordinates if using a SingleCellExperiment.
#' @param ... Other options.
#'
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
#' # spicyTest <- spicy(diabetesData, condition = "stage", subject = "case")
#'
#' # Test all pairwise combination of cell types with random effect of patient.
#' # spicy(diabetesData, condition = "condition", subject = "subject")
#'
#' @aliases
#' spicy
#' spicy,spicy-method
#' @importFrom mgcv gam ti
#' @importFrom BiocParallel SerialParam
#' @importFrom scam scam
#' @importFrom rlang .data
#' @importFrom tibble column_to_rownames
spicy <- function(cells,
                  condition = NULL,
                  subject = NULL,
                  covariates = NULL,
                  from = NULL,
                  to = NULL,
                  alternateResult = NULL,
                  verbose = TRUE,
                  weights = TRUE,
                  weightsByPair = FALSE,
                  weightFactor = 1,
                  window = "convex",
                  window.length = NULL,
                  BPPARAM = BiocParallel::SerialParam(),
                  sigma = NULL,
                  Rs = NULL,
                  minLambda = 0.05,
                  edgeCorrect = TRUE,
                  includeZeroCells = FALSE,
                  imageID = "imageID",
                  cellType = "cellType",
                  spatialCoords = c("x", "y"),
                  ...) {
  if (is(cells, "SingleCellExperiment") || is(cells, "SpatialExperiment")) {
    cells <- extractSpicyInfo(cells,
      imageID = imageID,
      cellType = cellType,
      spatialCoords = spatialCoords,
      condition = condition,
      subject = subject,
      covariates = covariates
    )
  } # TODO: Decrecate segementedCells

  if (!is(cells, "SegmentedCells")) {
    stop("cells needs to be a SegmentedCells object") # TODO: Make this error informative.
  }


  if (is.null(from) || is.null(to)) {
    if (is.null(from)) {
      from <- as.character(unique(cellType(cells)))
    }
    if (is.null(to)) {
      to <- as.character(unique(cellType(cells)))
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


  if (any((!to %in% cellType(cells)) | (!from %in% cellType(cells)))) {
    stop("to and from need to be cell type in your SegmentedCells")
  }

  nCells <- table(imageID(cells), cellType(cells))
  
  conditionVector <- as.data.frame(imagePheno(cells))[condition][,1]
  
  if(!is.factor(conditionVector)) {
    conditionVector <- as.factor(conditionVector)
    warning(paste0("Coercing condition into factor. Using ", condition, " = ", levels(conditionVector)[1], " as base comparison group. If this is not the desired base group, please convert cells$", condition ," into a factor and change the order of levels(cells$", condition,") so that the base group is at index 1."))
  }
  
  ## Find pairwise associations

  if (is.null(alternateResult)) {
    pairwiseAssoc <- getPairwise(cells,
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
    pairwiseAssoc <- as.data.frame(pairwiseAssoc)
    pairwiseAssoc <- pairwiseAssoc[labels]
  }


  if (!is.null(alternateResult) && !isKonditional(alternateResult)) {
    pairwiseAssoc <- alternateResult
    
    weights <- FALSE
    message("Cell count weighting set to FALSE for alternate results")
  }

  comparisons <- data.frame(from = m1, to = m2, labels = labels)

  if (!is.null(alternateResult) && isKonditional(alternateResult)) {
    pairwiseAssoc <- alternateResult |>
      dplyr::select(.data$imageID, .data$test, .data$konditional) |>
      tidyr::pivot_wider(
        names_from = .data$test, values_from = .data$konditional
      ) |>
      tibble::column_to_rownames("imageID")

    labels <- names(pairwiseAssoc)

    comparisons <- data.frame(labels) |>
      tidyr::separate(
        col = labels,
        into = c("fromName", "to", "parent"),
        sep = "__"
      ) |>
      dplyr::mutate(
        labels = paste(.data$fromName, .data$to, .data$parent, sep = "__")
      ) |>
      dplyr::mutate(from = paste(.data$fromName, .data$parent, sep = "__")) |>
      dplyr::select(-.data$parent)

    m1 <- comparisons$fromName
    m2 <- comparisons$to

    weights <- FALSE
    message("Cell count weighting set to FALSE for Konditional results")
  }

  weightFunction <- getWeightFunction(
    pairwiseAssoc, nCells, m1, m2, BPPARAM, weights, weightsByPair, weightFactor
  )

  pairwiseAssoc <- as.list(data.frame(pairwiseAssoc))
  names(pairwiseAssoc) <- labels


  ## Linear model
  if (is.null(subject) && !is.null(condition)) {
    if (verbose) {
      message("Testing for spatial differences across conditions")
    }

    MoreArgs2 <-
      list(
        cells = cells,
        condition = condition,
        covariates = covariates,
        cellCounts = table(imageID(cells), cellType(cells)),
        pheno = as.data.frame(imagePheno(cells))
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

    df <- cleanLM(linearModels, BPPARAM = BPPARAM)
  }


  ## Mixed effects model
  if ((!is.null(subject)) && !is.null(condition)) {
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
        cellCounts = table(imageID(cells), cellType(cells)),
        pheno = as.data.frame(imagePheno(cells))
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


    mixed.lmer <- lapply(mixed.lmer, function(x) {
      if (is(x, "lmerModLmerTest")) {
        if (x@devcomp$cmp["REML"] == -Inf) {
          return(NA)
        }
      }
      x
    })

    df <- cleanMEM(mixed.lmer, BPPARAM = BPPARAM)
  }
  
  df$condition <- conditionVector
  
  if(!is.null(subject)) {
    df$subject <- as.data.frame(imagePheno(cells))[subject][,1]    
  }

  df$pairwiseAssoc <- pairwiseAssoc
  df$comparisons <- comparisons

  df$weights <- weightFunction
  df$nCells <- nCells
  
  df$imageIDs <- as.data.frame(imagePheno(cells))["imageID"][,1]
  df$alternateResult <- ifelse(is.null(alternateResult), FALSE, TRUE)
  
  df <- methods::new("SpicyResults", df)
  df
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
    if (is(lmer, "lmerMod")) {
      coef <- as.data.frame(t(summary(lmer)$coef))
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
#' @param cells A SegmentedCells or data frame that contains at least the
#'     variables x and y, giving the location coordinates of each cell, and
#'     cellType.
#' @param from The 'from' cellType for generating the L curve.
#' @param to The 'to' cellType for generating the L curve.
#' @param window
#'     Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length
#'     A tuning parameter for controlling the level of concavity
#'     when estimating concave windows.
#' @param Rs
#'     A vector of the radii that the measures of association should be
#'     calculated.
#' @param sigma
#'     A numeric variable used for scaling when fitting inhomogeneous L-curves.
#' @param minLambda
#'     Minimum value for density for scaling when fitting inhomogeneous
#'     L-curves.
#' @param edgeCorrect A logical indicating whether to perform edge correction.
#' @param includeZeroCells A logical indicating whether to include cells with
#' zero counts in the pairwise association calculation.
#' @param BPPARAM A BiocParallelParam object.
#' @param imageID
#'     The imageID if using a SingleCellExperiment or SpatialExperiment.
#' @param cellType
#'     The cellType if using a SingleCellExperiment or SpatialExperiment.
#' @param spatialCoords
#'     The spatialCoords if using a SingleCellExperiment or SpatialExperiment.
#' @return Statistic from pairwise L curve of a single image.
#'
#'
#' @examples
#' data("diabetesData")
#' pairAssoc <- getPairwise(diabetesData[1, ])
#' @export
#' @importFrom BiocParallel bplapply
getPairwise <- function(
    cells,
    from = NULL,
    to = NULL,
    window = "convex",
    window.length = NULL,
    Rs = c(20, 50, 100),
    sigma = NULL,
    minLambda = 0.05,
    edgeCorrect = TRUE,
    includeZeroCells = TRUE,
    BPPARAM = BiocParallel::SerialParam(),
    imageID = "imageID",
    cellType = "cellType",
    spatialCoords = c("x", "y")) {
  cells2 <- prepCellSummary(
    cells, spatialCoords, cellType, imageID,
    bind = FALSE
  )


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




#' Get proportions from a SegmentedCells, SingleCellExperiment,
#' SpatialExperiment or data.frame.
#'
#' @param cells
#'     SegmentedCells, SingleCellExperiment, SpatialExperiment or data.frame
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

  if (is(cells, "SegmentedCells")) {
    cellSummary <- cellSummary(cells, bind = TRUE)
    df <- as.data.frame(cellSummary[, c(imageID, feature)])
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

    names(spatialData)[names(spatialData) == "covariates"] <- covariates

    formula <- "spatAssoc ~ condition + (1|subject)"

    if (!is.null(covariates)) {
      formula <-
        paste("spatAssoc ~ condition + (1|subject)",
          paste(covariates, collapse = "+"),
          sep = "+"
        )
    }

    spatialData$weights <- weightFunction

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


    spatialData <-
      data.frame(
        "spatAssoc" = spatAssoc,
        "condition" = pheno[, condition],
        "covariates" = pheno[, covariates]
      )

    names(spatialData)[names(spatialData) == "covariates"] <- covariates

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




#' @importFrom spatstat.explore density.ppp
#' @importFrom spatstat.geom closepairs nearest.valid.pixel area ppp
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join
#' @importFrom rlang .data
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
      marks = data$cellType
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
  data <- data[use, ]
  X <- X[use, ]


  p <- spatstat.geom::closepairs(X, max(Rs), what = "ijd", distinct = FALSE)

  p$j <- data$cellID[p$j]
  p$i <- data$cellID[p$i]

  cT <- data$cellType
  names(cT) <- data$cellID

  p$d <- cut(p$d, Rs, labels = Rs[-1], include.lowest = TRUE)

  # inhom density
  p$wt <- rep(1, length(p$d))
  if (!is.null(sigma)) {
    np <- spatstat.geom::nearest.valid.pixel(X$x, X$y, den)
    w <- den$v[cbind(np$row, np$col)]
    names(w) <- data$cellID
    p$wt <- 1 / w[p$j] * mean(w)
    rm(np)
  }


  lam <- table(data$cellType) / spatstat.geom::area(X)

  p$cellTypeJ <- cT[p$j]
  p$cellTypeI <- cT[p$i]
  p$i <- factor(p$i, levels = data$cellID)

  if (edgeCorrect) {
    edge <- sapply(Rs[-1], function(x) borderEdge(X, x), simplify = FALSE)
    edge <- do.call("cbind", edge)
    edge <- as.data.frame(edge)
    colnames(edge) <- Rs[-1]
    edge$i <- data$cellID
    edge <- tidyr::pivot_longer(edge, -.data$i, names_to = "d")
    p <- dplyr::left_join(as.data.frame(p), edge, c("i", "d"))
  } else {
    p <- as.data.frame(p)
    p$value <- 1
  }


  p$d <- factor(p$d, levels = Rs[-1])

  p <- p[p$i != p$j, ]

  use <- p$cellTypeI %in% from & p$cellTypeJ %in% to
  p <- p[use, ]

  r <- inhomL(p, lam, X, Rs)

  wt <- r$wt
  names(wt) <- paste(r$cellTypeI, r$cellTypeJ, sep = "__")

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


#' @importFrom data.table as.data.table setkey CJ .SD ":="
#' @importFrom spatstat.geom area
inhomL <-
  function(p, lam, X, Rs) {
    r <- data.table::as.data.table(p)
    r$wt <- r$wt / r$value / as.numeric(lam[r$cellTypeJ]) / as.numeric(lam[r$cellTypeI]) / spatstat.geom::area(X) # nolint
    r <- r[, j := NULL] # nolint
    r <- r[, value := NULL] # nolint
    r <- r[, i := NULL] # nolint
    data.table::setkey(r, d, cellTypeI, cellTypeJ) # nolint
    r <- r[data.table::CJ(d, cellTypeI, cellTypeJ, unique = TRUE)][, lapply(.SD, sum), by = .(d, cellTypeI, cellTypeJ)][is.na(wt), wt := 0] # nolint
    r <- r[, wt := cumsum(wt), by = list(cellTypeI, cellTypeJ)] # nolint
    r <- r[, list(wt = sum(sqrt(wt / pi))), by = .(cellTypeI, cellTypeJ)] # nolint
    r$wt <- r$wt - sum(Rs) # nolint

    r <- as.data.frame(r)

    r
  }




#' @importFrom spatstat.geom
#'     union.owin border inside.owin solapply intersect.owin area discs
borderEdge <- function(X, maxD) {
  W <- X$window
  bW <- spatstat.geom::union.owin(
    spatstat.geom::border(W, maxD, outside = FALSE),
    spatstat.geom::border(W, 2, outside = TRUE)
  )
  inB <- spatstat.geom::inside.owin(X$x, X$y, bW)
  e <- rep(1, X$n)
  if (any(inB)) {
    circs <- spatstat.geom::discs(X[inB], maxD, separate = TRUE)
    circs <- spatstat.geom::solapply(
      circs, spatstat.geom::intersect.owin, X$window
    )
    areas <- unlist(lapply(circs, spatstat.geom::area)) / (pi * maxD^2)
    e[inB] <- areas
  }

  e
}

#' @importFrom scam scam
#' @importFrom stats quantile
calcWeights <- function(rS, M1, M2, nCells, weightFactor) {
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
  weightFunction <- scam::scam(
    log10(resSqToWeight + 1) ~ s(log10(count1ToWeight + 1), bs = "mpd") + s(log10(count2ToWeight + 1), bs = "mpd") # nolint
  ) # , optimizer = "nlm.fd")

  z1 <- suppressWarnings(stats::predict(weightFunction, data.frame(
    count1ToWeight = as.numeric(count1),
    count2ToWeight = as.numeric(count2)
  )))
  w <- 1 / pmax(z1, stats::quantile(z1[z1 > 0.1], 0.01, na.rm = TRUE))
  w <- w / mean(w, na.rm = TRUE)
  w^weightFactor
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
    weightFactor) {
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
      MoreArgs = list(nCells = nCells, weightFactor), SIMPLIFY = FALSE
    )
  } else {
    weightFunction <- calcWeights(m1, m2, rS = resSq, nCells, weightFactor)
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
  if (is.data.frame(cells)) {
    cells <- SegmentedCells(cells,
      spatialCoords = spatialCoords,
      cellTypeString = cellType,
      imageIDString = imageID
    )

    cellSummary <- cellSummary(cells, bind = bind)
  }

  if (is(cells, "SingleCellExperiment")) {
    cells <- colData(cells)
    cells <- SegmentedCells(cells,
      spatialCoords = spatialCoords,
      cellTypeString = cellType,
      imageIDString = imageID
    )

    cellSummary <- cellSummary(cells, bind = bind)
  }

  if (is(cells, "SpatialExperiment")) {
    cells <- cbind(colData(cells), spatialCoords(cells))
    cells <- SegmentedCells(cells,
      spatialCoords = spatialCoords,
      cellTypeString = cellType,
      imageIDString = imageID
    )

    cellSummary <- cellSummary(cells, bind = bind)
  }

  if (is(cells, "SegmentedCells")) {
    cellSummary <- cellSummary(cells, bind = bind)
  }

  cellSummary
}


extractSpicyInfo <- function(cells,
                             imageID = imageID,
                             cellType = cellType,
                             spatialCoords = spatialCoords,
                             condition = condition,
                             subject = subject,
                             covariates = covariates) {
  extra <- c(condition, subject, covariates)

  if (is(cells, "SpatialExperiment")) {
    cells <- cbind(colData(cells), spatialCoords(cells))
    colnames(cells)[colnames(cells) %in% extra] <- paste0(
      "phenotype_", colnames(cells)[colnames(cells) %in% extra]
    )
    cells <- SegmentedCells(cells,
      spatialCoords = spatialCoords,
      cellTypeString = cellType,
      imageIDString = imageID
    )
  }

  if (is(cells, "SingleCellExperiment")) {
    cells <- colData(cells)
    colnames(cells)[colnames(cells) %in% extra] <- paste0(
      "phenotype_", colnames(cells)[colnames(cells) %in% extra]
    )
    cells <- SegmentedCells(cells,
      spatialCoords = spatialCoords,
      cellTypeString = cellType,
      imageIDString = imageID
    )
  }

  cells
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
#' props <- getProp(diabetesData)
#' condition <- imagePheno(diabetesData)$stage
#' names(condition) <- imagePheno(diabetesData)$imageID
#' condition <- condition[condition %in% c("Long-duration", "Onset")]
#' test <- colTest(props[names(condition), ], condition)
#' @export
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom stats wilcox.test t.test
#' @importFrom S4Vectors as.data.frame
#' @importFrom ClassifyR colCoxTests
colTest <- function(
    df, condition, type = NULL, feature = NULL, imageID = "imageID") {
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

  df <- data.frame(imageID = results$imageID,
                   condition = results$condition)
  
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

