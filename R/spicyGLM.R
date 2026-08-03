#' `Calculates pairwise spatial associations between cell types across images
#' and fit generalized linear model (GLM) models to test for condition effects.
#' 
#' @param cells A \code{SpatialExperiment}, \code{SingleCellExperiment}, or \code{data.frame}. 
#' The dataframe must have rows as markers and columns as cells.
#' containing single-cell or spatial data with cell metadata and coordinates.
#' @param condition  A character specifying which column in \code{cells} which contains the condition or grouping variable. 
#' @param subject A character specifying which column in \code{cells} which contains the patient/donor ID.
#' @param imageID A character specifying which column in \code{cells} which contains image/sample ID.
#' @param cellType A character specifying which column in \code{cells} which contains the cell types.
#' @param spatialCoords A character vector of length 2 specifying the columns for x and y coordinates if using a \code{SingleCellExperiment} object.
#' @param r Radius around each reference cell to consider for counting neighboring cells.
#' @param from Character vector of reference cell types. If NULL, all cell types are used.
#' @param to Character vector of target cell types. If NULL, all cell types are used.
#' @param window Defines the spatial window for each image. Options: "convex", "concave", or "rectangle".
#' @param cores Number of cores to use for parallel computation. 
#' @param cr2Method Character specifying how to compute the CR2-adjusted variance and p-value. 
#'   \code{"fast"} (default) closely matches \code{clubSandwich} at substantially lower cost.  
#'   \code{"clubSandwich"} uses the \code{clubSandwich} package directly.
#' @param fastMethod Eigendecomposition routine used internally when
#'   \code{cr2Method = "fast"}. \code{"direct"} (default) is faster in practice; \code{"dpr1"} is a validated reference alternative.
#' @param estimator Character specifying how the condition contrast is fit.
#'   \code{"firth"} (default) uses Firth's bias-reduced Poisson estimator, which
#'   remains finite even when a cell-type pair has structural zeros in one
#'   condition. \code{"mle"} uses ordinary Poisson MLE, which is skipped
#'   (with a warning) for such pairs, since the log rate ratio is non-estimable.
#' @param firthBackend Character specifying how Firth's estimator is computed
#'   when \code{estimator = "firth"}. \code{"closed_form"} (default) uses a fast,
#'   exact closed-form solution valid for spicyGLM's design; \code{"brglm2"}
#'   uses \code{brglm2::brglmFit} as a general fallback.
#' @param computeLeverage Logical; if \code{TRUE}, also computes and attaches
#'   a per-patient CR2 leverage diagnostic table (\code{$leverage} on the
#'   returned object). Adds some compute cost per pair; \code{FALSE} by default.
#'   
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{condition}{Factor vector of the condition used in the GLM models.}
#'   \item{subject}{Factor vector of subjects/donors, if provided.}
#'   \item{comparisons}{Data frame with the reference and target cell types for each pair
#'     and a combined label (from__to).}
#'   \item{nCells}{Table of cell counts per image and cell type.}
#'   \item{GLMresults}{Data frame of Poisson GLM results for each cell type pair using CR2.}
#' }
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' kerenSPE = SpatialDatasets::spe_Keren_2018()
#' spicyResult = spicyGLM(
#'   cells = kerenSPE,
#'   condition = "tumour_type",
#'   subject = "DONOR_NO",
#'   imageID = "imageID",
#'   from = "CD8_T_cell",
#'   to = "Tumour",
#'   spatialCoords = c("x", "y"),
#'   r = 50,
#'   window = "convex",
#'   cores = 1)
#' }
#' 
#' @importFrom cli cli_inform
#' 

# Internal constructor for a consistent SpicyResults object (GLM version)
.new_spicyGLM_result <- function(cells,
                                 condition,
                                 subject = NULL,
                                 imageID = "imageID",
                                 cellType = "cellType",
                                 GLMresults = NULL,
                                 messages = character(0)) {
  
  # --- condition vector (per image) ---
  conditionVector <- as.data.frame(getImagePheno(cells, imageID = imageID))[[condition]]
  wasFactor <- is.factor(conditionVector)
  if (!wasFactor) conditionVector <- as.factor(conditionVector)
  conditionVector <- droplevels(conditionVector)
  
  # NOTE: releveling is already done in spicyGLM() with cli_inform;
  # but we also do it here to ensure consistent base behavior even for early returns.
  if (nlevels(conditionVector) >= 1) {
    conditionVector <- relevel(conditionVector, ref = levels(conditionVector)[1])
  }
  
  # --- subject vector (per image), optional ---
  subjectVector <- NULL
  if (!is.null(subject)) {
    subjectVector <- as.data.frame(getImagePheno(cells, imageID = imageID))[[subject]]
    if (!is.factor(subjectVector)) subjectVector <- as.factor(subjectVector)
    subjectVector <- droplevels(subjectVector)
  }
  
  # --- nCells table (image x cellType) ---
  nCellsTab <- table(getImageID(cells), getCellType(cells))
  
  # --- enforce fixed schema for GLMresults ---
  required_cols <- c("from", "to", "conditionRef", "conditionComp", "coef_ref", "coef_comp",
                     "logRateRatio", "rateRatio", "p.value", "estimator",
                     "mle_would_skip", "mle_skip_reason")
  char_cols <- c("from", "to", "conditionRef", "conditionComp", "estimator", "mle_skip_reason")
  if (is.null(GLMresults)) {
    GLMresults <- as.data.frame(
      setNames(
        lapply(required_cols, function(cn) if (cn %in% char_cols) character(0) else numeric(0)),
        required_cols
      ),
      stringsAsFactors = FALSE
    )
  } else {
    # Ensure it's a data.frame
    GLMresults <- as.data.frame(GLMresults)
    
    # If missing required columns, add them as NA (keeps downstream code stable)
    missing_cols <- setdiff(required_cols, colnames(GLMresults))
    for (nm in missing_cols) GLMresults[[nm]] <- NA
    
    # Keep only required columns (and in a stable order)
    GLMresults <- GLMresults[, required_cols, drop = FALSE]
  }
  
  # --- comparisons derived from GLMresults (stable even when 0 rows) ---
  comparisons <- data.frame(
    from = GLMresults$from,
    to = GLMresults$to,
    stringsAsFactors = FALSE
  )
  if (nrow(comparisons) == 0) {
    comparisons$labels <- character(0)
  } else {
    comparisons$labels <- paste0(comparisons$from, "__", comparisons$to)
  }
  
  out <- list(
    condition = conditionVector,
    nCells = nCellsTab,
    GLMresults = GLMresults,
    comparisons = comparisons,
    isGLM = TRUE
  )
  
  if (!is.null(subjectVector)) out$subject <- subjectVector
  if (length(messages) > 0) out$messages <- messages
  
  class(out) <- "SpicyResults"
  out
}


spicyGLM = function(cells,
                    condition,
                    subject = NULL,
                    imageID = "imageID",
                    cellType = "cellType",
                    spatialCoords = c("x", "y"),
                    r,
                    from = NULL,
                    to = NULL,
                    window = "convex",
                    cores = 1,
                    cr2Method = c("fast", "clubSandwich"),
                    fastMethod = c("direct", "dpr1"),
                    estimator = c("firth", "mle"),
                    firthBackend = c("closed_form", "brglm2"),
                    computeLeverage = FALSE) {
  
  cr2Method <- match.arg(cr2Method)
  fastMethod <- match.arg(fastMethod)
  estimator <- match.arg(estimator)
  firthBackend <- match.arg(firthBackend)
  
  # this is a wrapper function
  # check if cells is a dataframe, SingleCellExperiment, or SpatialExperiment
  checkCells(cells)
  
  # check condition is provided
  if (is.null(condition)) {
    stop("Please provide a condition.")
  }
  
  # check condition column exists in data
  checkCondition(cells, condition)
  
  # check cell types exist in data
  if (!is.null(from) && !is.null(to)) {
    if (any(!from %in% getCellType(cells)) ||
        any(!to   %in% getCellType(cells))) {
      stop("`from` or `to` cell types not found in data.")
    }
  }
  
  # check that appropriate window is provided
  if (!window %in% c("rectangle", "convex", "concave")) {
    stop("Invalid value for `window`. Use 'rectangle', 'convex', or 'concave'.")
  }
  
  # check spatialCoords exist in data in appropriate format
  checkCoords(cells, spatialCoords)
  
  # check cores is a number
  if (!is.numeric(cores)) {
    stop("Please provide the number of cores you wish to use.")
  }
  
  # check if subject/image ID have a one-to-one mapping
  df = colData(cells) |> as.data.frame()
  oneToOne = FALSE
  
  if (is.null(subject)) {
    oneToOne = TRUE
    message("No subject ID provided. Clustering by image ID instead of subject. If your data includes ",
            "multiple images from the same subject, provide subject ID to cluster by ",
            "subject instead.")
  } else if (nrow(as.data.frame(unique(df[, subject]))) == nrow(as.data.frame(unique(df[, imageID])))) {
    oneToOne = TRUE
    message("Your specified subject parameter has a one-to-one mapping with imageID. Clustering by image ID instead of subject.")
  } 
  
  # coerce condition to factor and use first level as base group
  if (!is.null(condition)) {
    conditionVector = as.data.frame(getImagePheno(cells, imageID = imageID))[[condition]]
    wasFactor = is.factor(conditionVector)
    
    if (!wasFactor) conditionVector = as.factor(conditionVector)
    conditionVector = droplevels(conditionVector)
    conditionVector = relevel(conditionVector, ref = levels(conditionVector)[1])
    
    cli_inform(paste0(
      if (!wasFactor) "Coercing condition into factor. " else "",
      "Dropping unused levels. Using ",
      condition, " = ", levels(conditionVector)[1],
      " as base comparison group. If this is not the desired base group, please reorder the condition factor."
    ))
  }
  
  base_out <- .new_spicyGLM_result(
    cells = cells,
    condition = condition,
    subject = subject,
    imageID = imageID,
    cellType = cellType
  )
  
  if (!is.null(from) && !is.null(to) && length(from) == 1 && length(to) == 1) {
    cat("Computing pairwise spatial metrics...\n")
    dfPair <- modelDataGen(cells = cells, condition = condition, subject = subject,
                           from = from, to = to, r = r, imageID = imageID,
                           cellType = cellType, spatialCoords = spatialCoords,
                           window = window, cores = 1, oneToOne = oneToOne)
    
    if (is.null(dfPair) || nrow(dfPair) == 0) {
      base_out$messages <- paste0("Pair ", from, "__", to, " skipped: no spatial metrics could be computed.")
      return(base_out)
    }
    
    cat("Fitting GLM model...\n")
    GLMresults <- buildGLM(dfPair, oneToOne = oneToOne, cr2Method = cr2Method,
                           fastMethod = fastMethod, estimator = estimator,
                           firthBackend = firthBackend, computeLeverage = computeLeverage)
    
    if (is.data.frame(GLMresults) && "reason" %in% colnames(GLMresults)) {
      base_out$skipped <- GLMresults
      base_out$messages <- GLMresults$message
      return(base_out)
    }
    
    if (is.null(GLMresults)) {
      base_out$messages <- paste0("Pair ", from, "__", to, " skipped: model not estimable (see warnings for details).")
      return(base_out)
    }
    
    GLMresults$from <- from
    GLMresults$to <- to
    
    GLMresults <- GLMresults |>
      dplyr::select(c("from", "to", "conditionRef", "conditionComp", "coef_ref", "coef_comp",
                      "logRateRatio", "rateRatio", "p.value", "estimator", "mle_would_skip",
                      "mle_skip_reason"))
    
    base_out$GLMresults <- GLMresults
    base_out$comparisons <- data.frame(from = GLMresults$from, to = GLMresults$to,
                                       labels = paste0(GLMresults$from, "__", GLMresults$to),
                                       stringsAsFactors = FALSE)
    
    return(base_out)
    
  } else {
    cat("Computing pairwise spatial metrics...\n")
    dfList = getPairwiseAssoc(cells = cells, condition = condition, subject = subject,
                              from = from, to = to, r = r, imageID = imageID,
                              cellType = cellType, spatialCoords = spatialCoords,
                              window = window, cores = cores)
    
    cat("Fitting GLM models for each cell type pair...\n")
    GLMresults = combineGLM(dfResult = dfList, oneToOne = oneToOne, cores = cores,
                            cr2Method = cr2Method, fastMethod = fastMethod,
                            estimator = estimator, firthBackend = firthBackend,
                            computeLeverage = computeLeverage)
    
    skipped <- attr(GLMresults, "skipped")
    if (!is.null(skipped) && nrow(skipped) > 0) {
      base_out$skipped <- skipped
      base_out$messages <- paste0(nrow(skipped), " pair(s) were skipped. See `$skipped` for details.")
    }
  }
  
  base_out$GLMresults <- GLMresults
  base_out$leverage <- attr(GLMresults, "leverage")
  base_out$comparisons <- data.frame(from = GLMresults$from, to = GLMresults$to,
                                     labels = paste0(GLMresults$from, "__", GLMresults$to),
                                     stringsAsFactors = FALSE)
  
  return(base_out)
}

#' @importFrom BiocParallel bplapply MulticoreParam SerialParam
#' @importFrom dplyr bind_rows
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam
#' @importFrom dplyr bind_rows
modelDataGen = function(cells, 
                        condition,
                        subject = NULL,
                        from, 
                        to, 
                        r, 
                        imageID,
                        cellType,
                        spatialCoords = c("x", "y"),
                        window = "convex", 
                        cores = 1,
                        oneToOne) {
  
  # this function generates pairwise metrics for a single cell type pair across all images
  # format data into a dataframe
  if (is(cells, "SpatialExperiment")) {
    coords = as.data.frame(spatialCoords(cells))
    colnames(coords)[1:2] = c("x", "y")
    
    df = data.frame(imageID = cells[[imageID]],
                    condition = cells[[condition]],
                    cellType = cells[[cellType]],
                    x = coords$x,
                    y = coords$y)
    
    if (!is.null(subject)) {
      df$subject = cells[[subject]]
    }
    
    
  } else if (is(cells, "SingleCellExperiment") | is(cells, "data.frame")) {
    x = spatialCoords[[1]]
    y = spatialCoords[[2]]
    
    df = data.frame(imageID = cells[[imageID]],
                    condition = cells[[condition]],
                    cellType = cells[[cellType]],
                    x = cells[[x]],
                    y = cells[[y]])
    
    if (!is.null(subject)) {
      df$subject = cells[[subject]]
    }
    
  }
  
  # Always compute spatial metrics per image.
  # Subject is kept as a column for downstream clustering (CR2),
  # but we should never mix images inside computeImage().
  dfSplit = split(df, df$imageID)
  
  ## same bplapply/SerialParam fix as combineGLM/getPairwiseAssoc
  worker <- function(dfImg) {
    computeImage(dfImg, r = r, window = window, from = from, to = to)
  }
  
  if (cores > 1) {
    BPPARAM = MulticoreParam(workers = cores)
    resultList = bplapply(dfSplit, worker, BPPARAM = BPPARAM)
  } else {
    resultList = lapply(dfSplit, worker)
  }
  
  # combine results
  finalTable = bind_rows(resultList)
  
  if (nrow(finalTable) == 0) {
    warning(paste0(
      "Skipping pair ", from, "__", to,
      ": all images lacked at least one of the required cell types ",
      "(", from, " or ", to, "), so no spatial associations could be computed."
    ), call. = FALSE)
    return(NULL)
  }
  
  return(finalTable)
  
}

#' Compute pairwise spatial associations between cell types
#' 
#' @param cells A \code{SpatialExperiment}, \code{SingleCellExperiment}, or \code{data.frame}
#' containing single-cell or spatial data with cell metadata and coordinates.
#' @param condition  A character specifying which column in \code{cells} which contains the condition or grouping variable. 
#' @param subject A character specifying which column in \code{cells} which contains the patient/donor ID.
#' @param r Radius around each reference cell to consider for counting neighboring cells.
#' @param imageID A character specifying which column in \code{cells} which contains image/sample ID.
#' @param cellType A character specifying which column in \code{cells} which contains the cell types.
#' @param spatialCoords A character vector of length 2 specifying the columns for x and y coordinates if using a \code{SingleCellExperiment} object.
#' @param from Character vector of reference cell types. If \code{NULL}, all cell types are used.
#' @param to Character vector of target cell types. If \code{NULL}, all cell types are used.
#' @param window Defines the spatial window for each image. Options: "convex", "concave", or "rectangle".
#' @param cores Number of cores to use for parallel computation. 
#' 
#' 
#' @examples
#' \dontrun{
#' kerenSPE = SpatialDatasets::spe_Keren_2018()
#' getPairwiseAssoc(cells = kerenSPE,
#'                  condition = "tumour_type",
#'                  subject = "DONOR_NO",
#'                  from = "CD8_T_cell",
#'                  to = "Tumour",
#'                  imageID = "imageID",
#'                  cellType = "cellType",
#'                  spatialCoords = c("x", "y"),
#'                  r = 50,
#'                  cores = 1)
#' }
#'
#' @export
#' 
#' @importFrom BiocParallel MulticoreParam SerialParam bplapply
#' @importFrom dplyr mutate
#' @return A named list of data frames. Each element corresponds to a cell type pair and contains spatial association metrics for each image. 
#' The names of the list elements are of the form "from__to". Since spicyGLM's
#' log rate ratio is direction-invariant (A__B and B__A give identical results),
#' only one direction per unordered pair is included; self-pairs (A__A) are
#' always included.
getPairwiseAssoc = function(cells,
                            condition,
                            subject = NULL,
                            r = NULL,
                            imageID,
                            cellType,
                            spatialCoords,
                            from = NULL,
                            to = NULL, 
                            window = "convex",
                            cores = 1) { 
  # this function computes pairwise metrics for all images - a wrapper for modelDataGen
  # check if cells is a dataframe, SingleCellExperiment, or SpatialExperiment
  checkCells(cells)
  
  # check condition is provided
  if (is.null(condition)) {
    stop("Please provide a condition.")
  }
  
  # check condition column exists in data
  checkCondition(cells, condition)
  
  # check cell types exist in data
  if (!is.null(from) && !is.null(to)) {
    if (any(!from %in% getCellType(cells)) ||
        any(!to   %in% getCellType(cells))) {
      stop("`from` or `to` cell types not found in data.")
    }
  }
  
  # check that appropriate window is provided
  if (!window %in% c("rectangle", "convex", "concave")) {
    stop("Invalid value for `window`. Use 'rectangle', 'convex', or 'concave'.")
  }
  
  # check spatialCoords exist in data in appropriate format
  checkCoords(cells, spatialCoords)
  
  # check cores is a number
  if (!is.numeric(cores)) {
    stop("Please provide the number of cores you wish to use.")
  }
  
  meta <- colData(cells) |> as.data.frame()
  
  oneToOne = FALSE
  
  # check mapping between image ID and subject
  if (is.null(subject)) {
    oneToOne = TRUE
  } else if (nrow(as.data.frame(unique(meta[, subject]))) == nrow(as.data.frame(unique(meta[, imageID])))) {
    oneToOne = TRUE
  } 
  
  if (!is.null(from) || !is.null(to)) {
    allTypes <- union(from, to)
    cellPairs <- getCellTypePairs(cells = NULL, cellType = NULL, includeSelf = TRUE,
                                  typesOverride = allTypes)
  } else {
    cellPairs = getCellTypePairs(cells, cellType = cellType, includeSelf = TRUE)
  }
  
  if (is.null(r)) {
    r = 50
  }

  ## same bplapply/SerialParam fix as combineGLM
  worker <- function(i) {
    from.i = cellPairs$from[i]
    to.i = cellPairs$to[i]
    
    df = modelDataGen(cells = cells,
                      condition = condition,
                      subject = subject,
                      from = from.i,
                      to = to.i,
                      r = r,
                      imageID = imageID,
                      cellType = cellType,
                      spatialCoords = spatialCoords,
                      window = window,
                      cores = 1,
                      oneToOne = oneToOne) 
    if (is.null(df) || nrow(df) == 0) {
      df = NULL
    } else {
      df = dplyr::mutate(df, from = from.i, to = to.i)
    }
    
    pairName = paste0(from.i, "__", to.i)
    list(pairName = pairName, data = df)
  }
  
  if (cores > 1) {
    BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE)
    resultList = bplapply(seq_len(nrow(cellPairs)), worker, BPPARAM = BPPARAM)
  } else {
    resultList = lapply(seq_len(nrow(cellPairs)), worker)
  }
  
  # convert to named list
  namedList = setNames(
    lapply(resultList, `[[`, "data"),
    sapply(resultList, `[[`, "pairName"))
  
  return(namedList)
}

#' @importFrom spatstat.geom owin convexhull.xy area.owin crosspairs
#' @importFrom concaveman concaveman
computeImage = function(dfImg,
                        from,
                        to,
                        r,
                        window = "convex") {
  # extract image metadata
  img = dfImg$imageID[1]
  conditionImg = dfImg$condition[1]
  typesImg = dfImg$cellType
  coordsImg = dfImg[, c("x", "y")]
  
  # generate unique cell IDs for each cell per type
  typeCounts = ave(seq_along(typesImg), typesImg, FUN = seq_along)
  imgCellID = paste0(typesImg, "_", img, "_", typeCounts)
  
  # TODO: cell-level weights
  
  idxFrom = which(typesImg == from)
  idxTo = which(typesImg == to)
  
  if(length(idxFrom) == 0 || length(idxTo) == 0){
    warning(paste0(
      "Skipping image ", img,
      " for pair ", from, "__", to,
      ": missing required cell type(s) in this image (",
      from, "=", length(idxFrom), ", ",
      to, "=", length(idxTo), ")."
    ), call. = FALSE)
    return(NULL)
  }
  
  # define the spatial window for the image
  if (window == "rectangle") {
    win = owin(xrange = range(coordsImg$x), yrange = range(coordsImg$y))
  } else if (window == "convex") {
    win = convexhull.xy(coordsImg)
  } else if (window == "concave") {
    hullCoords = concaveman(coordsImg)
    win = spatstat.geom::owin(poly = list(x = hullCoords[,1], y = hullCoords[,2]))
  } else {
    stop("Invalid value for `window`. Use 'rectangle', 'convex', or 'concave'.")
  }
  
  # compute density
  areaImg = spatstat.geom::area.owin(win)
  ## should NOT BE idxFrom
  dens = (length(idxTo) / areaImg) * (pi * r^2)
  
  # compute target counts per reference cell -- exact radius search via spatial
  # binning, replaces the dense fields::rdist matrix (validated: identical
  # counts, substantially faster, see benchmark)
  ptsFrom = spatstat.geom::ppp(coordsImg$x[idxFrom], coordsImg$y[idxFrom], window = win, check = FALSE)
  ptsTo   = spatstat.geom::ppp(coordsImg$x[idxTo],   coordsImg$y[idxTo],   window = win, check = FALSE)
  
  cp = spatstat.geom::crosspairs(ptsFrom, ptsTo, rmax = r)
  countsToFrom = tabulate(cp$i, nbins = length(idxFrom))
  
  # assemble result dataframe
  dfResult = data.frame(cellID = as.factor(imgCellID[idxFrom]),
                        from = from,
                        to = to,
                        imageID = as.factor(img),
                        n = countsToFrom,
                        condition = as.factor(conditionImg),
                        density = dens)
  
  if ("subject" %in% colnames(dfImg)) {
    dfResult$subject = as.factor(dfImg$subject[1])
  }
  
  return(dfResult)
}

.new_spicy_skip <- function(from, to, reason, message) {
  data.frame(
    from = as.character(from),
    to = as.character(to),
    reason = as.character(reason),
    message = as.character(message),
    stringsAsFactors = FALSE
  )
}

buildGLM = function(dfResultPairwise,
                    oneToOne,
                    cr2Method = c("fast", "clubSandwich"),
                    fastMethod = c("direct", "dpr1"),
                    estimator = c("firth", "mle"),
                    firthBackend = c("closed_form", "brglm2"),
                    computeLeverage = FALSE) {
  
  cr2Method <- match.arg(cr2Method)
  fastMethod <- match.arg(fastMethod)
  estimator <- match.arg(estimator)
  firthBackend <- match.arg(firthBackend)
  
  from = dfResultPairwise$from |> unique()
  to = dfResultPairwise$to |> unique()
  
  dfResultPairwise$condition = droplevels(factor(dfResultPairwise$condition))
  condition_levels <- levels(dfResultPairwise$condition)
  
  if (length(unique(dfResultPairwise$condition)) < 2) {
    msg <- paste0(
      "Skipping pair ", from, "__", to,
      ": no data available for one condition. '", from, "' or '", to,
      "' cells were not found in any image belonging to at least one ",
      "condition group, so there is no basis for a comparison at all ",
      "(this is different from zero neighbour counts -- the cell types ",
      "themselves are absent from that condition's images)."
    )
    warning(msg, call. = FALSE)
    return(.new_spicy_skip(from, to, reason = "one_condition_level", message = msg))
  }
  
  ## computed regardless of estimator, so Firth-fitted pairs still record
  ## whether MLE would have been non-estimable here
  mle_all_zero <- !any(dfResultPairwise$n > 0, na.rm = TRUE)
  
  pos_by_cond <- tapply(
    dfResultPairwise$n > 0, dfResultPairwise$condition,
    function(x) any(x, na.rm = TRUE)
  )
  mle_zero_conds <- names(pos_by_cond)[!pos_by_cond]
  mle_one_condition_zero <- length(mle_zero_conds) > 0 && !mle_all_zero
  
  mle_would_skip <- mle_all_zero || mle_one_condition_zero
  mle_skip_reason <- if (mle_all_zero) {
    "all_zero"
  } else if (mle_one_condition_zero) {
    "one_condition_zero"
  } else {
    NA_character_
  }
  
  if (estimator == "mle" && mle_all_zero) {
    msg <- paste0(
      "Skipping pair ", from, "__", to, ": neighbour counts are identically zero ",
      "(n = 0 for all reference cells in all images across both conditions). ",
      "Condition effects are non-identifiable under the Poisson log-link, and ",
      "this holds regardless of estimator -- both conditions carry no signal."
    )
    warning(msg, call. = FALSE)
    return(.new_spicy_skip(from, to, reason = "all_zero", message = msg))
  }
  
  if (estimator == "mle" && mle_one_condition_zero) {
    msg <- paste0(
      "Skipping pair ", from, "__", to, ": structural zeros by condition. ",
      "In condition(s) {", paste(mle_zero_conds, collapse = ", "), "} all ",
      "neighbour counts are zero, while other conditions have positive counts. ",
      "This places the Poisson condition contrast on the boundary of the ",
      "parameter space (log rate ratio -> +/-Inf) under ordinary MLE, so the ",
      "effect is non-estimable here. Set estimator = 'firth' to fit this pair ",
      "with Firth's bias-reduced estimator instead, which remains finite in ",
      "this case."
    )
    warning(msg, call. = FALSE)
    return(.new_spicy_skip(from, to, reason = "one_condition_zero", message = msg))
  }
  
  if (estimator == "mle") {
    GLMfit <- fit_mle_glm(dfResultPairwise)
  } else {
    backend <- firthBackend
    if (cr2Method == "clubSandwich" && backend == "closed_form") {
      message(
        "Pair ", from, "__", to, ": cr2Method = 'clubSandwich' requires a real ",
        "fitted model object; switching firthBackend to 'brglm2' for this pair."
      )
      backend <- "brglm2"
    }
    fitResult <- fit_pair(dfResultPairwise, estimator = "firth", backend = backend)
    GLMfit <- fitResult$fit
  }
  
  clusterVec = if (oneToOne) {
    dfResultPairwise$imageID
  } else {
    dfResultPairwise$subject
  }
  
  beta = stats::coef(GLMfit)
  
  if (length(beta) != 2) {
    msg <- paste0(
      "Skipping pair ", from, "__", to,
      ": fit did not return exactly two condition coefficients ",
      "(possible rank deficiency or dropped level)."
    )
    warning(msg, call. = FALSE)
    return(.new_spicy_skip(from, to, reason = "one_condition_dropped", message = msg))
  }
  
  logRR = unname(beta[2]) - unname(beta[1])
  tmp = sub("^condition", "", names(beta))
  
  leverage_tbl <- NULL 
  
  if (cr2Method == "fast") {
    
    cluster_ids <- unique(clusterVec)
    patients <- lapply(cluster_ids, build_patient,
                       cluster_vec = clusterVec,
                       df_result = dfResultPairwise,
                       fit = GLMfit,
                       condition_levels = condition_levels)
    
    group_of_patient <- sapply(patients, function(p) p$group)
    group_sizes <- table(group_of_patient)
    
    if (any(group_sizes < 2)) {
      msg <- paste0(
        "Skipping pair ", from, "__", to,
        ": group(s) ", paste(names(group_sizes)[group_sizes < 2], collapse = ", "),
        " contain only one patient/cluster. CR2's leverage correction is undefined ",
        "for a singleton group, since it works by measuring how much a patient's ",
        "data pulled the fit relative to everyone else in their group - with ",
        "no one else in the group, there's nothing to measure that against."
      )
      warning(msg, call. = FALSE)
      return(.new_spicy_skip(from, to, reason = "one_patient_per_group", message = msg))
    }
    
    V = vcovCR2_fast_multi(patients, method = fastMethod)
    
    leverage_tbl <- NULL
    if (computeLeverage) {
      S_g_vec  <- attr(V, "S_g")
      grp_idx  <- attr(V, "group")
      
      leverage_tbl <- data.frame(
        patient_id = as.character(cluster_ids),
        group      = condition_levels[grp_idx],
        N_i        = sapply(patients, function(p) sum(p$n_ij)),
        T_i        = sapply(patients, patient_total),
        stringsAsFactors = FALSE
      )
      leverage_tbl$S_g       <- S_g_vec[grp_idx]
      leverage_tbl$leverage  <- leverage_tbl$T_i / leverage_tbl$S_g
      leverage_tbl$from      <- from
      leverage_tbl$to        <- to
    }
    
    waldResult = waldTest_CR2_fast(logRR, V, patients, method = fastMethod)
    waldP = waldResult$p.value
    
  } else {
    
    V = vcovClubSandwichCluster(GLMfit, type = "CR2", cluster = clusterVec)
    
    L = matrix(c(-1, 1), nrow = 1)
    waldP = clubSandwich::Wald_test(GLMfit, L, V, tidy = TRUE)$p_val[1] |> as.numeric()
    
  }
  
  out = data.frame(conditionRef  = tmp[1],
                   conditionComp = tmp[2],
                   coef_ref      = unname(beta[1]),
                   coef_comp     = unname(beta[2]),
                   logRateRatio  = logRR,
                   rateRatio     = exp(logRR),
                   p.value       = waldP,
                   estimator     = estimator,
                   mle_would_skip = mle_would_skip,
                   mle_skip_reason = mle_skip_reason)
  
  attr(out, "leverage") <- leverage_tbl
  return(out)
}

#' @importFrom clubSandwich vcovCR
vcovClubSandwichCluster = function(fit,
                             type = c("naive", "CR0", "CR2", "CR3"),
                             cluster) {
  type = match.arg(type)

  if (type == "naive") return(stats::vcov(fit))
  if (missing(cluster)) stop("Provide 'cluster' for clustered vcov.")
  
  clubSandwich::vcovCR(fit, cluster = cluster, type = type)
}

getCellTypePairs = function(cells,
                            cellType = "cellType",
                            includeSelf = TRUE,
                            typesOverride = NULL) {
  
  cellTypes <- if (!is.null(typesOverride)) typesOverride else unique(cells[[cellType]])
  
  crossPairs <- combn(cellTypes, 2, simplify = FALSE)
  crossPairs <- do.call(rbind, lapply(crossPairs, function(x) data.frame(from = x[1], to = x[2])))
  
  if (includeSelf) {
    selfPairs <- data.frame(from = cellTypes, to = cellTypes)
    cellPairs <- rbind(crossPairs, selfPairs)
  } else {
    cellPairs <- crossPairs
  }
  
  cellPairs
}

#' @importFrom BiocParallel MulticoreParam SerialParam bplapply
#' @importFrom dplyr bind_rows
combineGLM = function(dfResult,
                      oneToOne,
                      cores = 1,
                      cr2Method = c("fast", "clubSandwich"),
                      fastMethod = c("direct", "dpr1"),
                      estimator = c("firth", "mle"),
                      firthBackend = c("closed_form", "brglm2"),
                      computeLeverage = FALSE) {
  
  cr2Method <- match.arg(cr2Method)
  fastMethod <- match.arg(fastMethod)
  estimator <- match.arg(estimator)
  firthBackend <- match.arg(firthBackend)
  
  worker <- function(pairName) {
    tryCatch({
      dfPair = dfResult[[pairName]]
      
      if (is.null(dfPair) || nrow(dfPair) == 0) {
        return(NULL)
      }
      
      from_to = strsplit(pairName, "__")[[1]]
      from = from_to[1]
      to = from_to[2]
      
      modelFit = buildGLM(dfPair, oneToOne = oneToOne, cr2Method = cr2Method,
                          fastMethod = fastMethod, estimator = estimator,
                          firthBackend = firthBackend, computeLeverage = computeLeverage)
      
      if (is.null(modelFit)) return(NULL)
      
      if (!("from" %in% names(modelFit))) modelFit$from <- from
      if (!("to" %in% names(modelFit))) modelFit$to <- to
      
      list(fit = modelFit, leverage = attr(modelFit, "leverage"))
      
    }, error = function(e) {
      structure(
        list(pair = pairName, message = conditionMessage(e), call = conditionCall(e)),
        class = "spicy_error"
      )
    })
  }
  
  ## bplapply's promise-evaluation machinery misbehaves for this package even
  ## under SerialParam (cores == 1); plain lapply is equivalent when serial
  ## and avoids the issue entirely. Only route through BiocParallel when
  ## genuine parallelism is requested.
  if (cores > 1) {
    BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE)
    resultList = bplapply(names(dfResult), worker, BPPARAM = BPPARAM)
  } else {
    resultList = lapply(names(dfResult), worker)
  }
  
  fitList      <- lapply(resultList, `[[`, "fit")
  leverageList <- lapply(resultList, `[[`, "leverage")
  
  combined     <- bind_rows(fitList)
  leverage_all <- bind_rows(leverageList)
  
  if (!("reason" %in% colnames(combined))) {
    skipped <- data.frame(from = character(0), to = character(0),
                          reason = character(0), message = character(0),
                          stringsAsFactors = FALSE)
  } else {
    skipped <- combined[!is.na(combined$reason), c("from","to","reason","message"), drop = FALSE]
    combined <- combined[is.na(combined$reason), , drop = FALSE]
  }
  
  combined = combined |> dplyr::select(c("from", "to", "conditionRef", "conditionComp",
                                         "coef_ref", "coef_comp", "logRateRatio", "rateRatio",
                                         "p.value", "estimator", "mle_would_skip", "mle_skip_reason"))
  
  combined = combined |> dplyr::mutate(p.adj = p.adjust(p.value, method = "fdr")) |>
    dplyr::arrange(p.adj)
  
  attr(combined, "skipped") <- skipped
  attr(combined, "leverage") <- leverage_all
  return(combined)
}

.showSpicyGLMResults <- function(object) {
  df <- object$GLMresults
  
  message("Number of cell type pairs tested: ", nrow(df))
  
  if (nrow(df) == 0) {
    message("No valid models were fit (all requested pairs were skipped).")
    if (!is.null(object$messages)) message("Note: ", object$messages)
    return(invisible(NULL))
  }
  
  sigPairs <- sum(df$p.value < 0.05, na.rm = TRUE)
  message("\nNumber of differentially localised cell type pairs: ", sigPairs)
}


setMethod(
  "show", methods::signature(object = "SpicyResults"), function(object) {
    .showSpicyGLMResults(object)
  }
)