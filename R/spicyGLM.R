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
#' 
#' @return A list with the following elements:
#' \describe{
#'   \item{condition}{Factor vector of the condition used in the GEE models.}
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
  required_cols <- c("from", "to", "conditionRef", "conditionComp", "coef", "rateRatio", "p.value")
  if (is.null(GLMresults)) {
    GLMresults <- data.frame(
      from = character(0),
      to = character(0),
      conditionRef = character(0),
      conditionComp = character(0),
      coef = numeric(0),
      rateRatio = numeric(0),
      p.value = numeric(0),
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
                    cores = 1) {
  
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
    message("No subject ID provided. Clustering by image ID instead of subject.")
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
  
  # compute metrics and build model for single cell type pair
  if (!is.null(from) && !is.null(to) && length(from) == 1 && length(to) == 1) {
    cat("Computing pairwise spatial metrics...\n")
    dfPair <- modelDataGen(
      cells = cells,
      condition = condition,
      subject = subject,
      from = from,
      to = to,
      r = r,
      imageID = imageID,
      cellType = cellType,
      spatialCoords = spatialCoords,
      window = window,
      cores = 1,
      oneToOne = oneToOne
    )
      
    if (is.null(dfPair) || nrow(dfPair) == 0) {
      base_out$messages <- paste0("Pair ", from, "__", to, " skipped: no spatial metrics could be computed.")
      return(base_out)
    }
    
    cat("Fitting GLM model...\n")
    GLMresults <- buildGLM(dfPair, oneToOne = oneToOne)
    
    # buildGLM may return a "skip record" data.frame instead of results
    if (is.data.frame(GLMresults) && "reason" %in% colnames(GLMresults)) {
      base_out$skipped <- GLMresults
      base_out$messages <- GLMresults$message
      return(base_out)
    }
    
    # if something truly returned NULL (rare will not happen - also means skip record might not be working), keep the old fallback
    if (is.null(GLMresults)) {
      base_out$messages <- paste0(
        "Pair ", from, "__", to,
        " skipped: model not estimable (see warnings for details)."
      )
      return(base_out)
    }
    
    
    GLMresults$from <- from
    GLMresults$to <- to
    GLMresults <- GLMresults |>
      dplyr::select(c("from", "to", "conditionRef", "conditionComp", "coef", "rateRatio", "p.value"))
    
    # Fill the base object with results (keeps schema consistent)
    base_out$GLMresults <- GLMresults
    base_out$comparisons <- data.frame(
      from = GLMresults$from,
      to = GLMresults$to,
      labels = paste0(GLMresults$from, "__", GLMresults$to),
      stringsAsFactors = FALSE
    )
    
    return(base_out)
    
    # compute metrics and build model for multiple cell type pairs
  } else {
    cat("Computing pairwise spatial metrics...\n")
    dfList = getPairwiseAssoc(
      cells = cells,
      condition = condition,
      subject = subject,
      from = from, 
      to = to,
      r = r,
      imageID = imageID,
      cellType = cellType,
      spatialCoords = spatialCoords,
      window = window,
      cores = cores)
    
    cat("Fitting GLM models for each cell type pair...\n")
    GLMresults = combineGLM(dfResult = dfList, oneToOne = oneToOne, cores = cores)
    
    skipped <- attr(GLMresults, "skipped")
    if (!is.null(skipped) && nrow(skipped) > 0) {
      base_out$skipped <- skipped
      base_out$messages <- paste0(nrow(skipped), " pair(s) were skipped. See `$skipped` for details.")
    }
  } 
  
  # Fill base_out with multi-pair results
  base_out$GLMresults <- GLMresults
  base_out$comparisons <- data.frame(
    from = GLMresults$from,
    to = GLMresults$to,
    labels = paste0(GLMresults$from, "__", GLMresults$to),
    stringsAsFactors = FALSE
  )
  
  return(base_out)
}

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
                        cores = cores,
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
  
  
  if (cores > 1) {
    # parallel processing using parallel::mclapply (Unix only)
    BPPARAM = MulticoreParam(workers = cores)
  } else {
    BPPARAM = SerialParam()
  }
  
  # compute pairwise metrics for each image
  resultList = bplapply(dfSplit, 
                        FUN = computeImage, 
                        r = r,
                        window = window,
                        from = from,
                        to = to,
                        BPPARAM = BPPARAM)
  
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
#' @return A named list of data frames. Each element corresponds to a cell type pair and contains spatial association metrics for each image. 
#' The names of the list elements are of the form "from__to".
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
  
  # get cell type pairs
  if (!is.null(from) || !is.null(to)) {
    # expand all combinations of from × to
    cellPairs = expand.grid(from = from, to = to, stringsAsFactors = FALSE)
  } else {
    # compute all pairs from unique cell types
    cellPairs = getCellTypePairs(cells, cellType = cellType, includeSelf = TRUE)
  }
  
  if (is.null(r)) {
    r = 50
  }
  
  # parallelisation
  if (cores > 1) {
    BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE)
  } else {
    BPPARAM = SerialParam(progressbar = TRUE)
  }
  
  # compute for all pairs
  resultList = bplapply(
    seq_len(nrow(cellPairs)),
    FUN = function(i) {
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
    }, BPPARAM = BPPARAM)
  
  # convert to named list
  namedList = setNames(
    lapply(resultList, `[[`, "data"),
    sapply(resultList, `[[`, "pairName"))
  
  return(namedList)
}

#' @importFrom spatstat.geom owin convexhull.xy area.owin
#' @importFrom concaveman concaveman
#' @importFrom fields rdist
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
  
  # need image ID
  #subjectImg = if ("subject" %in% colnames(dfImg)) {
  #  dfImg$subject[1]
  #} else {
  #  NA
  #}
      
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
  dens = (length(idxFrom) / areaImg) * (pi * r^2)
  
  # compute target counts per reference cell
  D = fields::rdist(as.matrix(coordsImg[idxFrom, ]), as.matrix(coordsImg[idxTo, ]))
  countsToFrom = rowSums(D <= r)
    
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

#' @importFrom fixest feglm
buildGLM = function(dfResultPairwise, oneToOne) {
  # fit GLM model for one cell type pair
  from = dfResultPairwise$from |> unique()
  to = dfResultPairwise$to |> unique()
  
  dfResultPairwise$condition = droplevels(factor(dfResultPairwise$condition))
  
  if (length(unique(dfResultPairwise$condition)) < 2) {
    msg <- paste(
      "Skipping pair", from, "__", to,
      ": only one condition level after dropping unused levels"
    )
    warning(msg, call. = FALSE)
    return(.new_spicy_skip(from, to, reason = "one_condition_level", message = msg))
  }
  
  if (!any(dfResultPairwise$n > 0, na.rm = TRUE)) {
    msg <- paste0(
      "Skipping pair ", from, "__", to, ": neighbour counts are identically zero ",
      "(n = 0 for all reference cells in all images across both conditions). ",
      "Condition effects are non-identifiable under the Poisson log-link."
    )
    warning(msg, call. = FALSE)
    return(.new_spicy_skip(from, to, reason = "all_zero", message = msg))
  }
  
  
  pos_by_cond <- tapply(
    dfResultPairwise$n > 0,
    dfResultPairwise$condition,
    function(x) any(x, na.rm = TRUE)
  )
  
  if (any(!pos_by_cond)) {
    zero_conds <- names(pos_by_cond)[!pos_by_cond]
    
    msg <- paste0(
      "Skipping pair ", from, "__", to, ": structural zeros by condition. ",
      "In condition(s) {", paste(zero_conds, collapse = ", "), "} all neighbour counts are zero, ",
      "while other conditions have positive counts. ",
      "This places the Poisson condition contrast on the boundary of the parameter space ",
      "(log rate ratio -> +/-Inf), so the condition effect is non-estimable and ",
      "Wald/CR2 inference is invalid."
    )
    
    warning(msg, call. = FALSE)
    return(.new_spicy_skip(from, to, reason = "one_condition_zero", message = msg))
  }
  
  
  GLMfit = fixest::feglm(n ~ 0 + condition,
                         offset = log(dfResultPairwise$density),
                         family = "poisson",
                         data = dfResultPairwise,
                         data.save = TRUE)
  
  clusterVec = if (oneToOne) {
    GLMfit$data$imageID
  } else {
    GLMfit$data$subject
  }
  
  V = vcovClubSandwichCluster(GLMfit, type = "CR2", cluster = clusterVec)
  
  beta = stats::coef(GLMfit)
  
  if (length(beta) != 2) {
    msg <- paste0(
      "Skipping pair ", from, "__", to,
      ": GLM fit did not return exactly two condition coefficients ",
      "(possible rank deficiency or dropped level)."
    )
    warning(msg, call. = FALSE)
    return(.new_spicy_skip(from, to, reason = "one_condition_dropped", message = msg))
  }
  
  
  
  L = matrix(c(-1, 1), nrow = 1)
  
  waldP = clubSandwich::Wald_test(GLMfit, L, V, tidy = TRUE)$p_val[1] |> as.numeric()
  
  # log rate ratio?
  logRR = unname(beta[2]) - unname(beta[1])
  tmp = sub("^condition", "", names(beta))
  out = data.frame(conditionRef = tmp[1],
                   conditionComp = tmp[2],
                   coef = logRR,
                   rateRatio = exp(logRR),
                   p.value = waldP)
  
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
                            includeSelf = TRUE) { # whether to include self pairs, like T cell--T cell
  
  # get unique cell types
  cellTypes = unique(cells[[cellType]])
  
  # generate pairwise combinations
  if (includeSelf) {
    cellPairs = expand.grid(from = cellTypes, to = cellTypes, stringsAsFactors = FALSE)
  } else {
    combs = combn(cellTypes, 2, simplify = FALSE)
    cellPairs = do.call(rbind, lapply(combs, function(x) data.frame(from = x[1], to = x[2])))
  }
  
  return(cellPairs)
}

#' @importFrom BiocParallel MulticoreParam SerialParam bplapply
#' @importFrom dplyr bind_rows
combineGLM = function(dfResult,
                      oneToOne,
                      cores = 1) {
  
  # fit GLM model for all cell type pairs - wrapper for buildGLM
  BPPARAM = if (cores > 1) MulticoreParam(workers = cores, progressbar = TRUE) else SerialParam(progressbar = TRUE)
  
<<<<<<< HEAD
  resultList = bplapply(names(dfResult), function(pairName) {
    dfPair = dfResult[[pairName]]
    
    from_to = strsplit(pairName, "__")[[1]]
    from = from_to[1]
    to = from_to[2]
    
    modelFit = buildGLM(dfPair, oneToOne = oneToOne)
    
    if (is.null(modelFit)) {
      return(NULL)
    }
    
    modelFit$from = from
    modelFit$to = to
    
    return(modelFit)
  }, BPPARAM = BPPARAM)
    
=======
  resultList = bplapply(
    names(dfResult),
    function(pairName) {
      
      tryCatch({
        dfPair = dfResult[[pairName]]
        
        if (is.null(dfPair) || nrow(dfPair) == 0) {
          return(NULL)
        }
        
        from_to = strsplit(pairName, "__")[[1]]
        from = from_to[1]
        to = from_to[2]
        
        modelFit = buildGLM(dfPair, oneToOne = oneToOne)
        
        if (is.null(modelFit)) return(NULL)
        
        if (!("from" %in% names(modelFit))) modelFit$from <- from
        if (!("to" %in% names(modelFit))) modelFit$to <- to
    
        modelFit
        
      }, error = function(e) {
        structure(
          list(
            pair = pairName,
            message = conditionMessage(e),
            call = conditionCall(e)
          ),
          class = "spicy_error"
        )
      })
    },
    BPPARAM = BPPARAM
  )
>>>>>>> b27909f431b26ab75a3d1a9bb9b8b53f2dcda3c8
  
  combined = bind_rows(resultList)
  
  # Split successful fits vs skipped-pair records (returned by buildGLM)
  if (!("reason" %in% colnames(combined))) {
    # no skips occurred, create empty skipped table
    skipped <- data.frame(from = character(0), to = character(0),
                          reason = character(0), message = character(0),
                          stringsAsFactors = FALSE)
  } else {
    skipped <- combined[!is.na(combined$reason), c("from","to","reason","message"), drop = FALSE]
    combined <- combined[is.na(combined$reason), , drop = FALSE]
  }
  
  
  combined = combined |> dplyr::select(c("from", "to", "conditionRef", "conditionComp", "coef", "rateRatio", "p.value"))
  
  combined = combined |> dplyr::mutate(p.adj = p.adjust(p.value, method = "fdr")) |>
    dplyr::arrange(p.adj)
  
  attr(combined, "skipped") <- skipped

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