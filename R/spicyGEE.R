#' `Calculates pairwise spatial associations between cell types across images
#' and fit generalized estimating equation (GEE) models to test for condition effects.
#' 
#' @param cells A \code{SpatialExperiment}, \code{SingleCellExperiment}, or \code{data.frame}
#' containing single-cell or spatial data with cell metadata and coordinates.
#' @param condition  A character specifying which column in \code{cells} which contains the condition or grouping variable. 
#' @param subject A character specifying which column in \code{cells} which contains the patient/donor ID.
#' @param imageID A character specifying which column in \code{cells} which contains image/sample ID.
#' @param cellType A character specifying which column in \code{cells} which contains the cell types.
#' @param spatialCoords Character vector of length 2 specifying the columns for x and y coordinates.
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
#'   \item{GEEresults}{Data frame of GEE model results for each cell type pair.}
#' }
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' kerenSPE = SpatialDatasets::spe_Keren_2018()
#' spicyResult = spicyGEE(
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
spicyGEE = function(cells,
                    condition, 
                    subject,
                    imageID = "imageID",
                    cellType = "cellType",
                    spatialCoords = c("x", "y"),
                    r, 
                    from = NULL,
                    to = NULL,
                    window = "convex",
                    cores = 1) {
  
  if (!is.null(condition)) {
    conditionVector = as.data.frame(getImagePheno(cells))[[condition]]
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
  
  if (!is.null(from) && !is.null(to) && length(from) == 1 && length(to) == 1) {
    cat("Computing pairwise spatial metrics...\n")
    dfPair = modelDataGen(cells = cells,
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
    
    cat("Fitting GEE model...\n")
    GEEresults = buildGEE(dfPair)
    GEEresults$from = from
    GEEresults$to = to
    
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
    
    cat("Fitting GEE models for each cell type pair...\n")
    GEEresults = combineGEE(dfResult = dfList, cores = cores)
  }
  
  spicyGEEResult = list()
  spicyGEEResult$condition = conditionVector
  if (!is.null(subject)) spicyGEEResult$subject = as.data.frame(getImagePheno(cells))[[subject]]
  
  spicyGEEResult$comparisons = data.frame(from = from,
      to = to,
      labels = paste0(from, "__", to))
  
  spicyGEEResult$nCells = table(getImageID(cells), getCellType(cells))
  spicyGEEResult$GEEresults = GEEresults
  
  return(spicyGEEResult)
}

#' @importFrom BiocParallel bplapply MulticoreParam SerialParam
#' @importFrom dplyr bind_rows
modelDataGen = function(cells, 
                        condition,
                        subject,
                        from, 
                        to, 
                        r, 
                        imageID,
                        cellType,
                        spatialCoords = c("x", "y"),
                        window = "convex", 
                        cores = cores) {
  
  
  # format data into a dataframe
  if (is(cells, "SpatialExperiment")) {
    coords = as.data.frame(spatialCoords(cells))
    colnames(coords)[1:2] = c("x", "y")
    
    df = data.frame(subject  = cells[[subject]],
                    imageID  = cells[[imageID]],
                    condition = cells[[condition]],
                    cellType = cells[[cellType]],
                    x = coords$x,
                    y = coords$y)
    
  } else if (is(cells, "SingleCellExperiment") | is(cells, "data.frame")) {
    x = spatialCoords[[1]]
    y = spatialCoords[[2]]
    df = data.frame(subject  = cells[[subject]],
                    imageID  = cells[[imageID]],
                    condition = cells[[condition]],
                    cellType = cells[[cellType]],
                    x = cells[[x]],
                    y = cells[[y]])
  }
  
  # compute per image
  dfSplit = split(df, df$imageID)
  
  if (cores > 1) {
    # parallel processing using parallel::mclapply (Unix only)
    BPPARAM = MulticoreParam(workers = cores)
  } else {
    BPPARAM = SerialParam()
  }
    
  resultList = bplapply(dfSplit, 
                        FUN = computeImage, 
                        r = r,
                        window = window,
                        from = from,
                        to = to,
                        BPPARAM = BPPARAM)
  
  # combine results
  finalTable = bind_rows(resultList)
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
#' @param spatialCoords Character vector of length 2 specifying the columns for x and y coordinates.
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
                           subject,
                           r = NULL,
                           imageID,
                           cellType,
                           spatialCoords,
                           from = NULL,
                           to = NULL, 
                           window = "convex",
                           cores = 1) { 
  
  # get cell type pairs
  if (!is.null(from) && !is.null(to)) {
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
    BPPARAM = MulticoreParam(workers = cores)
  } else {
    BPPARAM = SerialParam()
  }
  
  # compute for all pairs
  resultList = bplapply(
    seq_len(nrow(cellPairs)),
    FUN = function(i) {
      from.i = cellPairs$from[i]
      to.i   = cellPairs$to[i]
      
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
                        cores = 1) |>
        mutate(from = from.i, to = to.i)
      
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
  subjectImg = dfImg$subject[1]
  conditionImg = dfImg$condition[1]
  typesImg = dfImg$cellType
  coordsImg = dfImg[, c("x", "y")]
  
  # generate unique cell IDs for each cell per type
  typeCounts = ave(seq_along(typesImg), typesImg, FUN = seq_along)
  imgCellID = paste0(typesImg, "_", img, "_", typeCounts)
  
  # TODO: cell-level weights
  
  idxFrom = which(typesImg == from)
  idxTo = which(typesImg == to)
  
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
  if (length(idxFrom) > 0 && length(idxTo) > 0) {
    D = fields::rdist(as.matrix(coordsImg[idxFrom, ]), as.matrix(coordsImg[idxTo, ]))
    countsToFrom = rowSums(D <= r)
    
    # assemble result dataframe
    dfResult = data.frame(cellID = as.factor(imgCellID[idxFrom]),
      imageID = as.factor(img),
      subject = as.factor(subjectImg),
      n = countsToFrom,
      condition = as.factor(conditionImg),
      density = dens)
    
    return(dfResult)
    
  } else {
    cat(sprintf("No reference or target cells found for image %s\n", img))
  }
  
  
}

#' @importFrom geepack geeglm
buildGEE = function(dfResultPairwise) {
  
  dfValid = dfResultPairwise[!is.na(dfResultPairwise$n), ]
  
  if (nrow(dfValid) == 0) {
    message("No valid data for this pair. Skipping GEE.")
    return(NULL)
  }
  
  from = dfResultPairwise$from |> unique()
  to = dfResultPairwise$to |> unique()

  GEEfit = tryCatch({
    geepack::geeglm(n ~ 1 + condition + offset(log(density)),
      id = dfValid$imageID,
      data = dfValid,
      family = poisson("log"),
      corstr = "independence")
  }, error = function(e) {
    message(paste("Model fitting failed due to insufficient number of cells:", from, "__", to, "-", e$message))
    return(NULL)
  })
  
  if (is.null(GEEfit)) return(NULL)
  
  coefs = summary(GEEfit)$coefficients
  
  modelFit = data.frame(estimate = coefs[2, "Estimate"],
    std.err = coefs[2, "Std.err"],
    wald = coefs[2, "Wald"],
    p.value = coefs[2, "Pr(>|W|)"])
  
  return(modelFit)
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
combineGEE = function(dfResult, 
                      cores = 1) {
  BPPARAM = if (cores > 1) MulticoreParam(workers = cores) else SerialParam()
  
  resultList = bplapply(names(dfResult), function(pairName) {
    message("Starting pair: ", pairName)
    dfPair = dfResult[[pairName]]
    
    from_to = strsplit(pairName, "__")[[1]]
    from = from_to[1]
    to = from_to[2]
    
    modelFit = buildGEE(dfPair)
    modelFit$from = from
    modelFit$to = to
  
    message("Finished pair: ", pairName)
    return(modelFit)
  }, BPPARAM = BPPARAM)
  
  bind_rows(resultList)
}