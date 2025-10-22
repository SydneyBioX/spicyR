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
                    cores = 1,
                    includeSelf = TRUE) {
  
  if (!is.null(condition)) {
    conditionVector = as.data.frame(getImagePheno(cells))[[condition]]
    wasFactor = is.factor(conditionVector)
    
    if (!wasFactor) conditionVector = as.factor(conditionVector)
    conditionVector = droplevels(conditionVector)
    conditionVector = relevel(conditionVector, ref = levels(conditionVector)[1])
    
    cli::cli_inform(paste0(
      if (!wasFactor) "Coercing condition into factor. " else "",
      "Dropping unused levels. Using ",
      condition, " = ", levels(conditionVector)[1],
      " as base comparison group. If this is not the desired base group, please reorder the condition factor."
    ))
  }
  
  if (!is.null(from) && !is.null(to) && length(from) == 1 && length(to) == 1) {
    cat("Computing pairwise spatial metrics...\n")
    dfPair = modelDataGen(
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
      cores = cores,
      includeSelf = includeSelf)
    
    cat("Fitting GEE models for each cell type pair...\n")
    GEEresults = combineGEE(dfResult = dfList, cores = cores)
  }
  
  spicyGEEResult = list()
  spicyGEEResult$condition = conditionVector
  if (!is.null(subject)) spicyGEEResult$subject = as.data.frame(getImagePheno(cells))[[subject]]
  
  if (exists("dfList")) {
    spicyGEEResult$comparisons = data.frame(
      from = GEEresults$from,
      to   = GEEresults$to,
      labels = paste0(GEEresults$from, "__", GEEresults$to)
    )
  } else {
    spicyGEEResult$comparisons = data.frame(
      from = from,
      to = to,
      labels = paste0(from, "__", to)
    )
  }
  
  spicyGEEResult$nCells = table(getImageID(cells), getCellType(cells))
  spicyGEEResult$GEEresults = GEEresults
  
  return(spicyGEEResult)
}



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

getPairwiseAssoc = function(cells,
                           condition,
                           subject,
                           from = NULL,
                           to = NULL, 
                           r = NULL,
                           imageID,
                           cellType,
                           spatialCoords,
                           window = "convex",
                           cores = 1,
                           includeSelf = TRUE) { 
  
  # get cell type pairs
  if (!is.null(from) && !is.null(to)) {
    # expand all combinations of from × to
    cellPairs = expand.grid(from = from, to = to, stringsAsFactors = FALSE)
  } else {
    # compute all pairs from unique cell types
    cellPairs = getCellTypePairs(cells, cellType = cellType, includeSelf = includeSelf)
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
    win = spatstat.geom::owin(xrange = range(coordsImg$x), yrange = range(coordsImg$y))
  } else if (window == "convex") {
    win = spatstat.geom::convexhull.xy(coordsImg)
  } else if (window == "concave") {
    hullCoords = concaveman::concaveman(coordsImg)
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
    dfResult = data.frame(cellID = if(length(idxFrom) > 0) as.factor(imgCellID[idxFrom]) else NA,
      imageID = as.factor(img),
      subject = as.factor(subjectImg),
      n = NA_integer_,
      condition = as.factor(conditionImg),
      density   = dens)
    
    return(dfResult)
  }
  
}

buildGEE = function(dfResultPairwise) {
  
  dfValid = dfResultPairwise[!is.na(dfResultPairwise$n), ]
  
  if (nrow(dfValid) == 0) {
    warning("No valid data (n not NA) for this pair. Skipping GEE.")
    return(NULL)
  }

  GEEfit = geepack::geeglm(
    n ~ 1 + condition + offset(log(density)),
    id = dfValid$imageID,
    data = dfValid,
    family = poisson("log"),
    corstr = "independence")
  
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

combineGEE = function(dfResult, cores = 1) {
  if (cores > 1) {
    BPPARAM = MulticoreParam(workers = cores)
  } else {
    BPPARAM = SerialParam()
  }
  
  # Closure captures dfResult
  FUN = function(pairName) {
    dfPair = dfResult[[pairName]]
    
    # extract from/to from the name
    from_to = strsplit(pairName, "__")[[1]]
    from = from_to[1]
    to   = from_to[2]
    
    modelFit = buildGEE(dfPair)
    # add pair info
    modelFit$from = from
    modelFit$to   = to
    
    return(modelFit)
  }
  
  # Only one bplapply call, using the closure
  resultList = bplapply(names(dfResult), FUN = FUN, BPPARAM = BPPARAM)
  combinedDF = bind_rows(resultList)
  return(combinedDF)
}
