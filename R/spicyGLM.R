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
#' @param computeDiagnostics Logical; if \code{TRUE}, also computes and attaches
#'   a QC diagnostics table set on the returned object: \code{$diagnostics =
#'   list(pair = <df>, patient = <df>, image = <df>, crossPair = list(patient = <df>,
#'   image = <df>))}. \code{pair} (Table 1) is a one-row-per-pair summary (patient
#'   count, Satterthwaite df, which patient dominates influence/point-estimate shift);
#'   \code{patient} (Table 2) and \code{image} (Table 3) merge pre-fit leverage share,
#'   post-fit CR2 variance-share influence, closed-form leave-one-out point-estimate
#'   shift, and within-pair relative/percentile-rank versions of leverage and
#'   influence, at the patient and image level respectively; \code{crossPair}
#'   (Table 4) aggregates those percentile ranks across every cell-type pair a
#'   patient (or patient/image) appears in, with a Wilson interval on the
#'   proportion flagged, for the leverage and influence pathways independently.
#'   Diagnostics are only assembled when \code{cr2Method = "fast"} and
#'   \code{estimator = "firth"} with \code{firthBackend = "closed_form"} --
#'   the point-estimate-shift formula is the exact closed-form Firth
#'   solution, and influence depends on the fast CR2 sandwich machinery, so
#'   other estimator/backend/cr2Method combinations leave \code{$diagnostics}
#'   \code{NULL} for that pair. Adds some compute cost per pair; \code{FALSE}
#'   by default.
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
#'   \item{diagnostics}{If \code{computeDiagnostics = TRUE}, \code{list(pair = <df>,
#'     patient = <df>, image = <df>, crossPair = list(patient = <df>, image = <df>))}
#'     of QC diagnostics across all pairs -- see \code{computeDiagnostics} above.}
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
  nCellsTab <- table(getImageID(cells, imageID = imageID), getCellType(cells, cellType = cellType))
  
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
                    computeDiagnostics = FALSE) {

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
    if (any(!from %in% getCellType(cells, cellType = cellType)) ||
        any(!to   %in% getCellType(cells, cellType = cellType))) {
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
  
  cellTypePresence <- computeCellTypePresence(
    cells = cells, condition = condition, imageID = imageID, cellType = cellType
  )
  
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
                           window = window, cores = 1, oneToOne = oneToOne, cellTypePresence = cellTypePresence)
    
    if (nrow(dfPair) == 0) {
      base_out$messages <- attr(dfPair, "skipMessage")
      base_out$skipped <- .new_spicy_skip(from, to, reason = attr(dfPair, "skipReason"),
                                          message = attr(dfPair, "skipMessage"))
      return(base_out)
    }
    
    cat("Fitting GLM model...\n")
    GLMresults <- buildGLM(dfPair, oneToOne = oneToOne, subject = subject, cr2Method = cr2Method,
                           fastMethod = fastMethod, estimator = estimator,
                           firthBackend = firthBackend, computeDiagnostics = computeDiagnostics,
                           cellTypePresence = cellTypePresence)
    
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
    
    diagnostics_tbl <- attr(GLMresults, "diagnostics")
    if (!is.null(diagnostics_tbl)) {
      diagnostics_tbl$patient <- computeRelativeDiagnostics(diagnostics_tbl, level = "patient")
      diagnostics_tbl$image   <- computeRelativeDiagnostics(diagnostics_tbl, level = "image")
      diagnostics_tbl$crossPair <- list(
        patient = crossPairDiagnostics(diagnostics_tbl$patient),
        image   = crossPairDiagnostics(diagnostics_tbl$image)
      )
    }

    GLMresults <- GLMresults |>
      dplyr::select(c("from", "to", "conditionRef", "conditionComp", "coef_ref", "coef_comp",
                      "logRateRatio", "rateRatio", "p.value", "estimator", "mle_would_skip",
                      "mle_skip_reason"))

    base_out$GLMresults <- GLMresults
    base_out$diagnostics <- diagnostics_tbl
    base_out$comparisons <- data.frame(from = GLMresults$from, to = GLMresults$to,
                                       labels = paste0(GLMresults$from, "__", GLMresults$to),
                                       stringsAsFactors = FALSE)
    
    return(base_out)
    
  } else {
    cat("Computing pairwise spatial metrics...\n")
    dfList = getPairwiseAssoc(cells = cells, condition = condition, subject = subject,
                              from = from, to = to, r = r, imageID = imageID,
                              cellType = cellType, spatialCoords = spatialCoords,
                              window = window, cores = cores, cellTypePresence = cellTypePresence)
    
    cat("Fitting GLM models for each cell type pair...\n")
    GLMresults = combineGLM(dfResult = dfList, oneToOne = oneToOne, subject = subject, cores = cores,
                            cr2Method = cr2Method, fastMethod = fastMethod,
                            estimator = estimator, firthBackend = firthBackend,
                            computeDiagnostics = computeDiagnostics, cellTypePresence = cellTypePresence)
    
    skipped <- attr(GLMresults, "skipped")
    if (!is.null(skipped) && nrow(skipped) > 0) {
      base_out$skipped <- skipped
      base_out$messages <- paste0(nrow(skipped), " pair(s) were skipped. See `$skipped` for details.")
    }
  }
  
  base_out$GLMresults <- GLMresults
  diagnostics_tbl <- attr(GLMresults, "diagnostics")
  if (!is.null(diagnostics_tbl)) {
    diagnostics_tbl$patient <- computeRelativeDiagnostics(diagnostics_tbl, level = "patient")
    diagnostics_tbl$image   <- computeRelativeDiagnostics(diagnostics_tbl, level = "image")
    diagnostics_tbl$crossPair <- list(
      patient = crossPairDiagnostics(diagnostics_tbl$patient),
      image   = crossPairDiagnostics(diagnostics_tbl$image)
    )
  }
  base_out$diagnostics <- diagnostics_tbl
  base_out$comparisons <- data.frame(from = GLMresults$from, to = GLMresults$to,
                                     labels = paste0(GLMresults$from, "__", GLMresults$to),
                                     stringsAsFactors = FALSE)
  
  return(base_out)
}

#' Diagnose why a cell type pair has no usable data for one or more condition groups
#'
#' Shared by modelDataGen() (fires before any spatial metrics exist) and
#' buildGLM() (fires when a specific condition group's rows are entirely
#' absent after the fact). Distinguishes three real causes: `from` absent,
#' `to` absent, or both individually present but never co-occurring in the
#' same image within that group.
#'
#' @param from,to Character; the reference/target cell types for this pair.
#' @param missingConditions Character vector of condition group(s) with no
#'   usable rows for this pair.
#' @param cellTypePresence The presence table from computeCellTypePresence().
#' @param context "spatial" (modelDataGen) or "comparison" (buildGLM) --
#'   only changes the closing clause of the message.
#'
#' @return A list: reason (a short code) and message (the full skip text).
diagnoseMissingConditions <- function(from, to, missingConditions, cellTypePresence,
                                      context = c("spatial", "comparison")) {
  context <- match.arg(context)
  closing <- if (context == "spatial") {
    "No spatial associations could be computed for this group."
  } else {
    "No basis for a comparison in this group."
  }
  
  perGroup <- lapply(missingConditions, function(g) {
    fromRow <- cellTypePresence[cellTypePresence$cellType == from & cellTypePresence$condition == g, ]
    toRow   <- cellTypePresence[cellTypePresence$cellType == to   & cellTypePresence$condition == g, ]
    
    fromPresent <- fromRow$n_images_with_type > 0
    toPresent   <- toRow$n_images_with_type > 0
    nTotal      <- fromRow$n_images_total
    
    if (!fromPresent && !toPresent) {
      list(reason = "both_absent", text = sprintf(
        'condition level "%s" has 0 images containing "%s" and 0 images containing "%s". %s',
        g, from, to, closing
      ))
    } else if (!fromPresent) {
      list(reason = "type_absent", text = sprintf(
        'condition level "%s" has 0 images where "%s" is present at all (0 of %d images in this group contain any %s cells). %s',
        g, from, nTotal, from, closing
      ))
    } else if (!toPresent) {
      list(reason = "type_absent", text = sprintf(
        'condition level "%s" has 0 images where "%s" is present at all (0 of %d images in this group contain any %s cells). %s',
        g, to, nTotal, to, closing
      ))
    } else {
      list(reason = "no_cooccurrence", text = sprintf(
        'condition level "%s" has both cell types present somewhere (%s in %d of %d images, %s in %d of %d images), but no single image in this group contains both simultaneously. %s',
        g, from, fromRow$n_images_with_type, nTotal, to, toRow$n_images_with_type, nTotal, closing
      ))
    }
  })
  
  reasons <- unique(vapply(perGroup, `[[`, character(1), "reason"))
  reasonCode <- if (length(reasons) == 1) reasons else "mixed"
  
  list(
    reason = reasonCode,
    message = paste0(
      "Skipping pair ", from, "__", to, ": ",
      paste(vapply(perGroup, `[[`, character(1), "text"), collapse = " ")
    )
  )
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
                        oneToOne,
                        cellTypePresence) {
  
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
    # every condition group lacked usable data for this pair -- diagnose all of them
    allConditions <- unique(cellTypePresence$condition)
    diag <- diagnoseMissingConditions(from, to, allConditions, cellTypePresence, context = "spatial")
    warning(diag$message, call. = FALSE)
    
    emptyResult <- finalTable
    attr(emptyResult, "skipReason")  <- diag$reason
    attr(emptyResult, "skipMessage") <- diag$message
    return(emptyResult)
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
                            cores = 1,
                            cellTypePresence) { 
  # this function computes pairwise metrics for all images - a wrapper for modelDataGen
  # check if cells is a dataframe, SingleCellExperiment, or SpatialExperiment
  checkCells(cells)
  
  # check condition is provided
  if (is.null(condition)) {
    stop("Please provide a condition.")
  }
  
  # check condition column exists in data
  checkCondition(cells, condition)
  
  if (!is.null(from) && !is.null(to)) {
    if (any(!from %in% getCellType(cells, cellType = cellType)) ||
        any(!to   %in% getCellType(cells, cellType = cellType))) {
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
                      oneToOne = oneToOne,
                      cellTypePresence = cellTypePresence)
    
    if (is.null(df)) {
      df = NULL
    } else if (nrow(df) == 0) {
      # keep df as-is (zero rows, but tagged with skipReason/skipMessage)
      df = df
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

#' Compute pre-fit CR2 leverage shares from raw counts and density
#'
#' Leverage share is exactly pre-fit computable (the `exp(beta_hat)` term
#' cancels in the T_i / S_g ratio within a condition group), so this needs
#' no model fit -- only per-image reference-cell counts and the density
#' offset from \code{computeImage()}.
#'
#' @param dfResultPairwise Data frame for one cell-type pair, one row per
#'   reference cell, with columns \code{imageID}, \code{condition},
#'   \code{density}, and \code{subject} if \code{oneToOne} is \code{FALSE}.
#' @param oneToOne Logical; if \code{TRUE} cluster by \code{imageID},
#'   otherwise by \code{subject}.
#' @param subject Unused; kept for signature consistency with
#'   \code{buildGLM()}.
#'
#' @return A list with \code{patient} (columns \code{patient_id}, \code{group},
#'   \code{T_i}, \code{S_g}, \code{l_i}) and \code{image} (columns
#'   \code{patient_id}, \code{image_id}, \code{group}, \code{n_ij},
#'   \code{density_ij}, \code{l_ij}) data frames.
computeLeverage <- function(dfResultPairwise, oneToOne, subject = NULL) {

  patient_id <- if (oneToOne) {
    as.character(dfResultPairwise$imageID)
  } else {
    as.character(dfResultPairwise$subject)
  }

  df <- data.frame(
    patient_id = patient_id,
    image_id   = as.character(dfResultPairwise$imageID),
    group      = as.character(dfResultPairwise$condition),
    density    = dfResultPairwise$density,
    stringsAsFactors = FALSE
  )

  image_tbl <- df |>
    dplyr::group_by(patient_id, image_id) |>
    dplyr::summarise(
      group      = dplyr::first(group),
      n_ij       = dplyr::n(),
      density_ij = dplyr::first(density),
      .groups = "drop"
    )

  patient_tbl <- image_tbl |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(
      group = dplyr::first(group),
      T_i   = sum(n_ij * density_ij),
      .groups = "drop"
    )

  group_tbl <- patient_tbl |>
    dplyr::group_by(group) |>
    dplyr::summarise(S_g = sum(T_i), .groups = "drop")

  patient_tbl <- patient_tbl |>
    dplyr::left_join(group_tbl, by = "group") |>
    dplyr::mutate(l_i = T_i / S_g) |>
    as.data.frame()

  image_tbl <- image_tbl |>
    dplyr::left_join(
      dplyr::select(patient_tbl, patient_id, T_i),
      by = "patient_id"
    ) |>
    dplyr::mutate(l_ij = (n_ij * density_ij) / T_i) |>
    dplyr::select(patient_id, image_id, group, n_ij, density_ij, l_ij) |>
    as.data.frame()

  list(patient = patient_tbl, image = image_tbl)
}

#' Compute closed-form leave-one-out shift in the fitted log rate ratio
#'
#' No CR2 machinery required -- pure aggregation of raw counts/densities,
#' exploiting the closed-form Firth solution \code{beta_g = log((Y_g+0.5)/D_g)}.
#' \code{GLMfit} is used only to verify the recomputed group coefficients
#' against the actual fit, not for any computation.
#'
#' @param dfResultPairwise Data frame for one cell-type pair, one row per
#'   reference cell, with columns \code{n}, \code{density}, \code{condition},
#'   \code{imageID}, and \code{subject} if \code{oneToOne} is \code{FALSE}.
#' @param oneToOne Logical; if \code{TRUE} cluster by \code{imageID},
#'   otherwise by \code{subject}.
#' @param subject Unused; kept for signature consistency with
#'   \code{computeLeverage()}/\code{buildGLM()}.
#' @param GLMfit The fitted model object for this pair, used only for the
#'   step-5 correctness check.
#'
#' @return A data frame with columns \code{patient_id}, \code{group},
#'   \code{y_i}, \code{d_i}, \code{delta_i}.
computePointEstimateShift <- function(dfResultPairwise, oneToOne, subject = NULL, GLMfit) {

  patient_id <- if (oneToOne) {
    as.character(dfResultPairwise$imageID)
  } else {
    as.character(dfResultPairwise$subject)
  }

  group_levels <- levels(droplevels(factor(dfResultPairwise$condition)))

  df <- data.frame(
    patient_id = patient_id,
    group      = as.character(dfResultPairwise$condition),
    n          = dfResultPairwise$n,
    density    = dfResultPairwise$density,
    stringsAsFactors = FALSE
  )

  patient_tbl <- df |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(
      group = dplyr::first(group),
      y_i   = sum(n),
      d_i   = sum(density),
      .groups = "drop"
    ) |>
    as.data.frame()

  group_tbl <- patient_tbl |>
    dplyr::group_by(group) |>
    dplyr::summarise(Y_g = sum(y_i), D_g = sum(d_i), .groups = "drop")

  Y_g <- stats::setNames(group_tbl$Y_g, group_tbl$group)[group_levels]
  D_g <- stats::setNames(group_tbl$D_g, group_tbl$group)[group_levels]

  beta_g_hat <- log((Y_g + 0.5) / D_g)
  names(beta_g_hat) <- paste0("condition", group_levels)

  fit_beta <- stats::coef(GLMfit)
  if (!isTRUE(all.equal(as.numeric(beta_g_hat), as.numeric(fit_beta)))) {
    stop(
      "computePointEstimateShift(): recomputed group coefficients do not match ",
      "coef(GLMfit) -- aggregation of n/density by patient/group disagrees with ",
      "the actual fit; do not trust delta_i until this is resolved.",
      call. = FALSE
    )
  }

  delta_i <- vapply(seq_len(nrow(patient_tbl)), function(idx) {
    g   <- patient_tbl$group[idx]
    y_i <- patient_tbl$y_i[idx]
    d_i <- patient_tbl$d_i[idx]

    # sole-contributor patients hit log(0.5 / 0) here -- the same degenerate
    # boundary as l_i -> 1 / CR2 singularity for a group of size 1.
    if (y_i == Y_g[[g]] && d_i == D_g[[g]]) {
      return(NA_real_)
    }

    beta_g_hat_minus_i <- log((Y_g[[g]] - y_i + 0.5) / (D_g[[g]] - d_i))
    unname(beta_g_hat[[paste0("condition", g)]] - beta_g_hat_minus_i)
  }, numeric(1))

  data.frame(
    patient_id = patient_tbl$patient_id,
    group      = patient_tbl$group,
    y_i        = patient_tbl$y_i,
    d_i        = patient_tbl$d_i,
    delta_i    = delta_i,
    stringsAsFactors = FALSE
  )
}

#' Assemble patient-level (Table 2) and image-level (Table 3) QC diagnostics
#'
#' Joins the outputs of \code{computeLeverage()}, \code{computeInfluence()},
#' and \code{computePointEstimateShift()} into two merged tables. Verifies
#' the three sources agree on patient identity and group assignment before
#' merging -- a disagreement means one of the three functions is using a
#' different patient/group derivation than the others.
#'
#' @param leverage_result \code{list(patient=, image=)} from \code{computeLeverage()}.
#' @param influence_result \code{list(patient=, image=)} from \code{computeInfluence()}.
#' @param shift_result Data frame from \code{computePointEstimateShift()}.
#'
#' @return \code{list(patient = <Table 2>, image = <Table 3>)}.
assembleDiagnosticsTables <- function(leverage_result, influence_result, shift_result) {

  ids_leverage  <- sort(leverage_result$patient$patient_id)
  ids_influence <- sort(influence_result$patient$patient_id)
  ids_shift     <- sort(shift_result$patient_id)
  stopifnot(identical(ids_leverage, ids_influence), identical(ids_leverage, ids_shift))

  group_leverage  <- stats::setNames(leverage_result$patient$group, leverage_result$patient$patient_id)
  group_influence <- stats::setNames(influence_result$patient$group, influence_result$patient$patient_id)
  group_shift     <- stats::setNames(shift_result$group, shift_result$patient_id)
  stopifnot(
    identical(unname(group_leverage[ids_leverage]), unname(group_influence[ids_leverage])),
    identical(unname(group_leverage[ids_leverage]), unname(group_shift[ids_leverage]))
  )

  n_i_tbl <- leverage_result$image |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(n_i = sum(n_ij), .groups = "drop")

  residual_sums_tbl <- influence_result$image |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(
      raw_residual_sum      = sum(raw_residual_sum_ij),
      adjusted_residual_sum = sum(adjusted_residual_sum_ij),
      .groups = "drop"
    )

  table2 <- leverage_result$patient |>
    dplyr::inner_join(n_i_tbl, by = "patient_id") |>
    dplyr::inner_join(influence_result$patient, by = c("patient_id", "group")) |>
    dplyr::inner_join(residual_sums_tbl, by = "patient_id") |>
    dplyr::inner_join(shift_result, by = c("patient_id", "group")) |>
    dplyr::select(patient_id, group, n_i, T_i, S_g, l_i, raw_residual_sum,
                  adjusted_residual_sum, e_i, influence_i, y_i, d_i, delta_i) |>
    as.data.frame()

  n_ij_leverage <- leverage_result$image |>
    dplyr::arrange(patient_id, image_id) |>
    dplyr::pull(n_ij)
  n_ij_influence <- influence_result$image |>
    dplyr::arrange(patient_id, image_id) |>
    dplyr::pull(n_ij)
  stopifnot(isTRUE(all.equal(n_ij_leverage, n_ij_influence)))

  table3 <- leverage_result$image |>
    dplyr::inner_join(
      dplyr::select(influence_result$image, -n_ij),
      by = c("patient_id", "image_id", "group")
    ) |>
    dplyr::select(patient_id, image_id, group, n_ij, density_ij, l_ij,
                  raw_residual_sum_ij, adjusted_residual_sum_ij, e_ij, influence_ij) |>
    as.data.frame()

  list(patient = table2, image = table3)
}

#' Assemble the one-row pair-level QC summary (Table 1)
#'
#' @param waldResult Return value of \code{waldTest_CR2_fast()}; needs \code{$df} (nu).
#' @param table2 The merged patient-level table from \code{assembleDiagnosticsTables()}.
#'
#' @return A one-row data frame: \code{n_patients | nu | max_influence |
#'   patient_with_max_influence | max_abs_delta_logRR | patient_with_max_delta_logRR}.
computePairSummary <- function(waldResult, table2) {
  nu <- waldResult$df
  n_patients <- nrow(table2)

  max_influence <- max(table2$influence_i)
  patient_with_max_influence <- table2$patient_id[which.max(table2$influence_i)]

  max_abs_delta_logRR <- max(abs(table2$delta_i), na.rm = TRUE)
  patient_with_max_delta_logRR <- table2$patient_id[which.max(abs(table2$delta_i))]

  data.frame(
    n_patients = n_patients,
    nu = nu,
    max_influence = max_influence,
    patient_with_max_influence = patient_with_max_influence,
    max_abs_delta_logRR = max_abs_delta_logRR,
    patient_with_max_delta_logRR = patient_with_max_delta_logRR,
    stringsAsFactors = FALSE
  )
}

#' Compute within-pair relative leverage/influence and percentile ranks
#'
#' Rescales pre-fit leverage share and post-fit variance-share influence
#' against their own natural comparison scope (leverage: within-group or
#' within-patient; influence: pair-wide), so values become comparable across
#' cell-type pairs with different patient/image counts. Computes both the
#' leverage and influence pathways in one pass, since \code{crossPairDiagnostics()}
#' needs both.
#'
#' @param diagnostics The stacked, multi-pair \code{list(pair=, patient=, image=)}
#'   from \code{buildGLM()}/\code{combineGLM()}'s \code{computeDiagnostics = TRUE} path.
#' @param level \code{"patient"} or \code{"image"}; selects which of
#'   \code{diagnostics$patient} / \code{diagnostics$image} to operate on.
#'
#' @return A data frame -- the input to \code{crossPairDiagnostics()}.
computeRelativeDiagnostics <- function(diagnostics, level = c("patient", "image")) {
  level <- match.arg(level)

  if (level == "patient") {
    out <- diagnostics$patient |>
      dplyr::group_by(from, to, group) |>
      dplyr::mutate(
        n_patients_per_group = dplyr::n(),
        rel_l_i = l_i * n_patients_per_group,
        percentile_rank_leverage = dplyr::percent_rank(rel_l_i)
      ) |>
      dplyr::group_by(from, to) |>
      dplyr::mutate(
        n_patients_in_pair = dplyr::n(),
        rel_influence_i = influence_i * n_patients_in_pair,
        percentile_rank_influence = dplyr::percent_rank(rel_influence_i)
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        percentile_rank_leverage  = ifelse(is.nan(percentile_rank_leverage),  NA_real_, percentile_rank_leverage),
        percentile_rank_influence = ifelse(is.nan(percentile_rank_influence), NA_real_, percentile_rank_influence)
      ) |>
      dplyr::select(patient_id, group, from, to, n_patients_per_group, n_patients_in_pair,
                    l_i, influence_i, delta_i, rel_l_i, percentile_rank_leverage,
                    rel_influence_i, percentile_rank_influence,
                    n_i, T_i, S_g, raw_residual_sum, adjusted_residual_sum, e_i, y_i, d_i) |>
      as.data.frame()
  } else {
    out <- diagnostics$image |>
      dplyr::group_by(from, to, patient_id) |>
      dplyr::mutate(
        n_images_for_patient = dplyr::n(),
        rel_l_ij = l_ij * n_images_for_patient,
        percentile_rank_leverage = dplyr::percent_rank(rel_l_ij)
      ) |>
      dplyr::group_by(from, to) |>
      dplyr::mutate(
        n_images_in_pair = dplyr::n(),
        rel_influence_ij = influence_ij * n_images_in_pair,
        percentile_rank_influence = dplyr::percent_rank(rel_influence_ij)
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        percentile_rank_leverage  = ifelse(is.nan(percentile_rank_leverage),  NA_real_, percentile_rank_leverage),
        percentile_rank_influence = ifelse(is.nan(percentile_rank_influence), NA_real_, percentile_rank_influence)
      ) |>
      dplyr::select(patient_id, image_id, group, from, to,
                    n_images_for_patient, n_images_in_pair,
                    l_ij, influence_ij, rel_l_ij, percentile_rank_leverage,
                    rel_influence_ij, percentile_rank_influence,
                    density_ij, raw_residual_sum_ij, adjusted_residual_sum_ij, e_ij) |>
      as.data.frame()
  }

  out
}

#' Aggregate relative diagnostics across cell-type pairs (Table 4)
#'
#' For each patient (or patient/image), across every cell-type pair it
#' appears in, computes how often it lands in the top \code{topPercent} of
#' the within-pair percentile rank, for both the leverage and influence
#' pathways independently, with a Wilson score interval around each
#' proportion to discount thin evidence (few pairs present).
#'
#' @param relative_table Output of \code{computeRelativeDiagnostics()}.
#'   Patient- vs. image-level is inferred from whether an \code{image_id}
#'   column is present.
#' @param topPercent Numeric in (0,1); the top fraction of the within-pair
#'   percentile rank counted as "flagged". Default 0.05.
#'
#' @return A data frame, one row per patient (or patient/image), sorted by
#'   \code{wilson_lower_influence} descending.
#'
#' @importFrom binom binom.wilson
crossPairDiagnostics <- function(relative_table, topPercent = 0.05) {

  isImageLevel <- "image_id" %in% names(relative_table)
  groupVars <- if (isImageLevel) c("patient_id", "image_id") else "patient_id"

  result <- relative_table |>
    dplyr::mutate(
      flagged_influence = percentile_rank_influence >= (1 - topPercent),
      flagged_leverage  = percentile_rank_leverage  >= (1 - topPercent)
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(groupVars))) |>
    dplyr::summarise(
      n_pairs_present = dplyr::n(),
      n_pairs_flagged_influence = sum(flagged_influence, na.rm = TRUE),
      n_pairs_flagged_leverage  = sum(flagged_leverage,  na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      prop_flagged_influence = n_pairs_flagged_influence / n_pairs_present,
      prop_flagged_leverage  = n_pairs_flagged_leverage  / n_pairs_present
    )

  wilson_inf <- binom::binom.wilson(result$n_pairs_flagged_influence, result$n_pairs_present)
  wilson_lev <- binom::binom.wilson(result$n_pairs_flagged_leverage,  result$n_pairs_present)

  result$wilson_lower_influence <- wilson_inf$lower
  result$wilson_upper_influence <- wilson_inf$upper
  result$wilson_lower_leverage  <- wilson_lev$lower
  result$wilson_upper_leverage  <- wilson_lev$upper

  result |>
    dplyr::arrange(dplyr::desc(wilson_lower_influence)) |>
    as.data.frame()
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

#' Compute per-cell-type presence across condition groups
#'
#' For each cell type and each condition group, records how many images in
#' that group contain at least one cell of that type, out of how many images
#' total belong to that group. This is the shared diagnostic used by both
#' modelDataGen() and buildGLM() to explain *why* a pair has no testable
#' data in a given condition group, rather than just reporting that it doesn't.
#'
#' @param cells A SpatialExperiment, SingleCellExperiment, or data.frame.
#' @param condition Character; condition/grouping column name.
#' @param imageID Character; image/sample ID column name.
#' @param cellType Character; cell type column name.
#'
#' @return A data frame with columns: cellType, condition, n_images_with_type,
#'   n_images_total. One row per (cellType, condition) combination, including
#'   zero rows for cell types entirely absent from a condition group.
computeCellTypePresence <- function(cells, condition, imageID, cellType) {
  
  df <- data.frame(
    imageID   = as.character(cells[[imageID]]),
    condition = as.character(cells[[condition]]),
    cellType  = as.character(cells[[cellType]]),
    stringsAsFactors = FALSE
  )
  
  # total images per condition group
  imagesPerCondition <- df |>
    dplyr::distinct(imageID, condition) |>
    dplyr::count(condition, name = "n_images_total")
  
  # images per (cellType, condition) that actually contain that type
  imagesWithType <- df |>
    dplyr::distinct(cellType, imageID, condition) |>
    dplyr::count(cellType, condition, name = "n_images_with_type")
  
  # full grid: every cell type x every condition, filling zeros where a
  # cell type never appears in a given condition group at all
  full <- expand.grid(
    cellType  = unique(df$cellType),
    condition = unique(df$condition),
    stringsAsFactors = FALSE
  )
  
  full |>
    dplyr::left_join(imagesWithType, by = c("cellType", "condition")) |>
    dplyr::mutate(n_images_with_type = ifelse(is.na(n_images_with_type), 0L, n_images_with_type)) |>
    dplyr::left_join(imagesPerCondition, by = "condition")
}

buildGLM = function(dfResultPairwise,
                    oneToOne,
                    subject = NULL,
                    cr2Method = c("fast", "clubSandwich"),
                    fastMethod = c("direct", "dpr1"),
                    estimator = c("firth", "mle"),
                    firthBackend = c("closed_form", "brglm2"),
                    computeDiagnostics = FALSE,
                    cellTypePresence) {

  cr2Method <- match.arg(cr2Method)
  fastMethod <- match.arg(fastMethod)
  estimator <- match.arg(estimator)
  firthBackend <- match.arg(firthBackend)
  
  from = dfResultPairwise$from |> unique()
  to = dfResultPairwise$to |> unique()
  
  dfResultPairwise$condition = droplevels(factor(dfResultPairwise$condition))
  condition_levels <- levels(dfResultPairwise$condition)
  
  if (length(unique(dfResultPairwise$condition)) < 2) {
    allConditions <- unique(cellTypePresence$condition)
    missingConditions <- setdiff(allConditions, condition_levels)
    
    diag <- diagnoseMissingConditions(from, to, missingConditions, cellTypePresence,
                                      context = "comparison")
    warning(diag$message, call. = FALSE)
    return(.new_spicy_skip(from, to, reason = diag$reason, message = diag$message))
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
  
  diagnostics <- NULL

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
      failingGroups <- names(group_sizes)[group_sizes < 2]
      unit <- if (!is.null(subject)) "patient" else "image"
      
      perGroupMsgs <- vapply(failingGroups, function(gIdx) {
        gLabel <- condition_levels[as.integer(gIdx)]
        idx <- which(group_of_patient == as.integer(gIdx))
        
        survivorText <- if (length(idx) == 1) {
          p <- patients[[idx]]
          survivorID <- as.character(cluster_ids[idx])
          nImg <- length(p$n_ij)
          sprintf(
            'only %s "%s" has any data for this pair (%d qualifying image%s)',
            unit, survivorID, nImg, if (nImg == 1) "" else "s"
          )
        } else {
          sprintf("no %s has any data for this pair", unit)
        }
        
        fromRow <- cellTypePresence[cellTypePresence$cellType == from & cellTypePresence$condition == gLabel, ]
        toRow   <- cellTypePresence[cellTypePresence$cellType == to   & cellTypePresence$condition == gLabel, ]
        
        sprintf(
          'in "%s", %s -- every other %s in this group has zero images where %s and %s co-occur. For context, %s appears in %d of %d images in this group, and %s appears in %d of %d.',
          gLabel, survivorText, unit, from, to,
          from, fromRow$n_images_with_type, fromRow$n_images_total,
          to, toRow$n_images_with_type, toRow$n_images_total
        )
      }, character(1))
      
      subjectCaveat <- if (is.null(subject)) {
        paste0(
          " Note: no `subject` was provided, so images could not be grouped by patient -- ",
          "if these images do in fact come from multiple donors, that's fine, but if they ",
          "share a donor, consider passing `subject` for correct clustering."
        )
      } else {
        ""
      }
      
      msg <- paste0(
        "Skipping pair ", from, "__", to, ": ",
        paste(perGroupMsgs, collapse = " "),
        " CR2's leverage correction needs at least two independent ", unit, "s per group ",
        "to measure how much any one ", unit, "'s data pulled the fit relative to the rest ",
        "of the group; with only one, there's nothing to compare it against.",
        subjectCaveat
      )
      warning(msg, call. = FALSE)
      return(.new_spicy_skip(from, to, reason = "one_patient_per_group", message = msg))
    }
    
    V = vcovCR2_fast_multi(patients, method = fastMethod)

    waldResult = waldTest_CR2_fast(logRR, V, patients, method = fastMethod)
    waldP = waldResult$p.value

  } else {
    
    V = vcovClubSandwichCluster(GLMfit, type = "CR2", cluster = clusterVec)
    
    L = matrix(c(-1, 1), nrow = 1)
    waldP = clubSandwich::Wald_test(GLMfit, L, V, tidy = TRUE)$p_val[1] |> as.numeric()
    
  }

  if (computeDiagnostics && cr2Method == "fast" && estimator == "firth" && firthBackend == "closed_form") {
    leverage_result <- computeLeverage(dfResultPairwise, oneToOne, subject)
    influence_result <- computeInfluence(patients, waldResult, condition_levels, V)
    shift_result <- computePointEstimateShift(dfResultPairwise, oneToOne, subject, GLMfit)

    merged <- assembleDiagnosticsTables(leverage_result, influence_result, shift_result)
    table1 <- computePairSummary(waldResult, merged$patient)
    table1 <- cbind(data.frame(from = from, to = to, stringsAsFactors = FALSE), table1)

    merged$patient$from <- from;  merged$patient$to <- to
    merged$image$from <- from;    merged$image$to <- to

    diagnostics <- list(pair = table1, patient = merged$patient, image = merged$image)
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

  attr(out, "diagnostics") <- diagnostics
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
                      subject = NULL,
                      cores = 1,
                      cr2Method = c("fast", "clubSandwich"),
                      fastMethod = c("direct", "dpr1"),
                      estimator = c("firth", "mle"),
                      firthBackend = c("closed_form", "brglm2"),
                      computeDiagnostics = FALSE,
                      cellTypePresence) {

  cr2Method <- match.arg(cr2Method)
  fastMethod <- match.arg(fastMethod)
  estimator <- match.arg(estimator)
  firthBackend <- match.arg(firthBackend)
  
  worker <- function(pairName) {
    tryCatch({
      dfPair = dfResult[[pairName]]
      from_to = strsplit(pairName, "__")[[1]]
      from = from_to[1]
      to = from_to[2]
      
      if (is.null(dfPair)) {
        return(NULL)  
      }
      
      if (nrow(dfPair) == 0) {
        skipMessage <- attr(dfPair, "skipMessage")
        skipReason  <- attr(dfPair, "skipReason")
        if (is.null(skipMessage)) {
          # defensive fallback in case attributes are missing for any reason
          skipMessage <- paste0("Pair ", from, "__", to, ": no spatial metrics could be computed.")
          skipReason  <- "no_spatial_metrics"
        }
        skipRow <- .new_spicy_skip(from, to, reason = skipReason, message = skipMessage)
        return(list(fit = skipRow, diagnostics = NULL))
      }
      
      modelFit = buildGLM(dfPair, oneToOne = oneToOne, subject = subject, cr2Method = cr2Method,
                          fastMethod = fastMethod, estimator = estimator,
                          firthBackend = firthBackend, computeDiagnostics = computeDiagnostics,
                          cellTypePresence = cellTypePresence)
      
      if (is.null(modelFit)) return(NULL)
      
      if (!("from" %in% names(modelFit))) modelFit$from <- from
      if (!("to" %in% names(modelFit))) modelFit$to <- to
      
      list(fit = modelFit, diagnostics = attr(modelFit, "diagnostics"))
      
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
  
  fitList         <- lapply(resultList, `[[`, "fit")
  diagnosticsList <- lapply(resultList, `[[`, "diagnostics")

  combined <- bind_rows(fitList)
  diagnostics_all <- list(
    pair    = bind_rows(lapply(diagnosticsList, `[[`, "pair")),
    patient = bind_rows(lapply(diagnosticsList, `[[`, "patient")),
    image   = bind_rows(lapply(diagnosticsList, `[[`, "image"))
  )
  
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
  attr(combined, "diagnostics") <- diagnostics_all
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