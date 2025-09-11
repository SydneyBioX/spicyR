# check if alternativeResults is kontextual
isKontextual <- function(kontextualResult) {
  
    return("kontexutal"%in% names(kontextualResult))
}

#' A format SummarizedExperiment and data.frame objects into a canonical form.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr %>%
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom methods is
#' @importFrom dplyr select any_of rename mutate group_by ungroup
#' @importFrom cli cli_abort cli_inform
#' @noRd
.format_data <- function(
    cells, imageIDCol, cellTypeCol, spatialCoordCols, verbose = FALSE) {
    if (is(cells, "data.frame")) {
        # pass
    } else if (is(cells, "SpatialExperiment")) {
        spatialCoords <- names(spatialCoords(cells))
        cells <- cells %>%
            colData() %>%
            data.frame() %>%
            dplyr::select(-dplyr::any_of(c("x", "y"))) %>%
            cbind(spatialCoords(cells) %>% data.frame())
    } else if (is(cells, "SingleCellExperiment")) {
        cells <- cells %>%
            colData() %>%
            data.frame()
    } else {
        temp <- tryCatch(
            expr = {
                as.data.frame(cells)
            }, error = NULL
        )
        if (is.null(temp)) {
            cli::cli_abort(
                c(
                    "x" = "{.var cells} is an unsupported class: {.cls {class(cells)}}. \n", # nolint
                    "i" = "data.frame (or coercible), SingleCellExperiment and SpatialExperiment are currently supported." # nolint
                )
            )
        }
        cells <- temp
    }

    for (col in c(
        imageIDCol = imageIDCol,
        cellTypeCol = cellTypeCol,
        spatialCoordCols_x = spatialCoordCols[1],
        spatialCoordCols_y = spatialCoordCols[2]
    )) {
        if (!col %in% colnames(cells)) {
            cli::cli_abort(c(
                "x" = "Specified {.var {names(col)}} ({.emph {col}}) is not in {.var cells}.", # nolint
                "i" = "{.code names(cells)}: {names(cells)}"
            ))
        }
    }
  
  spatialCoords <- names(spatialCoords(cells))
  cells <- cells %>%
    colData() %>%
    data.frame() %>%
    dplyr::select(-dplyr::any_of(c("x", "y"))) %>%
    cbind(spatialCoords(cells) %>% data.frame())
  
    needed <- cells %>%
      select(!!imageIDCol, !!cellTypeCol, spatialCoordCols[1], spatialCoordCols[2]) %>%
      dplyr::rename(
        imageID = !!imageIDCol,
        cellType = !!cellTypeCol,
        x = spatialCoordCols[1],
        y = spatialCoordCols[2]
      )
  
    cells <- cells %>%
      dplyr::select(-dplyr::any_of(c("imageID", "cellType", "x", "y"))) %>%
        cbind(needed)
    
    # Check if imageID is a factor, if not make it as factor whilst keeping order
    if(class(cells$imageID) != "factor") {
      cells <- cells %>%
        dplyr::mutate(imageID = factor(imageID, levels = unique(imageID)))
    }
    
    # Check if imageID is a factor, if not make it as factor whilst keeping order
    if(class(cells$cellType) != "factor") {
      cells <- cells %>%
        dplyr::mutate(cellType = factor(cellType, levels = unique(cellType)))
    }
    
    # Order the cells by imageID. This is important for when the data is split downstream
    cells <- cells[order(cells$imageID),]  
    
    
    # create cellID if it does not exist
    if (is.null(cells$cellID)) {
        if (verbose) {
            cli::cli_inform(
                "No column called cellID. Creating one."
            )
        }
        cells <- dplyr::mutate(
            cells,
            cellID = paste0("cell", "_", dplyr::row_number())
        )
    }

    # create imageCellID if it does not exist
    if (is.null(cells$imageCellID)) {
        if (verbose) {
            cli::cli_inform(
                "No column called imageCellID. Creating one."
            )
        }
        cells <- cells %>%
            dplyr::group_by(imageID) %>%
            dplyr::mutate(
                imageCellID = paste0(imageID, "_", dplyr::row_number())
            ) %>%
            dplyr::ungroup()
    }
    
    cells
}

getCellSummary <- function(
    data,
    imageID = NULL,
    bind = TRUE) {
    data %>%
        dplyr::filter(
            if (!is.null(!!imageID)) imageID == !!imageID else TRUE
        ) %>%
        dplyr::select(imageID, cellID, imageCellID, x, y, cellType) %>%
        # dplyr::mutate(imageID = factor(imageID, levels = unique(imageID))) %>%
        S4Vectors::DataFrame() %>%
        {
            if (bind) . else S4Vectors::split(., .$imageID)
        }
}


getImageID <- function(x, imageID = NULL) {
    x %>%
        dplyr::filter(
            if (!is.null(!!imageID)) imageID == !!imageID else TRUE
        ) %>%
        dplyr::pull(imageID)
}

getCellType <- function(x, imageID = NULL) {
    x %>%
        dplyr::filter(
            if (!is.null(!!imageID)) imageID == !!imageID else TRUE
        ) %>%
        dplyr::pull(cellType)
}

getImagePheno <- function(x,
                          imageID = NULL,
                          bind = TRUE,
                          expand = FALSE) {
    x <- x[!duplicated(x$imageID),]
    x <- x %>%
      dplyr::mutate(imageID = factor(imageID, levels = unique(imageID)))
    # x <- x %>%
    #     dplyr::filter(
    #         if (!is.null(!!imageID)) imageID == !!imageID else TRUE
    #     ) %>%
    #     dplyr::select(-cellID, -imageCellID, -x, -y, -cellType) %>%
    #     dplyr::mutate(imageID = as.factor(imageID)) %>%
    #     {
    #         if (expand) . else dplyr::distinct(.)
    #     } %>%
    #     S4Vectors::DataFrame()
    # if (expand) rownames(x) <- x$imageID
    x
}

#' A function to handle validity of argumemts/check for deprecated arguments
#'
#' @importFrom lifecycle deprecate_warn
#' @importFrom methods is
#' @noRd
argumentChecks = function(function_name, user_vals) {
  
  rlang::local_options(lifecycle_verbosity = "warning")  
  
  # handle deprecated arguments
  handle_deprecated = function(old_arg, new_arg, user_vals) {
    if (old_arg %in% names(user_vals) && !(new_arg %in% names(user_vals))) {
      deprecate_warn("1.18.0", paste0(function_name, "(", old_arg, ")"), paste0(function_name, "(", new_arg, ")"))
      assign(new_arg, user_vals[[old_arg]], envir = sys.frame(sys.parent(1)))
    }
  }
  
  handle_deprecated("imageIDCol", "imageID", user_vals)
  handle_deprecated("cellTypeCol", "cellType", user_vals)
  handle_deprecated("spatialCoordCols", "spatialCoords", user_vals)
  handle_deprecated("nCores", "cores", user_vals)
  handle_deprecated("BPPARAM", "cores", user_vals)
  handle_deprecated("Rs", "r", user_vals)
  
  # enforce mutually exclusive arguments
  check_exclusive = function(arg_set, user_vals) {
    provided_args = intersect(arg_set, names(user_vals)) 
    if (length(provided_args) > 1) {
      stop(paste("Please specify only one of", paste(shQuote(arg_set), collapse = ", "), "\n"))
    }
  }
  
  check_exclusive(c("cellTypeCol", "cellType"), user_vals)
  check_exclusive(c("imageIDCol", "imageID"), user_vals)
  check_exclusive(c("spatialCoordCols", "spatialCoords"), user_vals)
  check_exclusive(c("cores", "nCores", "BPPARAM"), user_vals)
  check_exclusive(c("Rs", "r"), user_vals)

  # validity checks for cores/nCores/BPPARAM
  if ("nCores" %in% names(user_vals)) {
    deprecate_warn("1.18.0", paste0(function_name, "(nCores)"), paste0(function_name, "(cores)"))
    
    if (is(user_vals$nCores, "numeric")) {
      assign("cores", simpleSeg:::generateBPParam(cores = user_vals$nCores), envir = parent.frame())
    } else if (is(user_vals$nCores, "MulticoreParam") || is(user_vals$nCores, "SerialParam")) {
      assign("cores", user_vals$nCores, envir = sys.frame(sys.parent(1)))
    } else {
      stop("'nCores' must be either a numeric value, or a MulticoreParam or SerialParam object.\n")
    }
    
  } else if ("BPPARAM" %in% names(user_vals)) {
    deprecate_warn("1.18.0", paste0(function_name, "(BPPARAM)"), paste0(function_name, "(cores)"))
    
    if (is(user_vals$BPPARAM, "MulticoreParam") || is(user_vals$BPPARAM, "SerialParam")) {
      assign("cores", user_vals$BPPARAM, envir = sys.frame(sys.parent(1)))
    } else {
      stop("'BBPARAM' must be a MulticoreParam or SerialParam object.")
    }
    
  } else if ("cores" %in% names(user_vals)) {
    if (is(user_vals$cores, "numeric")) {
      assign("cores", simpleSeg:::generateBPParam(cores = user_vals$cores))
    } else if (is(user_vals$cores, "MulticoreParam") || is(user_vals$cores, "SerialParam")) {
      assign("cores", user_vals$cores, sys.frame(sys.parent(1)))
    } else {
      stop("'cores'  must be either a numeric value, or a MulticoreParam or SerialParam object.\n")
    }
  }
}
