# check if alternativeResults is kontextual
isKontextual <- function(kontextualResult) {

    return("kontexutal"%in% names(kontextualResult))
}

## ---- optional dependencies ---------------------------------------------------------------------

## Stop with an informative error unless the suggested package `pkg` is installed.
.need <- function(pkg, why = NULL) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package '", pkg, "' is needed", if (!is.null(why)) paste0(" ", why) else "",
             ". Install it with BiocManager::install(\"", pkg, "\").", call. = FALSE)
    }
    invisible(TRUE)
}

## methods::is() for an S4 input class from a suggested package. The package that defines the object's
## class is loaded first (if installed), so that is() can see the class hierarchy.
.is_class <- function(x, what) {
    if (isS4(x)) {
        pkg <- attr(class(x), "package")
        if (!is.null(pkg) && !isNamespaceLoaded(pkg)) requireNamespace(pkg, quietly = TRUE)
    }
    methods::is(x, what)
}

## colData() of a SummarizedExperiment (or subclass) as a data.frame, as data.frame(colData(x)) gives it.
.col_data <- function(x) {
    .need("SummarizedExperiment", "to use SummarizedExperiment-based inputs")
    data.frame(SummarizedExperiment::colData(x))
}

## ---- parallel helpers --------------------------------------------------------------------------

## Number of workers from `cores` (a number) or a BiocParallel param object (`cores` or `BPPARAM`), for
## backward compatibility with code that passed BiocParallel params.
.n_workers <- function(cores = 1, BPPARAM = NULL) {
    x <- if (!is.null(BPPARAM)) BPPARAM else cores
    if (is.null(x)) return(1L)
    if (is.numeric(x)) return(max(1L, as.integer(x[1])))
    if (.is_class(x, "BiocParallelParam")) {
        .need("BiocParallel", "to use a BiocParallel param object; or pass `cores` as a number")
        return(max(1L, as.integer(BiocParallel::bpnworkers(x))))
    }
    stop("`cores` must be a number (or a BiocParallel param object).", call. = FALSE)
}

## lapply()/mapply() over `cores` forked workers (parallel::mclapply), sequential with one core or on
## Windows.
.par_lapply <- function(X, FUN, ..., cores = 1L) {
    cores <- .n_workers(cores)
    if (cores > 1L && .Platform$OS.type != "windows" && length(X) > 1L) {
        out <- parallel::mclapply(X, FUN, ..., mc.cores = cores)
        .stop_on_worker_error(out)
    } else {
        lapply(X, FUN, ...)
    }
}

.par_mapply <- function(FUN, ..., MoreArgs = NULL, SIMPLIFY = FALSE, cores = 1L) {
    cores <- .n_workers(cores)
    if (cores > 1L && .Platform$OS.type != "windows") {
        out <- parallel::mcmapply(FUN, ..., MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY, mc.cores = cores)
        .stop_on_worker_error(out)
    } else {
        mapply(FUN, ..., MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY)
    }
}

## mclapply() returns errors as "try-error" values; raise the first, as bplapply() did.
.stop_on_worker_error <- function(out) {
    if (is.list(out)) {
        bad <- vapply(out, inherits, logical(1), what = "try-error")
        if (any(bad)) stop(attr(out[[which(bad)[1]]], "condition"))
    }
    out
}

## ---- small data helpers ------------------------------------------------------------------------

## Split "a__b__c" labels into columns `into` (missing pieces are NA, extra pieces dropped), as
## tidyr::separate(sep = "__") does.
.split_labels <- function(x, into, sep = "__") {
    parts <- strsplit(as.character(x), sep, fixed = TRUE)
    out <- lapply(seq_along(into), function(i) {
        vapply(parts, function(p) if (length(p) >= i) p[i] else NA_character_, "")
    })
    names(out) <- into
    as.data.frame(out, stringsAsFactors = FALSE)
}

## Row-bind data.frames whose columns may differ (missing columns filled with NA), as dplyr::bind_rows().
.bind_rows <- function(x) {
    x <- Filter(Negate(is.null), x)
    if (!length(x)) return(data.frame())
    x <- lapply(x, function(d) as.data.frame(as.list(d), check.names = FALSE, stringsAsFactors = FALSE))
    cols <- unique(unlist(lapply(x, names)))
    x <- lapply(x, function(d) {
        for (cl in setdiff(cols, names(d))) d[[cl]] <- NA
        d[cols]
    })
    out <- do.call(rbind, x)
    rownames(out) <- NULL
    out
}

## ---- deprecation -------------------------------------------------------------------------------

## A warning that `what` is deprecated in favour of `with`, as lifecycle::deprecate_warn() gave it.
.deprecate_warn <- function(when, what, with = NULL) {
    warning("`", what, "` was deprecated in spicyR ", when, ".",
            if (!is.null(with)) paste0(" Please use `", with, "` instead.") else "",
            call. = FALSE)
}

#' A format SummarizedExperiment and data.frame objects into a canonical form.
#'
#' @importFrom methods is
#' @importFrom S4Vectors as.data.frame
#' @noRd
.format_data <- function(
    cells, imageIDCol, cellTypeCol, spatialCoordCols, verbose = FALSE) {
    if (is.data.frame(cells)) {
        # pass
    } else if (.is_class(cells, "SpatialExperiment")) {
        .need("SpatialExperiment", "to use a SpatialExperiment")
        cd <- .col_data(cells)
        cd <- cd[, setdiff(colnames(cd), c("x", "y")), drop = FALSE]
        cells <- cbind(cd, data.frame(SpatialExperiment::spatialCoords(cells)))
    } else if (.is_class(cells, "SummarizedExperiment")) {
        cells <- .col_data(cells)
    } else {
        temp <- tryCatch(as.data.frame(cells), error = function(e) NULL)
        if (is.null(temp)) {
            stop("`cells` is an unsupported class: ", paste(class(cells), collapse = "/"), ". ",
                 "data.frame (or coercible), SingleCellExperiment and SpatialExperiment are currently supported.",
                 call. = FALSE)
        }
        cells <- temp
    }

    cols <- c(
        imageIDCol = imageIDCol,
        cellTypeCol = cellTypeCol,
        spatialCoordCols_x = spatialCoordCols[1],
        spatialCoordCols_y = spatialCoordCols[2]
    )
    for (i in seq_along(cols)) {
        if (!cols[[i]] %in% colnames(cells)) {
            stop("Specified `", names(cols)[i], "` (", cols[[i]], ") is not in `cells`. ",
                 "names(cells): ", paste(names(cells), collapse = ", "), call. = FALSE)
        }
    }

    needed <- data.frame(
        imageID = cells[[imageIDCol]],
        cellType = cells[[cellTypeCol]],
        x = cells[[spatialCoordCols[1]]],
        y = cells[[spatialCoordCols[2]]],
        stringsAsFactors = FALSE
    )
    keep <- setdiff(colnames(cells), c("imageID", "cellType", "x", "y"))
    cells <- cbind(as.data.frame(cells)[, keep, drop = FALSE], needed)

    # imageID and cellType as factors, levels in order of first appearance
    if (!is.factor(cells$imageID)) {
        cells$imageID <- factor(cells$imageID, levels = unique(cells$imageID))
    }
    if (!is.factor(cells$cellType)) {
        cells$cellType <- factor(cells$cellType, levels = unique(cells$cellType))
    }

    # Order the cells by imageID. This is important for when the data is split downstream
    cells <- cells[order(cells$imageID), , drop = FALSE]

    # create cellID if it does not exist
    if (is.null(cells$cellID)) {
        if (verbose) message("No column called cellID. Creating one.")
        cells$cellID <- paste0("cell", "_", seq_len(nrow(cells)))
    }

    # create imageCellID (row number within image) if it does not exist
    if (is.null(cells$imageCellID)) {
        if (verbose) message("No column called imageCellID. Creating one.")
        within <- stats::ave(seq_len(nrow(cells)), cells$imageID, FUN = seq_along)
        cells$imageCellID <- paste0(cells$imageID, "_", within)
        rownames(cells) <- NULL
    }

    cells
}

getCellSummary <- function(
    data,
    imageID = NULL,
    bind = TRUE) {
    if (!is.null(imageID)) data <- data[which(data$imageID == imageID), , drop = FALSE]
    data <- as.data.frame(data)[, c("imageID", "cellID", "imageCellID", "x", "y", "cellType"), drop = FALSE]
    data <- S4Vectors::DataFrame(data)
    if (bind) data else S4Vectors::split(data, data$imageID)
}


getImageID <- function(x, imageID = NULL) {
    if (!is.null(imageID)) x <- x[which(x$imageID == imageID), , drop = FALSE]
    x$imageID
}

getCellType <- function(x, imageID = NULL) {
    if (!is.null(imageID)) x <- x[which(x$imageID == imageID), , drop = FALSE]
    x$cellType
}

getImagePheno <- function(x,
                          imageID = NULL,
                          bind = TRUE,
                          expand = FALSE) {
    x <- x[!duplicated(x$imageID), ]
    x$imageID <- factor(x$imageID, levels = unique(x$imageID))
    x
}

#' A function to handle validity of argumemts/check for deprecated arguments
#'
#' Deprecated arguments are copied to their new names in the calling function's frame (`env`).
#' `cores`, `nCores` and `BPPARAM` may be a number or a BiocParallel param; the caller converts them
#' to a worker count with .n_workers().
#'
#' @importFrom methods is
#' @noRd
argumentChecks = function(function_name, user_vals, env = parent.frame()) {

  warned <- character()
  warn_once <- function(old_arg, new_arg) {
    if (!old_arg %in% warned) {
      .deprecate_warn("1.18.0", paste0(function_name, "(", old_arg, ")"), paste0(function_name, "(", new_arg, ")"))
      warned <<- c(warned, old_arg)
    }
  }

  # handle deprecated arguments
  handle_deprecated = function(old_arg, new_arg, user_vals) {
    if (old_arg %in% names(user_vals) && !(new_arg %in% names(user_vals))) {
      warn_once(old_arg, new_arg)
      assign(new_arg, user_vals[[old_arg]], envir = env)
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

  is_param <- function(x) .is_class(x, "BiocParallelParam")

  # validity checks for cores/nCores/BPPARAM
  if ("nCores" %in% names(user_vals)) {
    if (!is.numeric(user_vals$nCores) && !is_param(user_vals$nCores)) {
      stop("'nCores' must be either a numeric value, or a MulticoreParam or SerialParam object.\n")
    }
  } else if ("BPPARAM" %in% names(user_vals)) {
    if (!is_param(user_vals$BPPARAM)) {
      stop("'BBPARAM' must be a MulticoreParam or SerialParam object.")
    }
  } else if ("cores" %in% names(user_vals)) {
    if (!is.numeric(user_vals$cores) && !is_param(user_vals$cores)) {
      stop("'cores'  must be either a numeric value, or a MulticoreParam or SerialParam object.\n")
    }
  }
}
