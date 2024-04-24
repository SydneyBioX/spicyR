# a function to convert SegmentedCells to
# to a SingleCellExperiment
.as_SingleCellExperiment <- function(data) {
    if (!is(data, "SegmentedCells")) {
        stop("This function only does SegmentedCells -> SingleCellExperiment")
    }

    colData <- as.data.frame(cellSummary(data)) |>
        dplyr::left_join(
            as.data.frame(imagePheno(data))
        )
    sce <- SingleCellExperiment::SingleCellExperiment(
        colData = colData,
    )
}

isKonditional <- function(konditionalResult) {
    colNames <- c(
        "imageID",
        "test",
        "original",
        "konditional",
        "r",
        "weightQuantile",
        "inhom",
        "edge",
        "includeZeroCells",
        "window",
        "window.length"
    )

    return(all(colNames %in% names(konditionalResult)))
}

#' A format SummarizedExperiment and data.frame objects into a canonical form.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr %>%
#' @importFrom SpatialExperiment spatialCoords
.format_data <- function(
    cells, imageIDCol, cellTypeCol, spatialCoordCols, verbose) {
    if (is(cells, "data.frame")) {
        # pass
    } else if (is(cells, "SpatialExperiment")) {
        spatialCoordCols <- names(spatialCoords(cells))
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

    cells <- cells %>%
        dplyr::rename(
            imageID = !!imageIDCol,
            cellType = !!cellTypeCol,
            x = spatialCoordCols[1],
            y = spatialCoordCols[2]
        )
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
