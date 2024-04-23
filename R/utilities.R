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

#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr %>%
#' @importFrom SpatialExperiment spatialCoords
.format_data <- function(cells, imageIDCol, cellTypeCol, spatialCoordCols) {
    if (is(cells, "SpatialExperiment")) {
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
        stop(
            "`cells` is an unsupported subclass of SummarizedExperiment. ",
            "SingleCellExperiment and SpatialExperiment",
            " are currently supported."
        )
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
        cli::cli_inform(
            "No column called cellID. Creating one."
        )
        cells <- dplyr::mutate(
            cells,
            cellID = paste0("cell", "_", dplyr::row_number())
        )
    }

    # create imageCellID if it does not exist
    if (is.null(cells$imageCellID)) {
        cli::cli_inform(
            "No column called imageCellID. Creating one."
        )
        cells <- cells %>%
            dplyr::group_by(imageID) %>%
            dplyr::mutate(
                imageCellID = paste0(imageID, "_", dplyr::row_number())
            ) %>%
            dplyr::ungroup()
    }
    cells
}
