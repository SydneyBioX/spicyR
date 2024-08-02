# check if alternativeResults is konditional
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
#' @importFrom methods is
#' @importFrom dplyr select any_of rename mutate group_by ungroup
#' @importFrom cli cli_abort cli_inform
#' @noRd
.format_data <- function(
    cells, imageIDCol, cellTypeCol, spatialCoordCols, verbose = FALSE) {
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

getCellSummary <- function(
    data,
    imageID = NULL,
    bind = TRUE) {
    data %>%
        dplyr::filter(
            if (!is.null(!!imageID)) imageID == !!imageID else TRUE
        ) %>%
        dplyr::select(imageID, cellID, imageCellID, x, y, cellType) %>%
        dplyr::mutate(imageID = factor(imageID, levels = unique(imageID))) %>%
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
