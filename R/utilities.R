
# a function to convert SegmentedCells to
# to a SingleCellExperiment
.as_SingleCellExperiment <- function(data) {
    if  (!is(data, "SegmentedCells")) {
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

isKonditional <- function(konditionalResult){
    
    colNames = c(
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
