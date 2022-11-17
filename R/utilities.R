
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