#' @export
#' @rdname SegmentedCells
#' @importClassesFrom S4Vectors DataFrame
setClass("SegmentedCells", contains = "DataFrame")

#' @export
#' @aliases SpicyResults,list,ANY-method
#' @rdname spicy
setClass("SpicyResults", contains = "list")
