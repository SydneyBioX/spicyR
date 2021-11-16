#' @export
#' @rdname SegmentedCells
#' @importClassesFrom S4Vectors DFrame
setClass("SegmentedCells", contains = "DFrame")

#' @export
#' @aliases SpicyResults,list,ANY-method
#' @rdname spicy
setClass("SpicyResults", contains = "list")
