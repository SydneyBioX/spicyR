#' @export
#' @rdname segmentedCells
#' @importClassesFrom S4Vectors DataFrame
setClass("segmentedCells", contains = "DataFrame")

#' @export
#' @aliases spicy,list,ANY-method
#' @rdname spicy
setClass("spicy", contains = "list")
