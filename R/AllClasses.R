#' @export
#' @rdname segmentedCells
#' @importClassesFrom S4Vectors DataFrame
setClass("segmentedCells", contains = "DataFrame")

#' @export
#' @rdname spicy-class
setClass("spicy", contains = "list")
