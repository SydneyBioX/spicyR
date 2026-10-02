#' Plots boxplot for a specified cell-cell relationship
#'
#' @param results The SpicyResults object returned by spicy() (either method).
#' @param from Cell type which you would like to compare to the to cell type.
#' @param to Cell type which you would like to compare to the from cell type.
#' @param rank Ranking of cell type in terms of p-value, the smaller the p-value
#'   the higher the rank.
#'
#' @return a ggplot2 boxplot. For \code{method = "cell"} results the y axis is the
#'   per-image excess (extra \code{to} cells per \code{from} cell); for
#'   \code{method = "image"} results it is the L function.
#'
#' @examples
#' data(spicyTest)
#' 
#' spicyBoxPlot(spicyTest,
#'              rank = 1)
#'
#' @export
#' @importFrom ggplot2 ggplot aes .data

spicyBoxPlot <- function(results,
                         from = NULL,
                         to = NULL,
                         rank = NULL) {
  
  if(is.null(c(from, to, rank))) {
    stop("Please specify either a pairwise relationship or rank")
  }
  
  pVal <- results$p.value
  
  if(is.null(rank)) {
    if(length(c(from, to)) == 1) {
      stop("Please specify both from and to parameters")
    }
    pairName <- paste0(from, "__", to)
  }
  
  if(!is.null(rank)) {
    pVal <- pVal[order(pVal[, 2]),]
    
    pairName <- rownames(pVal)[rank]
    from <- unlist(strsplit(pairName, split = "__"))[1]
    to <- unlist(strsplit(pairName, split = "__"))[2]
  }
  
  df <- data.frame(imageID = results$imageID, 
                   pairwiseAssoc = results$pairwiseAssoc[[pairName]],
                   condition = results$condition)
  
  if (identical(results$method, "cell")) {
    ylabel <- paste0("Excess (extra ", to, " cells per ", from, " cell)")
    title <- paste0(to, " cells around ", from, " cells")
  } else {
    ylabel <- if (isTRUE(results$alternateResult)) "Alternate Result" else "L Function"
    title <- paste0("L-function values between ", from, " cells and ", to, " cells")
  }
  
  ggplot2::ggplot(df, ggplot2::aes(x = .data$condition, y = .data$pairwiseAssoc, fill = .data$condition)) +
    ggplot2::geom_boxplot() +
    ggplot2::ggtitle(title) +
    ggplot2::xlab("Condition") + 
    ggplot2::ylab(ylabel) +
    ggplot2::theme_classic()
}
