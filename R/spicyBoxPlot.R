#' Plots boxplot for a specified cell-cell relationship
#'
#' @param results Data frame obtained from spicy.
#' @param from Cell type which you would like to compare to the to cell type.
#' @param to Cell type which you would like to compare to the from cell type.
#' @param rank Ranking of cell type in terms of p-value, the smaller the p-value
#'   the higher the rank.
#'
#' @return a ggplot2 boxplot
#'
#' @examples
#' data(spicyTest)
#' 
#' spicyBoxPlot(spicyTest,
#'              rank = 1)
#'
#' @export
#' @importFrom ggplot2 ggplot scale_colour_gradient2 geom_point scale_shape_manual guides labs scale_color_manual theme_classic theme element_text aes guide_legend element_blank guide_colourbar

spicyBoxPlot <- function(results,
                         from = NULL,
                         to = NULL,
                         rank = NULL) {
  
  if(is.null(c(from, to, rank))) {
    stop("Please specify either a pairwse relationship or rank")
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
  }
  
  df <- data.frame(imageID = results$imageID, 
                   pairwiseAssoc = results$pairwiseAssoc[pairName],
                   condition = results$condition)
  
  if(results$alternateResult) {
    ylabel <- "Alternate Result"
  } else {
    ylabel <- "L Function"
  }
  
  ggplot2::ggplot(df, ggplot2::aes(x = condition, y = .data[[pairName]], fill = condition, label = imageID)) +
    ggplot2::geom_boxplot() +
    # ggplot2::geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
    ggplot2::ggtitle("Boxplot of Pairwise Assocations across Conditions") +
    ggplot2::xlab("Condition") + 
    ggplot2::ylab(ylabel) +
    ggplot2::theme_classic()
}
