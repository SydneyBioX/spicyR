#' The evidence for pairs at each radius
#'
#' For a [spicy()] result with several radii, plots the p-value of each chosen pair at every radius, to show the
#' scale at which the evidence is strongest. The combined p-value over the radii is in [topPairs()].
#'
#' @param results The `SpicyResults` object returned by `spicy()` with several radii (`r`).
#' @param from,to The `from` and `to` cell types of the pairs to show: every pair between them.
#'
#' @return A ggplot.
#'
#' @examples
#' data("diabetesData")
#' res <- spicy(diabetesData, condition = "stage", subject = "case", r = c(25, 50, 75),
#'              from = c("Tc", "Th"), to = "beta")
#' radiusPlot(res, from = c("Tc", "Th"), to = "beta")
#'
#' @export
#' @concept plot
radiusPlot <- function(results, from, to) {
  tab <- results$radiusResults
  if (is.null(tab)) stop("radiusPlot() needs a result of spicy() with several radii.", call. = FALSE)
  tab <- .pickPairs(tab, from, to)
  tab$pair <- paste(tab$from, "\u2192", tab$to)
  what <- paste(rev(levels(results$condition)), collapse = " vs ")
  ggplot(tab, aes(.data$r, -log10(.data$p_value), colour = .data$pair)) +
    geom_line() +
    geom_point() +
    expand_limits(y = 0) +
    labs(x = "radius", y = paste0("-log10 p at each radius, ", what), colour = NULL) +
    theme_classic()
}
