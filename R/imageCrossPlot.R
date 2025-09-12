#' Pairwise cross-plot of spatial associations
#'
#' Builds a pairwise scatterplot of marker-marker associations for a single image
#' using the output of `getPairwise()`. 
#'
#' @param result A matrix or data frame produced by `getPairwise()` where
#'   rownames are image IDs and colnames are concatenated marker pairs in the
#'   form `"from__to"`. Must contain at least one row with name matching `image`.
#' @param image Character scalar. Which image (row in `result`) to plot. If
#'   `NULL`, defaults to the first rowname of `result`.
#' @param colourGradient Character vector of length 3 giving the low, mid, and
#'   high colours for a diverging palette used by `scale_colour_gradient2()`.
#'   Default is `c("#4575B4", "white", "#D73027")`.
#' @param marksToPlot Optional character vector of marker names. If supplied, the
#'   plot is restricted to rows and columns where both `from` and `to` are in
#'   this set.
#' @param limits Numeric length-2 vector giving the lower and upper caps applied
#'   to association values for colour mapping. Use `NULL` to avoid
#'   clamping. Point sizes are clamped to `c(0, max(abs(limits)))`.
#'
#' @details
#' The function expects `colnames(result)` to contain `"__"` separating marker
#' names; otherwise an error is thrown. Values are transformed into two aesthetics:
#' \itemize{
#'   \item \code{value}: the (possibly clamped) signed association used for colour.
#'   \item \code{size}: the (possibly clamped) absolute association used for point size.
#' }
#' The colourbar shows three ticks at \code{lower}, \code{0}, and \code{upper}
#' with labels \code{"<lower"}, \code{"0"}, and \code{">upper"}.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' # Minimal example (toy data)
#' set.seed(1)
#' mks <- c("A","B","C")
#' pairs <- as.vector(outer(mks, mks, paste, sep="__"))
#' res <- matrix(rnorm(length(pairs)*2), nrow = 2,
#'               dimnames = list(c("img1","img2"), pairs))
#'
#' p <- imageCrossPlot(res, image = "img1",
#'                        colourGradient = c("#4575B4","white","#D73027"),
#'                        marksToPlot = c("A","B","C"),
#'                        limits = c(-3, 3))
#' print(p)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_gradient2 labs theme_classic guide_axis guides
#' @importFrom tidyr separate
#' @importFrom dplyr arrange filter
#' @importFrom scales squish
#' @seealso \code{\link{getPairwise}}
#' @export
imageCrossPlot = function(result,
                          image = NULL,
                          colourGradient = c("#4575B4", "white", "#D73027"),
                          marksToPlot = NULL,
                          limits = NULL){
  
  
  if(!any(grep("__", colnames(result)))) {
    stop("It doesn't look like this is output from `getPairwise()`")
  }
  
  if(is.null(image))image = rownames(result)[1]
  
  if(!image%in%rownames(result)) stop(paste(image, "is not in rownames of `result`"))
  
  data <- data.frame(test = colnames(result), 
                     value = as.numeric(result[image,]), 
                     size = abs(as.numeric(result[image,])))
  
  plotData = data |>
    tidyr::separate(test,
                    into = c("from", "to"),
                    sep = "__") |>
    dplyr::arrange(to, from)
  
  
  if(!is.null(marksToPlot)) {
    plotData = plotData |> 
      filter(to %in% marksToPlot) |> 
      filter(from %in% marksToPlot)
  }
  
  if(is.null(limits)) limits <- signif(range(data$value),2)
  
  val_limits  <- limits                   # e.g. c(-10, 10)
  size_limits <- c(0, max(abs(limits)))   # e.g. c(0, 10)
  
  # colour scale
  val_breaks <- c(val_limits[1], 0, val_limits[2])
  val_labels <- c(
    paste0("< ", val_limits[1]),
    "0",
    paste0("> ", val_limits[2])
  )
  
  # size scale
  size_breaks <- size_limits
  size_labels <- c(
    paste0("< ", size_limits[1]),
    paste0("> ", size_limits[2])
  )
  
  plotData$size_clamped <- pmin(pmax(plotData$size, size_limits[1]), size_limits[2])
  
  plot <- ggplot(plotData, aes(x = to, y = from)) +
    geom_point(aes(size = size_clamped, colour = value)) +
    scale_colour_gradient2(
      low = colourGradient[[1]], mid = colourGradient[[2]], high = colourGradient[[3]],
      midpoint = 0,
      limits = val_limits,
      oob = scales::squish,
      breaks = val_breaks,
      labels = val_labels,
      name = "Spatial association"
    ) +
    guides(size = "none") + 
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(x = NULL, y = NULL) +
    theme_classic()
  
  return(plot)
  
}