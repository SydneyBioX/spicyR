#' Box plot of one pair, with a point per image
#'
#' Shows the per-image value of one pair in each condition: a box plot with the images as points behind it.
#' For the cell method the value is the excess (extra `to` cells per `from` cell beyond chance), and each point
#' is sized by how much the image contributes to the test (its weight in the frailty model, relative to the
#' average image of its condition). With `interactive = TRUE` the plot is a plotly widget: hover over a point to
#' see its image, patient, excess and weight, which helps to find images worth looking at with [plotImage()].
#'
#' @param results The SpicyResults object returned by spicy() (either method).
#' @param from The `from` cell type (the centre).
#' @param to The `to` cell type (counted around each `from` cell).
#' @param rank Alternatively, the rank of the pair by p-value (1 is the most significant).
#' @param interactive Return an interactive plotly widget instead of a ggplot (needs the plotly package).
#'
#' @return A ggplot, or with `interactive = TRUE` a plotly htmlwidget. Images in which the pair was not
#'   tested (no `from` cells, for example) are not shown.
#'
#' @examples
#' data(spicyTest)
#' spicyBoxPlot(spicyTest, rank = 1)
#'
#' data("diabetesData")
#' res <- spicy(diabetesData, condition = "stage", subject = "case", r = 50,
#'              from = "Tc", to = c("Th", "beta"))
#' spicyBoxPlot(res, from = "Tc", to = "Th")
#'
#' @export
#' @importFrom ggplot2 ggplot aes .data
spicyBoxPlot <- function(results,
                         from = NULL,
                         to = NULL,
                         rank = NULL,
                         interactive = FALSE) {
  if (is.null(c(from, to, rank))) stop("Please specify either a pairwise relationship or rank")
  pVal <- results$p.value
  if (is.null(rank)) {
    if (length(c(from, to)) == 1) stop("Please specify both from and to parameters")
    pairName <- paste0(from, "__", to)
  } else {
    pVal <- pVal[order(pVal[, 2]), ]
    pairName <- rownames(pVal)[rank]
    from <- unlist(strsplit(pairName, split = "__"))[1]
    to <- unlist(strsplit(pairName, split = "__"))[2]
  }
  if (is.null(results$pairwiseAssoc[[pairName]])) stop("pair ", from, " -> ", to, " not found in the results.")

  cell <- identical(results$method, "cell")
  df <- data.frame(imageID = as.character(results$imageID),
                   value = results$pairwiseAssoc[[pairName]],
                   condition = results$condition)
  df$subject <- if (!is.null(results$subject)) as.character(results$subject) else df$imageID
  w <- results$imageWeights[[pairName]]
  df$weight <- if (cell && !is.null(w)) w else NA_real_
  df <- df[is.finite(df$value) & !is.na(df$condition), , drop = FALSE]
  # relative weight: 1 for an image of average weight in its condition
  df$relative <- stats::ave(df$weight, df$condition, FUN = function(z) z / mean(z, na.rm = TRUE))
  sized <- any(is.finite(df$relative))
  if (cell) {
    ylabel <- paste0("Extra ", to, " per ", from, "\n(beyond chance)")
    title <- paste0(to, " around ", from)
  } else {
    ylabel <- if (isTRUE(results$alternateResult)) "Alternate Result" else "L Function"
    title <- paste0("L-function values between ", from, " cells and ", to, " cells")
  }
  df$hover <- paste0("image: ", df$imageID,
                     if (!identical(df$subject, df$imageID)) paste0("<br>patient: ", df$subject) else "",
                     "<br>", if (cell) "excess: " else "value: ", signif(df$value, 3),
                     if (sized) paste0("<br>relative weight: ", signif(df$relative, 2)) else "")

  jitter <- ggplot2::position_jitter(width = 0.2, height = 0, seed = 1)
  pt <- if (sized) ggplot2::aes(size = .data$relative, text = .data$hover) else ggplot2::aes(text = .data$hover)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$condition, y = .data$value)) +
    suppressWarnings(ggplot2::geom_point(pt, position = jitter, colour = "grey45", alpha = 0.5, stroke = 0)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = .data$condition), outlier.shape = NA, alpha = 0.35, width = 0.55) +
    ggplot2::ggtitle(title) +
    ggplot2::labs(x = NULL, y = ylabel, fill = NULL, size = "Relative\nweight") +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_classic()
  if (sized) p <- p + ggplot2::scale_size_area(max_size = 4)
  if (cell) p <- p + ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey50")
  if (!interactive) return(p)
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("interactive = TRUE needs the plotly package: install.packages(\"plotly\")", call. = FALSE)
  plotly::ggplotly(p, tooltip = "text")
}
