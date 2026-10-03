#' Plots result of signifPlot.
#'
#' @param results A spicy results object
#' @param type Whether to make a bubble plot (\code{"bubble"}) or a heatmap (any other value). Note: For
#'   survival results a bubble plot will be used.
#' @param fdr TRUE if FDR correction is used.
#' @param breaks Vector of 3 numbers giving breaks used in legend. The first
#'     number is the minimum, the second is the maximum, the third is the
#'     number of breaks.
#' @param comparisonGroup A string specifying the name of the outcome group to compare with the base group.
#' @param colours Vector of colours to use to colour legend.
#' @param marksToPlot Vector of marks to include in plot.
#' @param cutoff significance threshold for circles in bubble plot.
#' @param contextColours Used for \code{\link[Statial]{Kontextual}} results. A named list specifying the colours for each context.
#'  By default the Tableau colour palette is used.
#' @param contextLabels Used for \code{\link[Statial]{Kontextual}} results. A named list to change the default labels for each context.
#'
#' @return a ggplot object
#'
#' @examples
#' data(spicyTest)
#'
#' p <- signifPlot(spicyTest, breaks = c(-3, 3, 0.5))
#' p
#' signifPlot(spicyTest, type = "heatmap")
#'
#' @export
#' @importFrom grDevices colorRampPalette
#' @importFrom stats p.adjust
#' @importFrom ggplot2
#'     ggplot scale_colour_gradient2 geom_point scale_shape_manual guides labs
#'     scale_color_manual theme_classic theme element_text aes guide_legend
#'     element_blank guide_colourbar waiver .data
#' @importFrom stats setNames
signifPlot <- function(results,
                       fdr = FALSE,
                       type = "bubble",
                       breaks = NULL,
                       comparisonGroup = NULL,
                       colours = c("#4575B4", "white", "#D73027"),
                       marksToPlot = NULL,
                       cutoff = 0.05,
                       contextColours = NULL,
                       contextLabels = waiver()) {

  if (is.null(comparisonGroup)) {
    coef <- 2
  } else {
    coef <- which(levels(results$condition) == comparisonGroup)
  }
  # every cell type in the tested pairs (from and to)
  marks <- unique(c(
    as.character(results$comparisons$from),
    as.character(results$comparisons$to)
  ))


  if (is.null(marksToPlot)) marksToPlot <- marks


  # the cell method's survival results (spicy() or Statial::kontextualTest()): the log hazard ratio per SD of the excess
  if (identical(results$method, "cell") && !is.null(results$survivalOutcome) && !"survivalResults" %in% names(results)) {
    tab <- results$cellResults
    results$survivalResults <- data.frame(test = rownames(tab), coef = tab$log_hr_sd, p.value = tab$p_value)
  }

  if("survivalResults" %in% names(results)) {
    return(
      survBubble(result = results,
                 fdr = fdr,
                 cutoff = cutoff,
                 colourGradient = colours,
                 marksToPlot = marksToPlot,
                 contextColours = contextColours,
                 contextLabels = contextLabels)
    )
  }


  if (type == "bubble") {
    return(
      bubblePlot(
        results, fdr, breaks, coef,
        colours = colours, cutoff = cutoff, marksToPlot = marksToPlot,
        contextColours = contextColours, contextLabels = contextLabels
      )
    )
  }

  heatmapPlot(results, fdr, breaks, coef, colours = colours, marksToPlot = marksToPlot)
}

# Heatmap of signed log10 p-values (positive when the coefficient is positive), rows `from` and columns
# `to`, as spicyR 1.x drew with pheatmap: colours from `colours` over `breaks`, values beyond the breaks
# shown at the end colours.
heatmapPlot <- function(results, fdr, breaks, coef, colours, marksToPlot) {
  if (is.null(breaks)) breaks <- c(-3, 3, 0.5)
  breaks <- seq(from = breaks[1], to = breaks[2], by = breaks[3])
  pal <- grDevices::colorRampPalette(colours)(length(breaks))

  pVal <- results$p.value[, coef]

  if (min(pVal, na.rm = TRUE) == 0) {
    pVal[pVal == 0] <-
      pVal[pVal == 0] + 10^floor(log10(min(pVal[pVal > 0], na.rm = TRUE)))
  }

  if (fdr) {
    pVal <- stats::p.adjust(pVal, method = "fdr")
  }

  isGreater <- results$coefficient[, coef] > 0

  pVal <- log10(pVal)

  pVal[which(isGreater)] <- abs(pVal[which(isGreater)])

  from <- as.character(results$comparisons$from)
  to <- as.character(results$comparisons$to)
  df <- data.frame(from = from, to = to, value = pVal)
  df <- df[df$from %in% marksToPlot & df$to %in% marksToPlot, , drop = FALSE]
  rowMarks <- intersect(marksToPlot, unique(df$from))
  colMarks <- intersect(marksToPlot, unique(df$to))
  df$from <- factor(df$from, levels = rev(rowMarks))
  df$to <- factor(df$to, levels = colMarks)

  limits <- range(breaks)
  squish <- function(x, range = limits, only.finite = TRUE) pmin(pmax(x, range[1]), range[2])

  ggplot2::ggplot(df, ggplot2::aes(x = .data$to, y = .data$from, fill = .data$value)) +
    ggplot2::geom_tile(colour = "grey60") +
    ggplot2::scale_fill_gradientn(
      colours = pal, values = (breaks - limits[1]) / diff(limits), limits = limits,
      breaks = breaks, oob = squish, na.value = "grey90"
    ) +
    ggplot2::scale_x_discrete(position = "bottom", guide = ggplot2::guide_axis(angle = 90)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL, fill = if (fdr) "signed\nlog10 FDR" else "signed\nlog10 p") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

# The Tableau 10 and Tableau 20 palettes (as in ggthemes), the default context colours.
.tableau10 <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948",
                "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC")
.tableau20 <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", "#8CD17D",
                "#B6992D", "#F1CE63", "#499894", "#86BCB6", "#E15759", "#FF9D9A",
                "#79706E", "#BAB0AC", "#D37295", "#FABFD2", "#B07AA1", "#D4A6C8",
                "#9D7660", "#D7B5A6")

# Default colours for the contexts of a Kontextual result.
.context_colours <- function(parents) {
  palette <- if (length(parents) > 10) .tableau20 else .tableau10
  stats::setNames(rep_len(palette, length(parents)), parents)
}

# Facet labeller for contexts: `contextLabels` (a named vector or list) renames contexts.
.context_labeller <- function(contextLabels) {
  if (inherits(contextLabels, "waiver") || is.null(contextLabels)) return(ggplot2::label_value)
  lab <- unlist(contextLabels)
  ggplot2::as_labeller(function(x) {
    out <- as.character(x)
    hit <- out %in% names(lab)
    out[hit] <- as.character(lab[out[hit]])
    out
  })
}

# Function to create x and y points for a half circle (left or right)
half_circle_coords = function(shape = "left", num_points = 100) {
  # Generate points from -pi/2 to pi/2 for vertical half circles
  t = seq(-pi / 2, pi / 2, length.out = num_points)
  if (shape == "left") {
    x = 0.5 - 0.5 * cos(t) # Shift to the left half
  } else {
    x = 0.5 + 0.5 * cos(t)  # Shift to the right half
  }
  y = 0.5 + 0.5 * sin(t) # Vertical component
  list(x = c(x,x[1])-mean(c(x,x[1]))+0.5, y = c(y,y[1]))
}

# Custom draw_key function to draw a left (base group) or right half circle in the legend
draw_key_half_circle = function(data, params, shape) {
  coords <- half_circle_coords(shape = if (data$shape == 16) "left" else "right")
  grid::grobTree(
    grid::polygonGrob(
      x = coords$x, y = coords$y,
      gp = grid::gpar(fill = "black", col = "black")
    )
  )
}

# Polygons for discs (or half discs) centred at (x0, y0) with radii r, angles measured clockwise from
# 12 o'clock as ggforce::geom_arc_bar() measures them: from = 0, to = pi is the right half. Half discs
# include the centre (pie wedges).
.disc_polygons <- function(x0, y0, r, from = 0, to = 2 * pi, n = 60) {
  if (!length(x0)) return(data.frame(.id = integer(), .x = numeric(), .y = numeric(), .row = integer()))
  a <- seq(from, to, length.out = n)
  full <- isTRUE(all.equal(to - from, 2 * pi))
  ang <- if (full) a else c(NA, a)
  m <- length(ang)
  row <- rep(seq_along(x0), each = m)
  aa <- rep(ang, times = length(x0))
  rr <- r[row]
  data.frame(
    .id = row,
    .x = x0[row] + ifelse(is.na(aa), 0, rr * sin(aa)),
    .y = y0[row] + ifelse(is.na(aa), 0, rr * cos(aa)),
    .row = row
  )
}

bubblePlot <- function(test,
                       fdr,
                       breaks,
                       coef,
                       colours = c("blue", "white", "red"),
                       cutoff = 0.05,
                       marksToPlot,
                       contextColours = NULL,
                       contextLabels = waiver()) {



  if (is.null(test$alternateResult)) {
    test$alternateResult <- FALSE
  }
  isCell <- identical(test$method, "cell")

  if (test$alternateResult || isCell) {
    # alternate results and the cell method's excess are plotted on their own scale
    groupA <- test$coefficient[, 1]
    groupB <- (test$coefficient[, 1] + test$coefficient[, coef])
  } else {
    groupA <- test$coefficient[, 1] * sqrt(pi) * 2 / sqrt(10) / 100
    groupB <- (
      test$coefficient[, 1] + test$coefficient[, coef]
    ) * sqrt(pi) * 2 / sqrt(10) / 100

  }


  cellTypeA <- factor(test$comparisons$from)
  cellTypeB <- factor(test$comparisons$to)


  pvalue = test$p.value[, coef]
  sig <- pvalue < cutoff
  sigLab <- paste0("p-value < ", cutoff)


  if (fdr) {
    pvalue = p.adjust(test$p.value[, coef], "fdr")
    sig <- pvalue < cutoff
    sigLab <- paste0("BH-adjusted p-value < ", cutoff)
  }

  size <- -log10(pvalue)


  df <- data.frame(
    cellTypeA, cellTypeB, groupA, groupB, size,
    stat = test$statistic[, coef], pvalue = pvalue,
    sig = factor(sig, levels= c("FALSE", "TRUE"))
  )
  rownames(df) <- rownames(test$statistic)

  isKontextual <- isTRUE(test$isKontextual)
  if (isKontextual) {
    df$parent = test$comparisons$parent
  }

  df <- df[df$cellTypeA %in% marksToPlot & df$cellTypeB %in% marksToPlot, ]

  df$cellTypeA <- droplevels(df$cellTypeA)
  df$cellTypeB <- droplevels(df$cellTypeB)

  df.shape <- data.frame(
    cellTypeA = c(NA, NA), cellTypeB = c(NA, NA), size = c(1, 1),
    condition = c(
      levels(test$condition)[1], levels(test$condition)[coef]
    )
  )

  if(is.null(breaks)) {
    groupAB <- c(groupA, groupB)

    limits <- c(min(groupAB, na.rm = TRUE), max(groupAB, na.rm = TRUE)) |>
      round(1)
    breaks <- seq(from = limits[1], to = limits[2], by = diff(limits) / 5)

  } else {
    limits <- c(breaks[1], breaks[2])
    breaks <- seq(from = breaks[1], to = breaks[2], by = breaks[3])
  }

  midpoint <- 0

  if(test$alternateResult && !isKontextual){
    midpoint <- (breaks[1] + breaks[length(breaks)])/2
  }


  df$groupA <- pmax(pmin(df$groupA, limits[2]), limits[1])
  df$groupB <- pmax(pmin(df$groupB, limits[2]), limits[1])

  labels <- round(breaks, 3)
  labels[1] <- "avoidance"
  labels[length(labels)] <- "attraction"


  if (isKontextual) {
    # positions within each context panel
    df$cellTypeB_numeric <- stats::ave(
      as.integer(df$cellTypeB), df$parent,
      FUN = function(v) as.integer(factor(v, levels = sort(unique(v))))
    )
    df$cellTypeB_id <- factor(paste(df$parent, df$cellTypeB, sep = "_"))
  } else {
    df$cellTypeB_numeric <- as.numeric(df$cellTypeB)
    df$cellTypeB_id <- df$cellTypeB
  }
  df$cellTypeA_numeric <- as.numeric(df$cellTypeA)
  df$radius <- pmax(df$size / max(df$size, na.rm = TRUE) / 2, 0.15)

  # half discs: the comparison group on the right, the base group on the left
  drawn <- df[is.finite(df$radius), , drop = FALSE]
  halfDiscs <- function(value, from, to) {
    poly <- .disc_polygons(drawn$cellTypeB_numeric, drawn$cellTypeA_numeric, drawn$radius, from, to)
    poly$.fill <- value[poly$.row]
    if (isKontextual) poly$parent <- drawn$parent[poly$.row]
    poly
  }
  right <- halfDiscs(drawn$groupB, 0, pi)
  left <- halfDiscs(drawn$groupA, pi, 2 * pi)

  df.shape$condition <- factor(df.shape$condition, levels = levels(test$condition))

  xLabels <- stats::setNames(as.character(df$cellTypeB), as.character(df$cellTypeB_id))

  plot = ggplot2::ggplot(df, ggplot2::aes(x = .data$cellTypeB_id, y = .data$cellTypeA)) +
    ggplot2::scale_fill_gradient2(
      low = colours[1], mid = colours[2], high = colours[3],
      midpoint = midpoint, breaks = breaks, labels = labels, limits = limits
    ) +
    ggplot2::geom_point(ggplot2::aes(col = sigLab), size = -1) +
    ggplot2::geom_point(ggplot2::aes(size = .data$size), x = 100000, y = 10000000) +
    ggplot2::geom_polygon(
      data = right, ggplot2::aes(x = .data$.x, y = .data$.y, group = .data$.id, fill = .data$.fill),
      colour = NA, inherit.aes = FALSE
    ) +
    ggplot2::geom_polygon(
      data = left, ggplot2::aes(x = .data$.x, y = .data$.y, group = .data$.id, fill = .data$.fill),
      colour = NA, inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = df.shape, ggplot2::aes(shape = .data$condition), x = 10000, y = 10000,
      key_glyph = draw_key_half_circle, inherit.aes = FALSE
    ) +
    ggplot2::scale_x_discrete(labels = function(b) unname(xLabels[as.character(b)]),
                              guide = ggplot2::guide_axis(angle = 45)) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "to (centre)", y = "from (counted)", size = if (fdr) "-log10 adjusted p-value" else "-log10 p-value",
      colour = NULL, fill = "Localisation", shape = "Condition"
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(order = 3,override.aes = list(size = 5)),
      size = ggplot2::guide_legend(order = 2),
      colour = ggplot2::guide_legend(
        order = 1, override.aes = list(size = 5, shape = 1, col = "black")
      )
    )

  # Plots black circle outlines, only if there are significant results.
  sigRows <- drawn[drawn$sig == "TRUE" & !is.na(drawn$sig), , drop = FALSE]
  if (nrow(sigRows) > 0) {
    circles <- .disc_polygons(sigRows$cellTypeB_numeric, sigRows$cellTypeA_numeric, sigRows$radius)
    if (isKontextual) circles$parent <- sigRows$parent[circles$.row]
    plot = plot +
      ggplot2::geom_polygon(
        data = circles, ggplot2::aes(x = .data$.x, y = .data$.y, group = .data$.id),
        fill = NA, colour = "black", inherit.aes = FALSE
      )
  }

  # Adds context panels if using Kontextual results: one panel per context, labelled in its strip,
  # with a band in the context's colour above the panel.
  if (isKontextual) {
    parents <- sort(unique(as.character(df$parent)))
    if (is.null(contextColours)) contextColours <- .context_colours(parents)
    contextColours <- unlist(contextColours)
    if (is.null(names(contextColours))) names(contextColours) <- parents[seq_along(contextColours)]
    top <- nlevels(df$cellTypeA) + 0.5
    # the fill scale is taken by the localisation, so each band is its own layer with a fixed colour
    for (p in intersect(parents, names(contextColours))) {
      plot <- plot + ggplot2::geom_rect(
        data = data.frame(parent = p, xmin = -Inf, xmax = Inf, ymin = top + 0.05, ymax = top + 0.3),
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
        fill = contextColours[[p]], inherit.aes = FALSE
      )
    }
    plot = plot +
      ggplot2::facet_grid(~parent, scales = "free_x", space = "free_x",
                          labeller = .context_labeller(contextLabels)) +
      ggplot2::theme(panel.spacing = ggplot2::unit(0.4, "lines"))
  }

  return(plot)

}

#' Plots survival results from spicy.
#'
#' @param result A spicyResults object that contains survival results.
#' @param fdr TRUE if FDR correction is used.
#' @param cutoff Significance threshold for circles in bubble plot.
#' @param colourGradient A vector of colours, used to define the low, medium, and high values for the colour scale.
#' @param marksToPlot Vector of marks to include in bubble plot.
#' @param contextColours Used for \code{\link[Statial]{Kontextual}} results. A named list specifying the colours for each context.
#'  By default the Tableau colour palette is used.
#' @param contextLabels Used for \code{\link[Statial]{Kontextual}} results. A named list to change the default labels for each context.
#'
#'
#' @return A ggplot object.
#'
#' @noRd
survBubble = function(result,
                      fdr = FALSE,
                      cutoff = 0.05,
                      colourGradient = c("#4575B4", "white", "#D73027"),
                      marksToPlot = NULL,
                      contextColours = NULL,
                      contextLabels = waiver()){


  if(!"survivalResults" %in% names(result)) {
    stop("Survival results are missing, please run spicy with survival outcomes.")
  }

  survivalResults = as.data.frame(result$survivalResults)
  isKontextual <- isTRUE(result$isKontextual)

  if (isKontextual) {
    plotData <- cbind(.split_labels(survivalResults$test, c("from", "to", "parent")),
                      survivalResults[setdiff(names(survivalResults), "test")])
    plotData <- plotData[order(plotData$parent, plotData$to, plotData$from, method = "radix"), , drop = FALSE]
    plotData$toParent <- paste(plotData$to, plotData$parent, sep = "__")
  } else {
    plotData <- cbind(.split_labels(survivalResults$test, c("from", "to")),
                      survivalResults[setdiff(names(survivalResults), "test")])
    plotData <- plotData[order(plotData$to, plotData$from, method = "radix"), , drop = FALSE]
  }

  if(!is.null(marksToPlot)) {
    plotData <- plotData[plotData$to %in% marksToPlot & plotData$from %in% marksToPlot, , drop = FALSE]
  }

  sigLab <- paste0("p-value < ", cutoff)

  if(fdr){
    plotData$p.value = p.adjust(plotData$p.value, "fdr")
    sigLab <- paste0("BH-adjusted p-value < ", cutoff)
  }

  plotData$sig <- plotData$p.value < cutoff
  plotData$logP <- -log10(plotData$p.value)
  plotData$size <- plotData$logP / max(plotData$logP, na.rm = TRUE)
  plotData$from <- factor(plotData$from)
  plotData$to <- factor(plotData$to)

  plot = ggplot2::ggplot(plotData, ggplot2::aes(x = .data$to, y = .data$from)) +
    ggplot2::geom_point(ggplot2::aes(size = pmax(.data$logP / 2, 0.15), colour = .data$coef)) +
    ggplot2::geom_point(data = plotData[which(plotData$sig), , drop = FALSE],
                        ggplot2::aes(size = pmax(.data$logP / 2, 0.15)),
                        shape = 21) +
    ggplot2::geom_point(ggplot2::aes(shape = sigLab), size = -1) +
    ggplot2::scale_colour_gradient2(low = colourGradient[[1]],
                                    mid = colourGradient[[2]],
                                    high = colourGradient[[3]],
                                    midpoint = 0) +
    ggplot2::scale_size(range = c(2, 6)) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45)) +
    ggplot2::labs(colour = "CoxPH \ncoefficient",
                  size = if (fdr) "-log10 adjusted p-value" else "-log10 p-value",
                  shape = NULL,
                  x = NULL,
                  y = NULL) +
    ggplot2::guides(shape = ggplot2::guide_legend(order = 1, override.aes = list(size=5, shape = 1, col = "black")),
                    colour = ggplot2::guide_colourbar(order = 2),
                    size = "none") +
    ggplot2::theme_classic()


  if (isKontextual) {
    # Adding context information to the plot: one panel per context, with a band in its colour
    parents <- sort(unique(plotData$parent))
    if (is.null(contextColours)) contextColours <- .context_colours(parents)
    contextColours <- unlist(contextColours)
    if (is.null(names(contextColours))) names(contextColours) <- parents[seq_along(contextColours)]
    top <- nlevels(plotData$from) + 0.5
    band <- data.frame(parent = parents, xmin = -Inf, xmax = Inf, ymin = top + 0.05, ymax = top + 0.3)
    fillLabels <- if (inherits(contextLabels, "waiver")) ggplot2::waiver() else {
      lab <- unlist(contextLabels)
      out <- parents
      out[parents %in% names(lab)] <- as.character(lab[parents[parents %in% names(lab)]])
      out
    }
    plot = plot +
      ggplot2::geom_rect(data = band,
                         ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin,
                                      ymax = .data$ymax, fill = .data$parent),
                         inherit.aes = FALSE) +
      ggplot2::facet_grid(~parent, scales = "free_x", space = "free_x",
                          labeller = .context_labeller(contextLabels)) +
      ggplot2::scale_fill_manual(values = contextColours, labels = fillLabels) +
      ggplot2::labs(fill = "Context") +
      ggplot2::guides(fill = ggplot2::guide_legend(order = 2)) +
      ggplot2::theme(panel.spacing = ggplot2::unit(0.4, "lines"))
  }


  return(plot)

}
