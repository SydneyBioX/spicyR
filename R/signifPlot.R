#' Plots result of signifPlot.
#'
#' @param results Data frame obtained from spicy.
#' @param type Where to make a bubble plot or heatmap.
#' @param fdr TRUE if FDR correction is used.
#' @param breaks Vector of 3 numbers giving breaks used in pheatmap. The first
#'     number is the minimum, the second is the maximum, the third is the
#'     number of breaks.
#' @param comparisonGroup A string specifying the name of the outcome group to compare with the base group.
#' @param colours Vector of colours to use in pheatmap.
#' @param marksToPlot Vector of marks to include in pheatmap.
#' @param cutoff significance threshold for circles in bubble plot.
#'
#' @return a pheatmap object
#'
#' @examples
#' data(spicyTest)
#'
#' p <- signifPlot(spicyTest, breaks = c(-3, 3, 0.5))
#' # plot includes unicode characters, do not use default pdf device
#' ggplot2::ggsave(p, filename = tempfile(), device = cairo_pdf)
#'
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom stats p.adjust
#' @importFrom ggplot2
#'     ggplot scale_colour_gradient2 geom_point scale_shape_manual guides labs
#'     scale_color_manual theme_classic theme element_text aes guide_legend
#'     element_blank guide_colourbar
#' @importFrom ggforce geom_arc_bar geom_circle
#' @importFrom grDevices colors
#' @importFrom stats setNames
signifPlot <- function(results,
                       fdr = FALSE,
                       type = "bubble",
                       breaks = NULL,
                       comparisonGroup = NULL,
                       colours = c("#4575B4", "white", "#D73027"),
                       marksToPlot = NULL,
                       cutoff = 0.05) {
  
  if (is.null(comparisonGroup)) {
    coef <- 2
  } else {
    coef <- which(levels(results$condition) == comparisonGroup)
  }
  marks <- unique(results$comparisons$from)
  
  if ("fromName" %in% names(results$comparisons)) {
    marks <- unique(c(
      results$comparisons$to,
      results$comparisons$from
    ))
  }
  
  
  if (is.null(marksToPlot)) marksToPlot <- marks
  
  if (type == "bubble") {
    return(
      bubblePlot(
        results, fdr, breaks, coef,
        colours = colours, cutoff = cutoff, marksToPlot = marksToPlot
      )
    )
  }
  
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
  
  pVal[isGreater] <- abs(pVal[isGreater])
  
  pVal <- matrix(pVal, nrow = length(marks))
  colnames(pVal) <- marks
  rownames(pVal) <- marks
  
  
  
  
  heatmap <- pheatmap::pheatmap(
    pVal[marksToPlot, marksToPlot],
    color = pal,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
  
  heatmap
}


bubblePlot <- function(
    test, fdr, breaks, coef,
    colours = c("blue", "white", "red"), cutoff = 0.05, marksToPlot) {
  size <- -log10(test$p.value[, coef])
  
  if (is.null(test$alternateResult)) {
    test$alternateResult <- FALSE
  }
  
  if (test$alternateResult) {
    groupA <- test$coefficient[, 1]
    groupB <- (test$coefficient[, 1] + test$coefficient[, coef])
  } else {
    groupA <- test$coefficient[, 1] * sqrt(pi) * 2 / sqrt(10) / 100
    groupB <- (
      test$coefficient[, 1] + test$coefficient[, coef]
    ) * sqrt(pi) * 2 / sqrt(10) / 100
    midpoint <- 0
  }
  cellTypeA <- factor(test$comparisons$from)
  cellTypeB <- factor(test$comparisons$to)
  
  sig <- test$p.value[, coef] < cutoff
  sigLab <- paste0("p-value < ", cutoff)
  if (fdr) {
    sig <- p.adjust(test$p.value[, coef], "fdr") < cutoff
    sigLab <- paste0("fdr < ", cutoff)
  }
  
  
  
  df <- data.frame(
    cellTypeA, cellTypeB, groupA, groupB, size,
    stat = test$statistic[, coef], pvalue = test$p.value[, coef],
    sig = factor(sig)
  )
  rownames(df) <- rownames(test$statistic)
  
  
  if ("fromName" %in% names(test$comparisons)) {
    df$cellTypeAName <- factor(test$comparisons$fromName)
    df <- df[df$cellTypeAName %in% marksToPlot & df$cellTypeB %in% marksToPlot, ] # nolint
  } else {
    df <- df[df$cellTypeA %in% marksToPlot & df$cellTypeB %in% marksToPlot, ]
  }
  
  df$cellTypeA <- droplevels(df$cellTypeA)
  df$cellTypeB <- droplevels(df$cellTypeB)
  
  shape.legend <- stats::setNames(
    c("\u25D6", "\u25D7"),
    c(levels(test$condition)[1], levels(test$condition)[coef])
  )
  
  df.shape <- data.frame(
    cellTypeA = c(NA, NA), cellTypeB = c(NA, NA), size = c(1, 1),
    condition = c(
      levels(test$condition)[1], levels(test$condition)[coef]
    )
  )
  
  
  if (test$alternateResult) {
    groupAB <- c(groupA, groupB)
    
    limits <- c(min(groupAB, na.rm = TRUE), max(groupAB, na.rm = TRUE))
    
    range <- limits[2] - limits[1]
    breaks <- c(
      limits[1] + range / 5,
      limits[1] + 2 * (range / 5),
      limits[1] + 3 * (range / 5),
      limits[1] + 4 * (range / 5),
      limits[2]
    )
    
    midpoint <- limits[1] + 3 * (range / 5)
  } else if (is.null(breaks)) {
    groupAB <- c(groupA, groupB)
    
    limits <- c(min(groupAB, na.rm = TRUE), max(groupAB, na.rm = TRUE)) |>
      round(1)
    
    breaks <- seq(from = limits[1], to = limits[2], by = diff(limits) / 5)
  } else {
    limits <- c(breaks[1], breaks[2])
    breaks <- seq(from = breaks[1], to = breaks[2], by = breaks[3])
  }
  
  
  df$groupA <- pmax(pmin(df$groupA, limits[2]), limits[1])
  df$groupB <- pmax(pmin(df$groupB, limits[2]), limits[1])
  
  pal <- grDevices::colorRampPalette(colours)(length(breaks)) # nolint
  
  labels <- round(breaks, 3)
  labels[1] <- "avoidance"
  labels[length(labels)] <- "attraction"
  
  if(.Platform$OS.type == "windows") {
    grDevices::windowsFonts(sans="Lucida Sans Unicode")
    extrafont::loadfonts(device="all", quiet = TRUE)
  } else if (.Platform$OS.type == "unix") {
    extrafont::font_import(pattern="DejaVuSans", prompt=FALSE)
    extrafont::loadfonts(device="postscript", quiet = TRUE)
    extrafont::loadfonts(device="pdf", quiet = TRUE)
  }
  
  ggplot2::ggplot(df, ggplot2::aes(x = cellTypeA, y = cellTypeB)) +
    ggplot2::scale_fill_gradient2(
      low = colours[1], mid = colours[2], high = colours[3],
      midpoint = midpoint, breaks = breaks, labels = labels, limits = limits
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = sig), size = 0) +
    ggplot2::geom_point(ggplot2::aes(size = size), x = 100000, y = 10000000) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "white", "TRUE" = "black"), labels = c("", sigLab)
    ) +
    ggforce::geom_arc_bar(
      ggplot2::aes(
        fill = groupB, r = pmax(size / max(size, na.rm = TRUE) / 2, 0.15),
        r0 = 0, x0 = as.numeric(cellTypeA), y0 = as.numeric(cellTypeB),
        start = 0, end = pi, x = NULL, y = NULL
      ),
      color = NA
    ) +
    ggforce::geom_arc_bar(
      ggplot2::aes(
        fill = groupA, r = pmax(size / max(size, na.rm = TRUE) / 2, 0.15),
        r0 = 0, x0 = as.numeric(cellTypeA), y0 = as.numeric(cellTypeB),
        start = pi, end = 2 * pi, x = NULL, y = NULL
      ),
      colour = NA
    ) +
    ggforce::geom_circle(
      data = df[df$sig == "TRUE", ], ggplot2::aes(
        r = pmax(size / max(size, na.rm = TRUE) / 2, 0.15),
        x0 = as.numeric(cellTypeA), y0 = as.numeric(cellTypeB),
        x = NULL, y = NULL
      ), colour = "black"
    ) +
    ggplot2::geom_point(
      data = df.shape, ggplot2::aes(shape = condition), x = 10000, y = 10000
    ) +
    ggplot2::scale_shape_manual(values = shape.legend) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0),
      legend.box.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = "Cell type i", y = "Cell type j", size = "-log10 p-value",
      colour = NULL, fill = "Localisation", shape = "Condition"
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(order = 3, override.aes = list(size = 5)),
      fill = ggplot2::guide_colourbar(order = 4),
      size = ggplot2::guide_legend(order = 2),
      colour = ggplot2::guide_legend(
        order = 1, override.aes = list(size = 5, shape = 1)
      )
    ) +
    coord_fixed()
}

#' Plots survival results from spicy.
#'
#' @param result A spicyResults object that contains survival results.
#' @param cutoff Significance threshold for circles in bubble plot.
#' @param colourGradient A vector of colours, used to define the low, medium, and high values for the colour scale.
#' @param marksToPlot Vector of marks to include in bubble plot.
#' @param kontextual A boolean indicating if the `result` is a Kontextual result from \code{\link[Statial]{Kontextual}}.
#' @param contextColours Only used if `kontextual = TRUE`. A named list specifying the colours for each context.
#'  By default the Tableau colour palette is used.
#' @param contextLabels Only used if `kontextual = TRUE`. A named list to change the default labels for each context.
#' 
#' @return A ggplot object.
#' @export
#' 
#' @importFrom ggh4x strip_themed facet_grid2
#' @import ggplot2
#' @import ggthemes
survBubble = function(result,
                      cutoff = 0.05,
                      colourGradient = c("#4575B4", "white", "#D73027"),
                      marksToPlot = NULL,
                      kontextual = FALSE,
                      contextColours = NULL,
                      contextLabels = waiver()){
  
  
  if(!"survivalResults" %in% names(result)) {
    stop("Survival results are missing, please run spicy with survival outcomes.")
  }
  
  survivalResults = result$survivalResults
  
  
  if(kontextual) {
    plotData = survivalResults |>
      separate(test,
               into = c("from", "to", "parent"),
               sep = "__") |>
      arrange(parent, to, from) |>
      mutate(toParent = paste(to, parent, sep = "__"))
    
    
    
    if(is.null(contextColours)) {
      # Defining colour palette for Context
      palette = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value
      
      if (length(unique(plotData$parent)) > 10) {
        palette = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 20`$value
      }
      
      contextColours = palette[1:length(unique(plotData$parent))]
      names(contextColours) = levels(plotData$parent)
    }
    
    # Assigning contextColours to the strip
    strip = strip_themed(background_x = elem_list_rect(fill = contextColours)
    )
    
  } else {
    plotData = survivalResults |>
      separate(test, into = c("from", "to"), sep = "__") |>
      arrange(to, from)
  }
  
  if(!is.null(marksToPlot)) {
    plotData = plotData |> 
      filter(to %in% marksToPlot) |> 
      filter(from %in% marksToPlot)
  }
  
  plotData = plotData |>
    mutate(
      sig = p.value < cutoff,
      logP = -log10(p.value),
      size = logP / max(logP, na.rm = TRUE),
      from = factor(from),
      to = factor(to)
    )
  
  plot = ggplot(plotData, aes(x = to, y = from)) +
    ggplot2::geom_point(ggplot2::aes(size = pmax(logP/2, 0.15), colour = coef)) +  
    ggplot2::geom_point(data = dplyr::filter(plotData, sig == TRUE),
                        aes(size = pmax(logP/2, 0.15)),
                        shape = 21) +
    ggplot2::geom_point(ggplot2::aes(shape = paste0("P < ", cutoff)), size = -1) + 
    scale_colour_gradient2(low = colourGradient[[1]],
                           mid = colourGradient[[2]],
                           high = colourGradient[[3]],
                           midpoint = 0) +
    scale_size(range = c(2, 6)) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(colour = "CoxPH \ncoefficient",
         shape = NULL,
         x = NULL,
         y = NULL) +
    ggplot2::guides(shape = ggplot2::guide_legend(order = 2, override.aes = list(size=5, shape = 1, col = "black")),
                    colour = ggplot2::guide_colourbar(order = 3),
                    size = "none") +
    theme_classic() 
  
  
  if(kontextual) {
    # Adding context information to the plot
    plot = plot + 
      geom_tile(aes(fill = parent), alpha = -1) +
      ggh4x::facet_grid2(~parent, scales = "free", space = "free", strip = strip) +
      scale_fill_manual(values = contextColours,
                        labels = contextLabels) +
      labs(fill = "Context") +
      ggplot2::guides(fill = guide_legend(order = 1, override.aes = list(alpha = 1))) +
      theme(strip.text = element_text(size = -1),
            strip.clip = "off",
            strip.background = element_rect(linewidth = NA),
            panel.spacing = unit(0.4,'lines'))
  }
  
  
  return(plot)
  
}
