#' Plots result of signifPlot.
#'
#' @param results A spicy results object
#' @param type Whether to make a bubble plot or heatmap. Note: For survival results a bubble plot will be used. 
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
#' @return a ggplot or pheatmap object
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
#' @import ggthemes
#' @importFrom ggh4x strip_themed facet_grid2
#' @importFrom grid grobTree polygonGrob gpar
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
  marks <- unique(results$comparisons$from)
  
  if ("fromName" %in% names(results$comparisons)) {
    marks <- unique(c(
      results$comparisons$to,
      results$comparisons$from
    ))
  }
  
  
  if (is.null(marksToPlot)) marksToPlot <- marks
  
  
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

# Function to create x and y points for a half circle (left or right) ---> CHANGE HERE
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

# Custom draw_key function to draw a left half circle in the legend using polygonGrob ----> CHANGE HERE
draw_key_half_circle = function(data, params, shape) {
  if(data$shape == 16){
    coords <- half_circle_coords(shape = "left")
    grid::grobTree(
      grid::polygonGrob(
        x = coords$x, y = coords$y,
        gp = grid::gpar(fill = "black", col = "black")
      )
    )
  } else{
    coords <- half_circle_coords(shape = "right")
    grid::grobTree(
      grid::polygonGrob(
        x = coords$x, y = coords$y,
        gp = grid::gpar(fill = "black", col = "black")
      )
    )
  }
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

# Custom draw_key function to draw a left half circle in the legend using polygonGrob 
draw_key_half_circle = function(data, params, shape) {
  if(data$shape == 16){
    coords <- half_circle_coords(shape = "left")
    grid::grobTree(
      grid::polygonGrob(
        x = coords$x, y = coords$y,
        gp = grid::gpar(fill = "black", col = "black")
      )
    )
  } else{
    coords <- half_circle_coords(shape = "right")
    grid::grobTree(
      grid::polygonGrob(
        x = coords$x, y = coords$y,
        gp = grid::gpar(fill = "black", col = "black")
      )
    )
  }
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
  
  if (test$alternateResult ) {
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
    sigLab <- paste0("fdr < ", cutoff)
  }
  
  size <- -log10(pvalue)
  
  
  df <- data.frame(
    cellTypeA, cellTypeB, groupA, groupB, size,
    stat = test$statistic[, coef], pvalue = pvalue,
    sig = factor(sig, levels= c("FALSE", "TRUE"))
  )
  rownames(df) <- rownames(test$statistic)
  
  
  if (isTRUE(test$isKontextual)) {
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
  
  if(test$alternateResult && !isTRUE(test$isKontextual)){
    midpoint <- (breaks[1] + breaks[length(breaks)])/2
  }
  
  
  df$groupA <- pmax(pmin(df$groupA, limits[2]), limits[1])
  df$groupB <- pmax(pmin(df$groupB, limits[2]), limits[1])
  
  pal <- grDevices::colorRampPalette(colours)(length(breaks)) # nolint
  
  labels <- round(breaks, 3)
  labels[1] <- "avoidance"
  labels[length(labels)] <- "attraction"


  if(isTRUE(test$isKontextual)) {
    df = df |> 
      group_by(parent) |>
      mutate(
        cellTypeB_numeric = as.numeric(droplevels(cellTypeB)),
        cellTypeB_id = factor(paste(parent, cellTypeB, sep = "_"))
      ) |>
      ungroup()
    
  } else {
    
    df = df |> 
      mutate(
        cellTypeB_numeric = as.numeric(cellTypeB),
        cellTypeB_id = cellTypeB
      )
  }
  
  plot = ggplot2::ggplot(df, ggplot2::aes(x = cellTypeB_id, y = cellTypeA)) +
    ggplot2::scale_fill_gradient2(
      low = colours[1], mid = colours[2], high = colours[3],
      midpoint = midpoint, breaks = breaks, labels = labels, limits = limits
    ) +
    ggplot2::geom_point(ggplot2::aes(col = sigLab), size = -1) +
    ggplot2::geom_point(ggplot2::aes(size = size), x = 100000, y = 10000000) +
    ggforce::geom_arc_bar(
      ggplot2::aes(
        fill = groupB, r = pmax(size / max(size, na.rm = TRUE) / 2, 0.15),
        r0 = 0, x0 = cellTypeB_numeric, y0 = as.numeric(cellTypeA),
        start = 0, end = pi, x = NULL, y = NULL
      ),
      color = NA
    ) +
    ggforce::geom_arc_bar(
      ggplot2::aes(
        fill = groupA, r = pmax(size / max(size, na.rm = TRUE) / 2, 0.15),
        r0 = 0, x0 = cellTypeB_numeric, y0 = as.numeric(cellTypeA),
        start = pi, end = 2 * pi, x = NULL, y = NULL
      ),
      colour = NA
    ) +
    ggplot2::geom_point(
      data = df.shape, ggplot2::aes(shape = condition), x = 10000, y = 10000,
      key_glyph = draw_key_half_circle 
    ) +
    scale_x_discrete(breaks = df$cellTypeB_id, labels = df$cellTypeB,
                     guide = guide_axis(angle = 45)) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Cell type j", y = "Cell type i", size = "-log10 p-value",
      colour = NULL, fill = "Localisation", shape = "Condition"
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(order = 3,override.aes = list(size = 5)),
      #fill = ggplot2::guide_colourbar(order = 4), # adding this line breaks the ggnewscales package
      size = ggplot2::guide_legend(order = 2),
      colour = ggplot2::guide_legend(
        order = 1, override.aes = list(size = 5, shape = 1, col = "black")
      )
    )
  
  
  # Plots black circle outlines, only if there are significant results.
  if(nrow(df[df$sig == "TRUE", ]) > 0) {
    plot = plot + 
      ggforce::geom_circle(
        data = df[df$sig == "TRUE", ], ggplot2::aes(
          r = pmax(size / max(size, na.rm = TRUE) / 2, 0.15),
          x0 = cellTypeB_numeric, y0 = as.numeric(cellTypeA),
          x = NULL, y = NULL
        ), colour = "black"
      ) 
  }
  
  # Adds context panels if using Kontextual results.
  if(isTRUE(test$isKontextual)) {
    if(is.null(contextColours)) {
      # Defining colour palette for Context
      palette = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value
      
      if (length(unique(df$parent)) > 10) {
        palette = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 20`$value
      }
      
      contextColours = palette[1:length(unique(df$parent))]
      names(contextColours) = levels(df$parent)
    }
    
    # Assigning contextColours to the strip
    strip = strip_themed(background_x = elem_list_rect(fill = contextColours))
    
    plot = plot +
      ggh4x::facet_grid2( ~ parent,
                          scales = "free",
                          space = "free",
                          strip = strip) +
      theme(
        strip.text = element_text(size = -1),
        strip.clip = "off",
        strip.background = element_rect(linewidth = NA),
        panel.spacing = unit(0.4, 'lines')
      ) +
      new_scale("fill") +
      geom_tile(aes(fill = parent), alpha = -1) +
      scale_fill_manual(values = contextColours, labels = contextLabels) +
      labs(fill = "Context") +
      ggplot2::guides(fill = guide_legend(order = 6, override.aes = list(alpha = 1)))
    
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
#' @importFrom ggh4x strip_themed facet_grid2
#' @import ggplot2
#' @import ggthemes
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
  
  survivalResults = result$survivalResults

  
  if(isTRUE(result$isKontextual)) {
    plotData = survivalResults |>
      tidyr::separate(test,
               into = c("from", "to", "parent"),
               sep = "__") |>
      dplyr::arrange(parent, to, from) |>
      dplyr::mutate(toParent = paste(to, parent, sep = "__"))
    
    
    
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
      tidyr::separate(test, into = c("from", "to"), sep = "__") |>
      dplyr::arrange(to, from)
  }
  
  if(!is.null(marksToPlot)) {
    plotData = plotData |> 
      filter(to %in% marksToPlot) |> 
      filter(from %in% marksToPlot)
  }
  
  sigLab <- paste0("p-value < ", cutoff)
  
  if(fdr){
    plotData$p.value = p.adjust(plotData$p.value, "fdr")
    sigLab <- paste0("fdr < ", cutoff)
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
    ggplot2::geom_point(ggplot2::aes(shape = sigLab), size = -1) + 
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
    ggplot2::guides(shape = ggplot2::guide_legend(order = 1, override.aes = list(size=5, shape = 1, col = "black")),
                    colour = ggplot2::guide_colourbar(order = 2),
                    size = "none") +
    theme_classic() 
  
  
  if(isTRUE(result$isKontextual)) {
    # Adding context information to the plot
    plot = plot + 
      geom_tile(aes(fill = parent), alpha = -1) +
      ggh4x::facet_grid2(~parent, scales = "free", space = "free", strip = strip) +
      scale_fill_manual(values = contextColours,
                        labels = contextLabels) +
      labs(fill = "Context") +
      ggplot2::guides(fill = guide_legend(order = 2, override.aes = list(alpha = 1))) +
      theme(strip.text = element_text(size = -1),
            strip.clip = "off",
            strip.background = element_rect(linewidth = NA),
            panel.spacing = unit(0.4,'lines'))
  }
  
  
  return(plot)
  
}
