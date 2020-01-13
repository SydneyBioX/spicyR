#' Hatching geom
#' 
#' The hatching geom is used to create hatching patterns for representation of spatial regions.
#'
#' @param mapping Set of aesthetic mappings created by aes() or aes_(). If specified and inherit.aes = TRUE (the default), it is combined with the default mapping at the top level of the plot. You must supply mapping if there is no plot mapping.
#' @param data The data to be displayed in this layer. There are three options: 
#' 
#' If NULL, the default, the data is inherited from the plot data as specified in the call to ggplot().
#' A data.frame, or other object, will override the plot data. All objects will be fortified to produce a data frame. See fortify() for which variables will be created. 
#' A function will be called with a single argument, the plot data. The return value must be a data.frame, and will be used as the layer data. A function can be created from a formula (e.g. ~ head(.x, 10)).
#' @param stat The statistical transformation to use on the data for this layer, as a string.
#' @param Position adjustment, either as a string, or the result of a call to a position adjustment function.
#' @param show.legend logical. Should this layer be included in the legends? NA, the default, includes if any aesthetics are mapped. FALSE never includes, and TRUE always includes. It can also be a named logical vector to finely select the aesthetics to display.
#' @param inherit.aes If FALSE, overrides the default aesthetics, rather than combining with them. This is most useful for helper functions that define both data and aesthetics and shouldn't inherit behaviour from the default plot specification, e.g. borders().
#' @param na.rm If FALSE, the default, missing values are removed with a warning. If TRUE, missing values are silently removed.
#' @param ling.spacing A integer indicating the spacing between hatching lines.
#' @param window Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length A tuning parameter for controlling the level of concavity when estimating concave windows.
#' @param nbp An integer tuning the granularity of the grid used when defining regions
#' @param line.width A numeric controlling the width of the hatching lines
#' 
#' @examples
#' 1+1
#' @export
#' @rdname lisa
#' @importFrom methods is
#' @importFrom BiocParallel bplapply
#' @import SegmentedCellExperiment
geom_hatching <- function(mapping = NULL, data = NULL, stat = "identity", position = "identity", 
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, line.spacing = 21, window = "square", 
    window.length = 0, nbp = 250, line.width = 1, ...) {
    
    layer(geom = GeomHatching, mapping = mapping, data = data, stat = stat, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, params = list(na.rm = na.rm, 
            line.spacing = line.spacing, window = window, window.length = window.length, 
            nbp = nbp, line.width = line.width, ...))
}


#' @export
scale_region <- function(aesthetics = "region", ..., guide = "legend") {
    discrete_scale("region", "region_d", palette = function(n) 1:n, ...)
}

#' @export
scale_region_manual <- function(..., values) {
    
    
    force(values)
    pal <- function(n) {
        if (n > length(values)) {
            stop("Insufficient values in manual scale. ", n, " needed but only ", 
                length(values), " provided.", call. = FALSE)
        }
        if (any(!values %in% 1:7)) {
            stop("values must be between 1 and 7")
        }
        
        values
    }
    discrete_scale("region", "manual", pal, ...)
    
}


ggname <- getFromNamespace("ggname", "ggplot2")

#' @importFrom grid grob polylineGrob gpar
draw_key_region <- function(data, params, size) {
    grobs <- grob()
    
    
    if (data$region == 1) {
        grobs <- polylineGrob(x = c(0, 0, 1, 1, 0), y = c(0, 1, 1, 0, 0), id = c(1, 
            1, 1, 1, 1), gp = gpar(col = 1, lwd = 1))
    }
    
    if (data$region == 2) {
        grobs <- polylineGrob(x = c(c(0, 0.5), c(0, 1), c(0.5, 1), c(0, 0, 1, 1, 
            0)), y = c(c(0.5, 1), c(0, 1), c(0, 0.5), c(0, 1, 1, 0, 0)), id = c(1, 
            1, 2, 2, 3, 3, 4, 4, 4, 4, 4), gp = gpar(col = 1, lwd = 1))
    }
    if (data$region == 3) {
        grobs <- polylineGrob(x = c(c(1, 0.5), c(1, 0), c(0.5, 0), c(0, 0, 1, 1, 
            0)), y = c(c(0.5, 1), c(0, 1), c(0, 0.5), c(0, 1, 1, 0, 0)), id = c(1, 
            1, 2, 2, 3, 3, 4, 4, 4, 4, 4), gp = gpar(col = 1, lwd = 1))
        
    }
    
    if (data$region == 4) {
        grobs <- polylineGrob(x = c(c(0, 1), c(0, 1), c(0, 0, 1, 1, 0)), y = c(c(0.33, 
            0.33), c(0.66, 0.66), c(0, 1, 1, 0, 0)), id = c(1, 1, 2, 2, 3, 3, 3, 
            3, 3), gp = gpar(col = 1, lwd = 1))
        
    }
    
    if (data$region == 5) {
        grobs <- polylineGrob(x = c(c(0.33, 0.33), c(0.66, 0.66), c(0, 0, 1, 1, 0)), 
            y = c(c(0, 1), c(0, 1), c(0, 1, 1, 0, 0)), id = c(1, 1, 2, 2, 3, 3, 3, 
                3, 3), gp = gpar(col = 1, lwd = 1))
        
    }
    
    
    if (data$region == 6) {
        grobs <- polylineGrob(x = c(c(1, 0.5), c(1, 0), c(0.5, 0), c(0, 0.5), c(0, 
            1), c(0.5, 1), c(0, 0, 1, 1, 0)), y = c(c(0.5, 1), c(0, 1), c(0, 0.5), 
            c(0.5, 1), c(0, 1), c(0, 0.5), c(0, 1, 1, 0, 0)), id = c(1, 1, 2, 2, 
            3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 7, 7), gp = gpar(col = 1, lwd = 1))
    }
    
    if (data$region == 7) {
        grobs <- polylineGrob(x = c(c(0, 1), c(0, 1), c(0.33, 0.33), c(0.66, 0.66), 
            c(0, 0, 1, 1, 0)), y = c(c(0.33, 0.33), c(0.66, 0.66), c(0, 1), c(0, 
            1), c(0, 1, 1, 0, 0)), id = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5), 
            gp = gpar(col = 1, lwd = 1))
        
    }
    grobs$name <- "region_key"
    grobs
    
    
}



hatchingLevels <- function(data, hatching = NULL) {
    
    if (!is.factor(data$region)) 
        data$region <- factor(data$region)
    regionLevels <- levels(data$region)
    if (!any(hatching %in% 1:7) & !is.null(hatching)) 
        stop("hatching must equal the number of regions and be <= 7.")
    if (all(regionLevels %in% names(hatching))) {
        hatching <- hatching[regionLevels]
    }
    if (is.null(hatching)) {
        hatching <- 1:length(regionLevels)
        names(hatching) <- regionLevels
    }
    if (length(hatching) == length(regionLevels)) {
        names(hatching) <- regionLevels
    } else {
        stop("hatching must be the same length as the number of regions.")
    }
    return(hatching)
}




GeomHatching <- ggplot2::ggproto("GeomHatching", ggplot2::GeomPoint, extra_params = c("na.rm", 
    "line.spacing", "window", "nbp", "line.width"), draw_panel = function(data, panel_params, 
    coord, na.rm = FALSE, line.spacing = 21, window = "square", window.length = 0, 
    nbp = 250, line.width = 1) {
    coords <- coord$transform(data, panel_params)
    
    if (is.factor(coords$region)) 
        coords$region <- as.numeric(coords$region)
    if (is.character(coords$region)) 
        coords$region <- as.numeric(as.factor(coords$region))
    
    ow <- makeWindow(coords, window, window.length)
    
    pp <- spatstat::ppp(coords$x, coords$y, window = ow, marks = coords$region)
    
    
    pp$region <- pp$marks
    grob <- plotRegions(pp, line.spacing, c(0, 1), c(0, 1), nbp = nbp, line.width = line.width)
    grob$name <- "geom_hatching"
    return(grob)
}, draw_key = draw_key_region, required_aes = c("x", "y", "region"), non_missing_aes = c("x", 
    "y", "region"), default_aes = ggplot2::aes(region = 0, size = 0.05, angle = 0, 
    alpha = 1))





