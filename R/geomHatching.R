scale_region <- function(aesthetics = "region", ..., guide = "legend") {
    discrete_scale("region", "region_d", palette = function(n) 1:n, ...)
}

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




GeomHatching <- ggproto("GeomHatching", GeomPoint, extra_params = c("na.rm", "line.spacing", 
    "window", "nbp", "line.width"), draw_panel = function(data, panel_params, coord, 
    na.rm = FALSE, line.spacing = 21, window = FALSE, window.length = 21, nbp = 250, line.width = 1) {
    coords <- coord$transform(data, panel_params)
    
    if (is.factor(coords$region)) 
        coords$region <- as.numeric(coords$region)
    if (is.character(coords$region)) 
        coords$region <- as.numeric(as.factor(coords$region))
    
    ow <- makeWindow(coords, window, window.length)
    
    ow <- makeWindow(coords, window) 
    
    pp <- ppp(coords$x, coords$y, window = ow, marks = coords$region)
    
    pp$region <- pp$marks
    grob <- plotRegions(pp, line.spacing, c(0, 1), c(0, 1), nbp = nbp, line.width = line.width)
    grob$name <- "geom_hatching"
    return(grob)
}, draw_key = draw_key_region, required_aes = c("x", "y", "region"), non_missing_aes = c("x", 
    "y", "region"), default_aes = aes(region = 0, size = 0.05, angle = 0, alpha = 1))


geom_hatching <- function(mapping = NULL, data = NULL, stat = "identity", position = "identity", 
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, line.spacing = 21, window = FALSE, window.length = 21,
    nbp = 250, line.width = 1, ...) {
    
    layer(geom = GeomHatching, mapping = mapping, data = data, stat = stat, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, params = list(na.rm = na.rm, 
            line.spacing = line.spacing, window = window, window.length = window.length, nbp = nbp, line.width = line.width, 
            ...))
}

