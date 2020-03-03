#' hatchingPlot
#' 
#' The hatchingPlot() function is used to create hatching patterns for representating spatial regions and cell-types.
#'
#' @param data A segmentedCells object that has region information 
#' @param imageID A vector of imageIDs to be plotted
#' @param line.spacing A integer indicating the spacing between hatching lines.
#' @param window Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length A tuning parameter for controlling the level of concavity when estimating concave windows.
#' @param nbp An integer tuning the granularity of the grid used when defining regions
#' 
#' @return A ggplot object
#' 
#' @examples
#' library(spicyR)
#' \dontrun{
#' # Given a segmentedCells object cellExp
#' # First constract LISA
#' 
#' lisaCurves <- lisa(cellExp)
#' 
#' # Then cluster the LISA curves
#' 
#' kM <- kmeans(lisaCurves,2)
#' region(cellExp) <- paste('region',kM$cluster,sep = '_')
#' 
#' # Then plot the regions
#' hatchingPlot(cellExp)
#' }
#' 
#' @export
#' @rdname lisa
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal facet_wrap
hatchingPlot <- function(data, imageID = NULL, window = "concave", line.spacing = 21, 
    nbp = 250, window.length = 0) {
    if (!is(data, "segmentedCells") | is.null(data$region)) 
        stop("Please provide a segmentedCells object with region information please.")
    
    if (is.null(imageID)) {
        df <- region(data[1, ], annot = TRUE)
        p <- ggplot(df, aes(x = x, y = y, colour = cellType)) + geom_point() + geom_hatching(aes(region = region), 
            show.legend = TRUE, window = "concave", line.spacing = 21, nbp = 250, 
            window.length = 0)
        q <- p + theme_minimal() + scale_region()
        return(q)
    }
    
    if (!is.null(imageID)) {
        if (any(!imageID %in% rownames(data))) 
            stop("Some of the imageIDs are not in your segmentedCells object")
        df <- region(data, imageID = imageID, annot = TRUE)
        p <- ggplot(df, aes(x = x, y = y, colour = cellType)) + geom_point() + facet_wrap(~imageID) + 
            geom_hatching(aes(region = region), show.legend = TRUE, window = "concave", 
                line.spacing = 21, nbp = 250, window.length = 0)
        q <- p + theme_minimal() + scale_region()
        q
    }
}


################################################################################
##
## geomHatching
##
################################################################################



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
#' @param line.spacing A integer indicating the spacing between hatching lines.
#' @param window Should the window around the regions be 'square', 'convex' or 'concave'.
#' @param window.length A tuning parameter for controlling the level of concavity when estimating concave windows.
#' @param nbp An integer tuning the granularity of the grid used when defining regions
#' @param line.width A numeric controlling the width of the hatching lines
#' 
#' 
#' @return A ggplot geom
#' 
#' @examples
#' library(spicyR)
#' \dontrun{
#' # Given a segmentedCells object cellExp
#' # First constract LISA
#' 
#' lisaCurves <- lisa(cellExp)
#' 
#' # Then cluster the LISA curves
#' 
#' kM <- kmeans(lisaCurves,2)
#' region(cellExp) <- paste('region',kM$cluster,sep = '_')
#' 
#' # Then plot the regions
#' p <- ggplot(df,aes(x = x,y = y, colour = cellType, region = region)) + 
#' geom_point() + 
#' facet_wrap(~imageID) +
#' geom_hatching()
#' }
#' 
#' @export
#' @rdname lisa
#' @importFrom methods is
#' @importFrom BiocParallel bplapply
#' @importFrom ggplot2 layer
geom_hatching <- function(mapping = NULL, data = NULL, stat = "identity", position = "identity", 
                          na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, line.spacing = 21, window = "square", 
                          window.length = 0, nbp = 250, line.width = 1, ...) {
    
    ggplot2::layer(geom = GeomHatching, mapping = mapping, data = data, stat = stat, position = position, 
          show.legend = show.legend, inherit.aes = inherit.aes, params = list(na.rm = na.rm, 
                                                                              line.spacing = line.spacing, window = window, window.length = window.length, 
                                                                              nbp = nbp, line.width = line.width, ...))
}


#' @export
#' @importFrom ggplot2 discrete_scale
scale_region <- function(aesthetics = "region", ..., guide = "legend") {
    discrete_scale("region", "region_d", palette = function(n) seq_len(n), ...)
}

#' @export
#' @importFrom ggplot2 discrete_scale
scale_region_manual <- function(..., values) {
    
    
    force(values)
    pal <- function(n) {
        if (n > length(values)) {
            stop("Insufficient values in manual scale. ", n, " needed but only ", 
                 length(values), " provided.", call. = FALSE)
        }
        if (any(!values %in% seq_len(7))) {
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
    if (!any(hatching %in% seq_len(7)) & !is.null(hatching)) 
        stop("hatching must equal the number of regions and be <= 7.")
    if (all(regionLevels %in% names(hatching))) {
        hatching <- hatching[regionLevels]
    }
    if (is.null(hatching)) {
        hatching <- seq_len(length(regionLevels))
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





################################################################################
##
## Create Hatchings
##
################################################################################



######## Plot the hatchings
#' @importFrom purrr map
#' @importFrom grid linesGrob gpar gList grobTree
plotRegions <- function(pp, line.spacing = 21, xrange = c(0, 1), yrange = c(0, 1), 
                        nbp = 250, line.width = 1) {
    
    rx <- xrange
    ry <- yrange
    width <- (rx[2] - rx[1])/line.spacing
    
    # Convert ppp to grid
    rG <- regionGrid(pp, nbp)
    
    # Convert grid to polygon
    tree <- purrr::map(as.character(unique(pp$region)), ~{
        
        rPoly <- regionPoly(rG, .)
        
        bdrys <- purrr::map(rPoly$bdry, ~{
            df <- do.call("cbind", .)
            df <- data.frame(rbind(df, df[1, ]))
            g <- linesGrob(x = df$x/rx[2], y = df$y/ry[2], gp = gpar(col = 1, lwd = line.width))
            return(g)
        })
        
        hatchFun <- switch(., `1` = hatchNull, `2` = hatch45, `3` = hatch315, `4` = hatch90, 
                           `5` = hatch180, `6` = hatchX, `7` = hatchPlus)
        return(c(bdrys, hatchFun(rPoly, width, rx, ry, line.width = line.width)))
    })
    
    g <- do.call("gList", (do.call("c", tree)))
    return(grobTree(g))
}






######## Map predicted regions to a regular grid
#' @importFrom class knn
#' @importFrom grid linesGrob gpar gList
#' @importFrom spatstat as.mask
regionGrid <- function(pp, nbp = 250) {
    ow <- pp$window
    m <- as.mask(ow, dimyx = c(nbp, nbp))$m
    x <- seq(from = ow$xrange[1], to = ow$xrange[2], length.out = nrow(m))
    y <- seq(from = ow$yrange[1], to = ow$yrange[2], length.out = ncol(m))
    grid <- expand.grid(x = x, y = y)
    grid <- data.frame(x = grid[, 1], y = grid[, 2])
    df <- as.data.frame(pp)
    k <- rep(NA, length(m))
    K <- knn(train = df[, c("x", "y")], test = grid[t(m)[seq_len(length(m))], ], cl = pp$region, 
             k = 1)
    k[t(m)[seq_len(length(m))]] <- as.character(K)
    data.frame(grid, regions = k)
}


######## Convert grid of regions into a polygon mask for a particular region.
#' @importFrom spatstat owin as.polygonal
regionPoly <- function(grid, region) {
    rx <- range(grid$x)
    ry <- range(grid$y)
    mat <- matrix(grid$regions == region, nrow = length(unique(grid$x)), ncol = length(unique(grid$y)), 
                  byrow = TRUE)
    mat[is.na(mat)] <- FALSE
    ow <- owin(xrange = rx, yrange = ry, mask = mat)
    return(as.polygonal(ow))
}


######## Create line grobs of the hatching
#' @importFrom purrr map map_dfr
#' @importFrom grid linesGrob

hatchingLines <- function(rPoly, allHatch, ordColumn, h90 = FALSE, xr, yr, rx, ry, 
                          line.width = 1) {
    
    purrr::map(seq_len(nrow(allHatch)), ~{
        
        if (h90) {
            hatch <- rbind(c(xr[1], allHatch[., "from"]), c(xr[2], allHatch[., "to"]))
        } else {
            hatch <- rbind(c(allHatch[., "from"], yr[1]), c(allHatch[., "to"], yr[2]))
        }
        
        
        intPoints <- purrr::map_dfr(rPoly$bdry, ~{
            df <- do.call("cbind", .)
            df <- rbind(df, df[1, ])
            colnames(df) <- c("x", "y")
            
            int <- purrr::map_dfr(seq_len(nrow(df) - 1), ~{
                x1 <- df[., ]
                x2 <- df[. + 1, ]
                return(data.frame(t(line.intersection(x1, x2, hatch[1, ], hatch[2, 
                ], interior.only = TRUE))))
            })
            
            x1 <- df[nrow(df), ]
            x2 <- df[1, ]
            int <- rbind(int, line.intersection(hatch[1, ], hatch[2, ], x1, x2, interior.only = TRUE))
            colnames(int) <- c("x", "y")
            return(int)
        })
        
        intPoints <- unique(round(intPoints, 9))
        intPoints <- intPoints[!rowSums(intPoints) %in% c("Inf", NA), ]
        linesH <- NULL
        
        if (length(unlist(intPoints)) > 2) {
            intPoints <- intPoints[order(intPoints[, ordColumn]), ]
            from <- sort(rep(seq(1, nrow(intPoints) - 1, by = 2), 2))
            pointSplit <- split(intPoints, from)
            
            linesH <- purrr::map(pointSplit, ~{
                return(linesGrob(x = .$x/rx[2], y = .$y/ry[2], gp = gpar(col = 1, 
                                                                         lwd = line.width)))
                
            })
        }
        return(linesH)
    })
}


######## No hatching

hatchNull <- function(rPoly, width, rx, ry, line.width = 1) {
    NULL
}

######## / hatching


hatch45 <- function(rPoly, width, rx, ry, line.width = 1) {
    xr <- range(rPoly$x)
    yr <- range(rPoly$y)
    from1 <- seq(from = xr[1], to = xr[2], by = width)
    to1 <- seq(from = xr[2], to = xr[2] + xr[2] - xr[1], by = width)
    from2 <- seq(from = xr[1], to = 2 * xr[1] - xr[2], by = -width)
    to2 <- seq(from = xr[2], to = xr[1], by = -width)
    allHatch <- data.frame(from = c(from1, from2), to = c(to1, to2))
    
    lines <- hatchingLines(rPoly, allHatch, ordColumn = 1, h90 = FALSE, xr, yr, rx, 
                           ry, line.width = line.width)
    
    return(do.call("c", lines))
}

######## \ hatching

hatch315 <- function(rPoly, width, rx, ry, line.width = 1) {
    xr <- range(rPoly$x)
    yr <- range(rPoly$y)
    from1 <- seq(from = xr[2], to = xr[1], by = -width)
    to1 <- seq(from = xr[1], to = xr[1] + xr[1] - xr[2], by = -width)
    from2 <- seq(from = xr[2], to = 2 * xr[2] - xr[1], by = width)
    to2 <- seq(from = xr[1], to = xr[2], by = width)
    allHatch <- data.frame(from = c(from1, from2), to = c(to1, to2))
    
    lines <- hatchingLines(rPoly, allHatch, ordColumn = 1, h90 = FALSE, xr, yr, rx, 
                           ry, line.width = line.width)
    
    return(do.call("c", lines))
    
}

######## | hatching

hatch180 <- function(rPoly, width, rx, ry, line.width = 1) {
    xr <- range(rPoly$x)
    yr <- range(rPoly$y)
    from1 <- seq(from = xr[1], to = xr[2], by = width)
    to1 <- seq(from = xr[1], to = xr[2], by = width)
    allHatch <- data.frame(from = c(from1), to = c(to1))
    
    lines <- hatchingLines(rPoly, allHatch, ordColumn = 2, h90 = FALSE, xr, yr, rx, 
                           ry, line.width = line.width)
    
    return(do.call("c", lines))
}


######## - hatching

hatch90 <- function(rPoly, width, rx, ry, line.width = 1) {
    xr <- range(rPoly$x)
    yr <- range(rPoly$y)
    from1 <- seq(from = yr[1], to = yr[2], by = width)
    to1 <- seq(from = yr[1], to = yr[2], by = width)
    allHatch <- data.frame(from = c(from1), to = c(to1))
    
    lines <- hatchingLines(rPoly, allHatch, ordColumn = 1, h90 = TRUE, xr, yr, rx, 
                           ry, line.width = line.width)
    
    return(do.call("c", lines))
    
}

######## x hatching

hatchX <- function(rPoly, width, rx, ry, line.width = 1) {
    h45 <- hatch45(rPoly, width, rx, ry, line.width = line.width)
    h315 <- hatch315(rPoly, width, rx, ry, line.width = line.width)
    return(c(h45, h315))
}

######## + hatching

hatchPlus <- function(rPoly, width, rx, ry, line.width = 1) {
    h90 <- hatch90(rPoly, width, rx, ry, line.width = line.width)
    h180 <- hatch180(rPoly, width, rx, ry, line.width = line.width)
    return(c(h90, h180))
}



####### Calculate intersection of two lines.

line.intersection <- function(P1, P2, P3, P4, interior.only = TRUE) {
    ## Modified from the retistruct package to address edge cases.
    P1 <- round(as.vector(P1), 10)
    P2 <- round(as.vector(P2), 10)
    P3 <- round(as.vector(P3), 10)
    P4 <- round(as.vector(P4), 10)
    dx1 <- P1[1] - P2[1]
    dx2 <- P3[1] - P4[1]
    dy1 <- P1[2] - P2[2]
    dy2 <- P3[2] - P4[2]
    D <- det(rbind(c(dx1, dy1), c(dx2, dy2)))
    if (is.na(D) | D == 0) {
        return(c(Inf, Inf))
    }
    D1 <- det(rbind(P1, P2))
    D2 <- det(rbind(P3, P4))
    X <- round(det(rbind(c(D1, dx1), c(D2, dx2)))/D, 10)
    Y <- round(det(rbind(c(D1, dy1), c(D2, dy2)))/D, 10)
    if (interior.only) {
        lambda1 <- -((X - P1[1]) * dx1 + (Y - P1[2]) * dy1)/(dx1^2 + dy1^2)
        lambda2 <- -((X - P3[1]) * dx2 + (Y - P3[2]) * dy2)/(dx2^2 + dy2^2)
        if (!((lambda1 >= 0) & (lambda1 <= 1) & (lambda2 >= 0) & (lambda2 <= 1))) {
            return(c(NA, NA))
        }
    }
    return(c(X, Y))
}
