######## Plot the hatchings

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
        cat("before hatchFun", ., "\n")
        hatchFun <- switch(., `1` = hatchNull, `2` = hatch45, `3` = hatch315, `4` = hatch90, 
            `5` = hatch180, `6` = hatchX, `7` = hatchPlus)
        return(c(bdrys, hatchFun(rPoly, width, rx, ry, line.width = line.width)))
    })
    
    g <- do.call("gList", (do.call("c", tree)))
    return(grobTree(g))
}






######## Map predicted regions to a regular grid

regionGrid <- function(pp, nbp = 250) {
    ow <- pp$window
    m <- as.mask(ow, dimyx = c(nbp, nbp))$m
    x <- seq(from = ow$xrange[1], to = ow$xrange[2], length.out = nrow(m))
    y <- seq(from = ow$yrange[1], to = ow$yrange[2], length.out = ncol(m))
    Grid <- expand.grid(x = x, y = y)
    Grid <- data.frame(x = Grid[, 1], y = Grid[, 2])
    df <- as.data.frame(pp)
    k <- rep(NA, length(m))
    K <- knn(train = df[, c("x", "y")], test = Grid[t(m)[1:length(m)], ], cl = pp$region, 
        k = 1)
    k[t(m)[1:length(m)]] <- as.character(K)
    return(data.frame(Grid, regions = k))
}


######## Convert grid of regions into a polygon mask for a particular region.

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

hatch90 <- function(rPoly, width, rx, ry, line.width = l) {
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
