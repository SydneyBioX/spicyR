pppGenerate <- function(cells) {
    pppCell <- spatstat::ppp(cells$x,
                             cells$y,
                             xrange = c(0,max(cells$x)),
                             yrange = c(0,max(cells$y)),
                             marks = cells$cellType)
    
    pppCell
}