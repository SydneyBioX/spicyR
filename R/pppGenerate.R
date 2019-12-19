#TODO: Add contingencies
#TODO: Make this function non-accessible?

pppGenerate <- function(cells) {
    pppCell <- spatstat::ppp(cells$x,
                             cells$y,
                             xrange = c(0,max(x)),
                             ytange = c(0,max(y)),
                             marks = cells$cellType)
    
    pppCell
}