Lplot <- function(cells, from, to, rinterval=c(0,100)) {
    #TODO: Produce error messages
    # e.g. Missing entries in 'cells', 'from' or 'to' not present in cellTypes
    #TODO: make plot with ggplot
    
    pppCell <- spatstat::ppp(cells$x,
                             cells$y,
                             xrange=c(0,max(x)),
                             yrange=c(0,max(y)),
                             marks=cells$cellType)
    
    plot(spatstat::Lcross(pppCell, from=from, to=to), xrange=rinterval)
}