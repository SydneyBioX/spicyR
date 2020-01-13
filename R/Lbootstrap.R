Lbootstrap <- function(x, 
                       from, 
                       to, 
                       imageID=1,
                       rinterval=c(0,100),
                       nsim=19,
                       global=FALSE,
                       alternative=c("two.sided", "less", "greater"),
                       seed=42) {
  
    #TODO: Produce error messages
    # e.g. Missing entries in 'cells', 'from' or 'to' not present in cellTypes
  
    if (is(x, "SegmentedCellExperiment")) {
      cells <- location(x)
      cells <- cells[cells$imageID==as.character(imageID),]
    } else if (is(x, "data.frame")) {
      cells <- x
      if (is.null(cells$imageID)) {
        cells$imageID <- rep(imageID, nrow(cells))
      }
    } else {
      stop("x is not a data.frame or SegmentedCellExperiment")
    }
    
    
    pppCell <- pppGenerate(cells)
    
    set.seed(seed)
    
    E <- spatstat::envelope(pppCell,
                            spatstat::Lcross,
                            nsim = nsim,
                            from = from,
                            to = to,
                            simulate = expression(spatstat::rlabel(pppCell)),
                            correction = "best",
                            global = global, 
                            savepatterns = TRUE)
    
    test <- spatstat::dclf.test(E,
                                rinterval = rinterval,
                                alternative = alternative)
    
    test
}