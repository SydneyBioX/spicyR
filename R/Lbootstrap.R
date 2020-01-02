Lbootstrap <- function(cells, 
                       from, 
                       to, 
                       rinterval=c(0,100),
                       nsim=199,
                       global=FALSE,
                       alternative=c("two.sided", "less", "greater"),
                       seed=42) {
  
    #TODO: Produce error messages
    # e.g. Missing entries in 'cells', 'from' or 'to' not present in cellTypes

    pppCell <- pppGenerate(cells)
    
    set.seed(seed)
    
    E <- spatstat::envelope(pppCell,
                            Lcross,
                            nsim = nsim,
                            from = from,
                            to = to,
                            simulate = expression(rlabel(pppCell)),
                            correction = "best",
                            global = global, 
                            savepatterns = TRUE)
    
    test <- spatstat::dclf.test(pppCell,
                                rinterval = rinterval,
                                alternative = alternative)
    
    test
}