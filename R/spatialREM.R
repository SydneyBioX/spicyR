spatialREM <- function(x, patient, condition, count1, count2, nsim=19) {
    
    spatialData <- data.frame(Pair = x,
                              Patient = patient,
                              Condition = condition)
    
    count1 <- sqrt(count1)
    count2 <- sqrt(count2)
    
    pairwiseVals <- x
    
    z <- (pairwiseVals-mean(pairwiseVals))^2
    z1 <- mgcv::gam(z~ti(x,y))
    w <- 1/sqrt(z1$fitted.values-min(z1$fitted.values)+1)
    w <- w/sum(w)
    
    mixed.lmer <- lme4::lmer(Pair ~ Group + (1|Patient), 
                             data = spatialData, 
                             weights = w)
    
    m <- fixef(mixed.lmer)[2]
}