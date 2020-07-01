## ----warning=FALSE, message=FALSE---------------------------------------------
# load required packages
library(spicyR)
library(ggplot2)

## ----message=FALSE------------------------------------------------------------
data("melanomaResponders")
melanomaResponders
location(melanomaResponders)
phenotype(melanomaResponders)

## -----------------------------------------------------------------------------
spicyTestPair <- spicy(melanomaResponders, 
                   condition = "condition", 
                   subject = "subject", 
                   from = "CD8+PD1+PDL1-", 
                   to = "CD8-PD1+PDL1+")
spicyTestPair
top(spicyTestPair)

## ----echo=FALSE---------------------------------------------------------------
data("spicyTest")

## ----eval = FALSE-------------------------------------------------------------
#  spicyTest <- spicy(melanomaResponders,
#                     condition = "condition",
#                     subject = "subject")
#  

## -----------------------------------------------------------------------------
spicyTest
top(spicyTest)  

## -----------------------------------------------------------------------------
signifPlot(spicyTest, breaks=c(-3, 3, 0.5))

## ----echo=FALSE---------------------------------------------------------------
data(spicyTestBootstrap)

## ----eval = FALSE-------------------------------------------------------------
#  spicyTestBootstrap <- spicy(melanomaResponders,
#                     condition = "condition",
#                     subject = "subject",
#                     nsim = 199)

## -----------------------------------------------------------------------------
spicyTestBootstrap

top(spicyTestBootstrap)  

signifPlot(spicyTestBootstrap, breaks=c(-3, 3, 0.5))

