## ----message=FALSE------------------------------------------------------------
# load required packages
library(spicyR)
library(tidyverse)


## ----eval=T-------------------------------------------------------------------
set.seed(51773)
x <- round(c(runif(200),runif(200)+1,runif(200)+2,runif(200)+3,
           runif(200)+3,runif(200)+2,runif(200)+1,runif(200)),4)
y <- round(c(runif(200),runif(200)+1,runif(200)+2,runif(200)+3,
             runif(200),runif(200)+1,runif(200)+2,runif(200)+3),4)
cellType <- factor(paste('c',rep(rep(c(1:2),rep(200,2)),4),sep = ''))
imageID <- rep(c('s1', 's2'),c(800,800))

cells <- data.frame(x, y, cellType, imageID)

ggplot(cells, aes(x,y, colour = cellType)) + geom_point() + facet_wrap(~imageID)



## -----------------------------------------------------------------------------

cellExp <- SegmentedCellExperiment(cells, cellTypeString = 'cellType')



## -----------------------------------------------------------------------------

lisaCurves <- lisa(cellExp)


## -----------------------------------------------------------------------------

kM <- kmeans(lisaCurves,2)
region(cellExp) <- paste('region',kM$cluster,sep = '_')

## -----------------------------------------------------------------------------
hatchingPlot(cellExp, image = c('s1','s2'))

## -----------------------------------------------------------------------------

df <- region(cellExp, annot = TRUE)

p <- ggplot(df,aes(x = x,y = y, colour = cellType, region = region)) + 
  geom_point() + 
  facet_wrap(~imageID) +
  geom_hatching(window = "concave", line.spacing = 11, nbp = 50, line.width = 1.5, window.length = 0.1) +
  theme_minimal()+ 
  scale_region_manual(values = 6:7, labels = c('ab','cd'))

p


