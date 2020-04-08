## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(BiocStyle)

## ----setup, message=FALSE-----------------------------------------------------
library(spicyR)
library(S4Vectors)

## -----------------------------------------------------------------------------

### Something that resembles cellProfiler data

set.seed(51773)

n = 10

cells <- data.frame(row.names = seq_len(n))
cells$ObjectNumber <- seq_len(n)
cells$ImageNumber <- rep(1:2,c(n/2,n/2))
cells$AreaShape_Center_X <- runif(n)
cells$AreaShape_Center_Y <- runif(n)
cells$AreaShape_round <- rexp(n)
cells$AreaShape_diameter <- rexp(n, 2)
cells$Intensity_Mean_CD8 <- rexp(n, 10)
cells$Intensity_Mean_CD4 <- rexp(n, 10)


## -----------------------------------------------------------------------------
cellExp <- segmentedCells(cells, cellProfiler = TRUE)
cellExp

## -----------------------------------------------------------------------------
loc <- location(cellExp)
head(loc)

location(cellExp) <- loc

## -----------------------------------------------------------------------------
intensities <- intensity(cellExp)
kM <- kmeans(intensities,2)
cellType(cellExp) <- paste('cluster',kM$cluster, sep = '')

loc <- location(cellExp)
head(loc)

## -----------------------------------------------------------------------------
isletFile <- system.file("extdata","isletCells.txt.gz", package = "spicyR")
cells <- read.table(isletFile, header = TRUE)

## -----------------------------------------------------------------------------
cellExp <- segmentedCells(cells, cellProfiler = TRUE)
cellExp

## -----------------------------------------------------------------------------
intensities <- intensity(cellExp)
kM <- kmeans(intensities,4)
cellType(cellExp) <- paste('cluster',kM$cluster, sep = '')

loc <- location(cellExp)
head(loc)

## ---- fig.width=5, fig.height= 6----------------------------------------------
plot(cellExp, imageID=1)

## -----------------------------------------------------------------------------
set.seed(51773)

n = 10

cells <- data.frame(row.names = seq_len(n))
cells$cellID <- seq_len(n)
cells$imageCellID <- rep(seq_len(n/2),2)
cells$imageID <- rep(1:2,c(n/2,n/2))
cells$x <- runif(n)
cells$y <- runif(n)
cells$shape_round <- rexp(n)
cells$shape_diameter <- rexp(n, 2)
cells$intensity_CD8 <- rexp(n, 10)
cells$intensity_CD4 <- rexp(n, 10)
cells$cellType <- paste('cluster',sample(1:2,n,replace = TRUE), sep = '_')


## -----------------------------------------------------------------------------

cellExp <- segmentedCells(cells, cellTypeString = 'cellType', intensityString = 'intensity_', morphologyString = 'shape_')
cellExp


## -----------------------------------------------------------------------------
morph <- morphology(cellExp)
head(morph)


## -----------------------------------------------------------------------------
phenoData <- DataFrame(imageID = c('1','2'), 
                       age = c(21,81), 
                       status = c('dead','alive'))
phenotype(cellExp) <- phenoData
phenotype(cellExp)
phenotype(cellExp, expand = TRUE)

## -----------------------------------------------------------------------------
set.seed(51773)

n = 10

cells <- data.frame(row.names = seq_len(n))
cells$x <- runif(n)
cells$y <- runif(n)
cellExp <- segmentedCells(cells)
cellExp


## -----------------------------------------------------------------------------
loc <- location(cellExp)
head(loc)


