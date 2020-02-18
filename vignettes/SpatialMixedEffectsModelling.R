## ----warning=F, message=F-----------------------------------------------------
# load required packages
library(spicyR)
library(SegmentedCellExperiment)
library(ggplot2)

## ----echo=F-------------------------------------------------------------------
cells <- readRDS("/Users/nick/Desktop/cellExpTest2.rds")

## -----------------------------------------------------------------------------
head(location(cells))
head(phenotype(cells))

## -----------------------------------------------------------------------------
pairwiseVals <- getPairwise(cells, from = "CD8+PD1+PDL1+", to = "CD8+PD1+PDL1+")

## -----------------------------------------------------------------------------
spatialData <- data.frame(Pairwise = pairwiseVals,
                          Condition = phenotype(cells)$condition)

ggplot(spatialData, aes(x=Condition, y=Pairwise, color = Condition)) +
    geom_boxplot() +
    theme_classic()

## -----------------------------------------------------------------------------
subject <- phenotype(cells)$subject
condition <- phenotype(cells)$condition

cellSplit <- location(cells, bind = FALSE)
count <- unlist(lapply(cellSplit, 
                       function(x) sum(x$cellType == "CD8+PD1+PDL1+")))
count <- as.numeric(count)

## -----------------------------------------------------------------------------
mixed.lmer <- spatialREM(pairwiseVals, 
                         subject = subject, 
                         condition = condition, 
                         count1 = count, 
                         count2 = count)

summary(mixed.lmer)

## -----------------------------------------------------------------------------
set.seed(51773)
pVal <- spatialREMBootstrap(mixed.lmer, nsim = 999)
pVal

## -----------------------------------------------------------------------------
df <- spatialREMMulti(cells, nsim = NULL)
head(df)

## -----------------------------------------------------------------------------
spatialREMMultiPlot(df,
                    fdr = FALSE,
                    breaks = c(-4,4,0.5),
                    col=c("blue", "white", "red"))

