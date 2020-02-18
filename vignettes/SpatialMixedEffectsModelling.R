## ----warning=F, message=F-----------------------------------------------------
# load required packages
library(spicyR)
library(SegmentedCellExperiment)
library(ggplot2)

## ----echo=F-------------------------------------------------------------------
cellFile <- system.file("extdata","responderData.rds", 
                        package = "spicyR")
cells <- readRDS(cellFile)

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
mixed.lmer <- spatialREM(cells,
                         pairwise = pairwiseVals, 
                         from = "CD8+PD1+PDL1+", 
                         to = "CD8+PD1+PDL1+")

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

