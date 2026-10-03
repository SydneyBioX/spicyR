# Engine of Statial's Kontextual test

The computation behind `Statial::kontextualTest()`; use that function,
whose help page describes the method. For each triple of `parentDf` it
compares Statial's Kontextual statistic in each image with its exact
expectation and variance when the `to` cells are a random choice among
the parent's cells, and combines the images as
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
does (frailty GEE, CR2 variance, label-clustering factor, abundance
adjustment).

## Usage

``` r
kontextualEngine(
  cells,
  parentDf,
  condition,
  subject = NULL,
  covariates = NULL,
  r = 50,
  from = NULL,
  to = NULL,
  imageID = "imageID",
  cellType = "cellType",
  spatialCoords = c("x", "y"),
  adjustAbundance = TRUE,
  variance = c("cr2", "hartung_knapp"),
  frailty = TRUE,
  labelClustering = TRUE,
  edgeCorrect = TRUE,
  window = c("convex", "rectangle"),
  ref = NULL
)
```

## Arguments

- cells:

  A SingleCellExperiment, SpatialExperiment or data frame with a row per
  cell.

- parentDf:

  A data frame of triples with columns `from`, `to`, `parent` (a list
  column of cell types containing `to`) and optionally `parent_name`, as
  made by `Statial::parentCombinations()`.

- condition:

  The column with each image's condition, or a `Surv` column of survival
  outcomes.

- subject, covariates, r, from, to, imageID, cellType, spatialCoords,
  adjustAbundance, variance, frailty, labelClustering, ref:

  As for
  [`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md);
  `r` is one radius.

- edgeCorrect:

  Correct the parent densities for the part of each disc outside the
  image's window.

- window:

  The window of each image for the edge correction: `"convex"` or
  `"rectangle"`.

## Value

A
[SpicyResults](https://sydneybiox.github.io/spicyR/dev/reference/SpicyResults-class.md)
object with one row per triple (`from`, `to`, `parent`).

## Examples

``` r
data("diabetesData")
parentDf <- data.frame(from = "alpha", to = c("Tc", "Th"), parent_name = "tcells")
parentDf$parent <- list(c("Tc", "Th"), c("Tc", "Th"))
res <- kontextualEngine(diabetesData, parentDf, condition = "stage", subject = "case", r = 50)
```
