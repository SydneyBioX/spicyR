# Engine of Statial's Kontextual test

The computation behind `Statial::kontextualTest()`; use that function,
whose help page describes the method. For each triple of `parentDf` it
compares Statial's Kontextual statistic in each image with its exact
expectation and variance when the `to` cells are a random choice among
the parent's cells, and combines the images as
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
does: each image's effect (observed minus expected, over Kontextual's
normaliser) with the image random effect within patients (Paule-Mandel),
then the frailty GEE across patients with each condition's CR2 variance
and Satterthwaite df, and enough effective patients per condition. For a
triple with `from == to` the statistic is quadratic in the labels, and
its variance is the exact variance of observed minus expected under that
random choice; an image where every choice gives the same value carries
no information and is left out.

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
  adjustAbundance = FALSE,
  variance = c("cr2", "hartung_knapp"),
  frailty = TRUE,
  labelClustering = FALSE,
  edgeCorrect = TRUE,
  window = c("convex", "rectangle"),
  ref = NULL,
  cellWeights = NULL
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
  variance, frailty, ref:

  As for
  [`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md);
  `r` is one radius (several radii are not supported).

- adjustAbundance:

  Adjust the test for the log share of the `to` type among the parent's
  cells in each image (default `FALSE`), as
  [`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)'s
  `adjustAbundance`; its effect is reported as `abundance_effect`.

- labelClustering:

  Whether to multiply each image's variance by a factor for clustering
  of the `to` labels. Off by default: the factor puts most of the weight
  on a few images, and without it the test stays at or below its nominal
  level when the conditions are shuffled between patients.

- edgeCorrect:

  Correct the parent densities for the part of each disc outside the
  image's window.

- window:

  The window of each image for the edge correction: `"convex"` or
  `"rectangle"`.

- cellWeights:

  Optional cell weights, a vector or a soft membership matrix as for
  [`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md):
  with a matrix, `from` and `to` (and the `from`, `to` of `parentDf`)
  name its columns, while `parent` still lists cell types. A `from` cell
  enters with its weight, the score of a candidate counts with its `to`
  weight, and the null shuffles the `to` weights over the candidates.
  Pairs of a label with itself are not tested, and `labelClustering` is
  not available.

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
