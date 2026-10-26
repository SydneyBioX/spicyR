# Package index

## Test

Every pair of cell types, between conditions or with survival.

- [`SpicyResults-class`](https://sydneybiox.github.io/spicyR/dev/reference/SpicyResults-class.md)
  [`SpicyResults`](https://sydneybiox.github.io/spicyR/dev/reference/SpicyResults-class.md)
  : The results of spicy()
- [`bind()`](https://sydneybiox.github.io/spicyR/dev/reference/bind.md)
  : The per-image values of every pair, or of one pair, as a data frame
- [`calibrationCheck()`](https://sydneybiox.github.io/spicyR/dev/reference/calibrationCheck.md)
  [`summary(`*`<spicyCalibration>`*`)`](https://sydneybiox.github.io/spicyR/dev/reference/calibrationCheck.md)
  [`plot(`*`<spicyCalibration>`*`)`](https://sydneybiox.github.io/spicyR/dev/reference/calibrationCheck.md)
  : Check that spicy() is calibrated on your data by shuffling the
  condition
- [`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
  : Test for changes in the co-localisation of cell types between
  conditions
- [`topPairs()`](https://sydneybiox.github.io/spicyR/dev/reference/topPairs.md)
  : The most significant pairs, or the results of chosen pairs

## Plot

- [`imageCrossPlot()`](https://sydneybiox.github.io/spicyR/dev/reference/imageCrossPlot.md)
  : Pairwise cross-plot of spatial associations

- [`plotImage()`](https://sydneybiox.github.io/spicyR/dev/reference/plotImage.md)
  :

  Plot images, showing their cell types or the `from` and `to` cells of
  a pair

- [`radiusPlot()`](https://sydneybiox.github.io/spicyR/dev/reference/radiusPlot.md)
  : The evidence for pairs at each radius

- [`signifPlot()`](https://sydneybiox.github.io/spicyR/dev/reference/signifPlot.md)
  : Plots result of signifPlot.

- [`spicyBoxPlot()`](https://sydneybiox.github.io/spicyR/dev/reference/spicyBoxPlot.md)
  : Box plot of one pair, with a point per image

## Per-image statistics (image method)

Inputs and helpers for `spicy(method = "image")`, the original spicyR
test.

- [`buildKnnGraph()`](https://sydneybiox.github.io/spicyR/dev/reference/buildKnnGraph.md)
  : Build a k-nearest-neighbour graph between cells
- [`colTest()`](https://sydneybiox.github.io/spicyR/dev/reference/colTest.md)
  : Perform a simple wilcoxon-rank-sum test or t-test on the columns of
  a data frame
- [`convPairs()`](https://sydneybiox.github.io/spicyR/dev/reference/convPairs.md)
  : Converts colPairs object into an abundance matrix based on number of
  nearby interactions for every cell type.
- [`getPairwise()`](https://sydneybiox.github.io/spicyR/dev/reference/getPairwise.md)
  : Get statistic from pairwise L curve of a single image.
- [`getPairwiseProp()`](https://sydneybiox.github.io/spicyR/dev/reference/getPairwiseProp.md)
  : Observed/expected proportion ratio for a fixed-k neighbourhood
- [`getProp()`](https://sydneybiox.github.io/spicyR/dev/reference/getProp.md)
  : Get proportions from a SummarizedExperiment.

## Large data sets

Cells stored on disk as a SpatialData Zarr store, which
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
reads image by image.

- [`writeCellStore()`](https://sydneybiox.github.io/spicyR/dev/reference/writeCellStore.md)
  : Write cells as a SpatialData Zarr store

## Data

- [`diabetesData`](https://sydneybiox.github.io/spicyR/dev/reference/diabetesData.md)
  : Imaging mass cytometry of human pancreas in type 1 diabetes (Damond
  et al. 2019)
- [`spicyTest`](https://sydneybiox.github.io/spicyR/dev/reference/spicyTest.md)
  : Results of the image-level test on diabetesData
