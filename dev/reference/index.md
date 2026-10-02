# Package index

## Test

Every pair of cell types, between conditions or with survival.

- [`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
  : Test for changes in the co-localisation of cell types between
  conditions
- [`SpicyResults-class`](https://sydneybiox.github.io/spicyR/dev/reference/SpicyResults-class.md)
  [`SpicyResults`](https://sydneybiox.github.io/spicyR/dev/reference/SpicyResults-class.md)
  : The results of spicy()
- [`topPairs()`](https://sydneybiox.github.io/spicyR/dev/reference/topPairs.md)
  : A table of the significant results from spicy tests
- [`bind()`](https://sydneybiox.github.io/spicyR/dev/reference/bind.md)
  : The per-image values of every pair, as a data frame

## Plot

- [`signifPlot()`](https://sydneybiox.github.io/spicyR/dev/reference/signifPlot.md)
  : Plots result of signifPlot.
- [`spicyBoxPlot()`](https://sydneybiox.github.io/spicyR/dev/reference/spicyBoxPlot.md)
  : Box plot of one pair, with a point per image
- [`plotImage()`](https://sydneybiox.github.io/spicyR/dev/reference/plotImage.md)
  : Plot one image, showing the \`from\` and \`to\` cells of a pair
- [`imageCrossPlot()`](https://sydneybiox.github.io/spicyR/dev/reference/imageCrossPlot.md)
  : Pairwise cross-plot of spatial associations

## Per-image statistics (image method)

Inputs and helpers for `spicy(method = "image")`, the original spicyR
test.

- [`getPairwise()`](https://sydneybiox.github.io/spicyR/dev/reference/getPairwise.md)
  : Get statistic from pairwise L curve of a single image.
- [`getPairwiseProp()`](https://sydneybiox.github.io/spicyR/dev/reference/getPairwiseProp.md)
  : Observed/expected proportion ratio for a fixed-k neighbourhood
- [`getProp()`](https://sydneybiox.github.io/spicyR/dev/reference/getProp.md)
  : Get proportions from a SummarizedExperiment.
- [`convPairs()`](https://sydneybiox.github.io/spicyR/dev/reference/convPairs.md)
  : Converts colPairs object into an abundance matrix based on number of
  nearby interactions for every cell type.
- [`colTest()`](https://sydneybiox.github.io/spicyR/dev/reference/colTest.md)
  : Perform a simple wilcoxon-rank-sum test or t-test on the columns of
  a data frame

## Data

- [`diabetesData`](https://sydneybiox.github.io/spicyR/dev/reference/diabetesData.md)
  : Imaging mass cytometry of human pancreas in type 1 diabetes (Damond
  et al. 2019)
- [`spicyTest`](https://sydneybiox.github.io/spicyR/dev/reference/spicyTest.md)
  : Results of the image-level test on diabetesData
