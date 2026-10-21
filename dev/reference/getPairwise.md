# Get statistic from pairwise L curve of a single image.

Get statistic from pairwise L curve of a single image.

## Usage

``` r
getPairwise(
  cells,
  imageID = "imageID",
  cellType = "cellType",
  spatialCoords = c("x", "y"),
  r = NULL,
  sigma = NULL,
  from = NULL,
  to = NULL,
  cores = 1,
  minLambda = 0.05,
  window = "convex",
  window.length = NULL,
  edgeCorrect = TRUE,
  includeZeroCells = FALSE,
  BPPARAM = NULL,
  imageIDCol = imageID,
  cellTypeCol = cellType,
  spatialCoordCols = spatialCoords,
  nCores = cores,
  Rs = r,
  cellWeights = NULL
)
```

## Arguments

- cells:

  A SummarizedExperiment that contains at least the variables x and y,
  giving the location coordinates of each cell, and cellType.

- imageID:

  The name of the imageID column if using a SingleCellExperiment or
  SpatialExperiment.

- cellType:

  The name of the cellType column if using a SingleCellExperiment or
  SpatialExperiment.

- spatialCoords:

  The names of the spatialCoords column if using a SingleCellExperiment.

- r:

  A vector of the radii that the measures of association should be
  calculated over.

- sigma:

  A numeric variable used for scaling when fitting inhomogenous
  L-curves.

- from:

  The 'from' cellType for generating the L curve.

- to:

  The 'to' cellType for generating the L curve.

- cores:

  Number of cores to use for parallel processing or a BiocParallel
  MulticoreParam or SerialParam object.

- minLambda:

  Minimum value density for scaling when fitting inhomogeneous L-curves.

- window:

  Should the window around the regions be 'square', 'convex' or
  'concave'.

- window.length:

  A tuning parameter for controlling the level of concavity when
  estimating concave windows.

- edgeCorrect:

  A logical indicating whether to perform edge correction.

- includeZeroCells:

  A logical indicating whether to include cells with zero counts in the
  pairwise association.

- BPPARAM:

  {DEPRECATED} A BiocParallel MulticoreParam or SerialParam object.

- imageIDCol:

  {DEPRECATED} The name of the imageID column if using a
  SingleCellExperiment or SpatialExperiment.

- cellTypeCol:

  {DEPRECATED} The name of the cellType column if using a
  SingleCellExperiment or SpatialExperiment.

- spatialCoordCols:

  {DEPRECATED} The names of the spatialCoords column if using a
  SingleCellExperiment.

- nCores:

  {DEPRECATED} Number of cores to use for parallel processing or a
  BiocParallel MulticoreParam or SerialParam object.

- Rs:

  {DEPRECATED} A vector of the radii that the measures of association
  should be calculated over.

- cellWeights:

  Optional cell weights in \\0, 1\\, in the order of `cells`: a numeric
  vector with one weight per cell for its own cell type (for example a
  quality score), or a numeric matrix (or sparse `Matrix`) of soft cell
  type membership with one row per cell and one named column per label,
  holding each cell's weight for that label (for example a classifier's
  confidence). With a matrix, `from` and `to` name its columns, and
  labels may overlap, for example a group and its subtypes. Each cell
  contributes to a label in proportion to its weight; a label's density
  is its summed weight per unit area. Weights of 1 give the values of
  the cell types. Needs a `"square"` or `"convex"` window and no
  `sigma`.

## Value

Statistic from pairwise L-curve of a single image.

## Examples

``` r
data("diabetesData")
# Subset by imageID for fast example
selected_cells <- diabetesData[
  , SummarizedExperiment::colData(diabetesData)$imageID == "A09"
]
pairAssoc <- getPairwise(selected_cells)
```
