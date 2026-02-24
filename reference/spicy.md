# Performs spatial tests on spatial cytometry data.

Performs spatial tests on spatial cytometry data.

## Usage

``` r
spicy(
  cells,
  condition,
  subject = NULL,
  covariates = NULL,
  imageID = "imageID",
  cellType = "cellType",
  spatialCoords = c("x", "y"),
  r = NULL,
  sigma = NULL,
  from = NULL,
  to = NULL,
  alternateResult = NULL,
  cores = 1,
  minLambda = 0.05,
  weights = TRUE,
  weightsByPair = FALSE,
  weightFactor = 1,
  window = "convex",
  window.length = NULL,
  edgeCorrect = TRUE,
  includeZeroCells = FALSE,
  verbose = FALSE,
  BPPARAM = NULL,
  imageIDCol = imageID,
  cellTypeCol = cellType,
  spatialCoordCols = spatialCoords,
  nCores = cores,
  Rs = r,
  ...
)
```

## Arguments

- cells:

  A SummarizedExperiment or data frame that contains at least the
  variables x and y, giving the location coordinates of each cell, and
  cellType.

- condition:

  A character specifying which column which contains the condition or
  \`Surv\` objects.

- subject:

  Vector of subject IDs corresponding to each image if cells is a data
  frame.

- covariates:

  Vector of covariate names that should be included in the mixed effects
  model as fixed effects.

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

  vector of cell types which you would like to compare to the to vector.

- to:

  vector of cell types which you would like to compare to the from
  vector.

- alternateResult:

  A pairwise association statistic between each combination of celltypes
  in each image.

- cores:

  Number of cores to use for parallel processing or a BiocParallel
  MulticoreParam or SerialParam object.

- minLambda:

  Minimum value density for scaling when fitting inhomogeneous L-curves.

- weights:

  logical indicating whether to include weights based on cell counts.

- weightsByPair:

  logical indicating whether weights should be calculated for each cell
  type pair.

- weightFactor:

  numeric that controls the convexity of the weight function.

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
  pairwise association calculation.

- verbose:

  logical indicating whether to output messages.

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

- ...:

  Other options

## Value

Data frame of p-values.

## Examples

``` r
data("diabetesData")

# Test with random effect for patient on a pairwise combination of cell
# types.
spicy(diabetesData,
  condition = "stage", subject = "case",
  from = "Tc", to = "Th"
)
#> Dropping unused levels. Using stage = Non-diabetic as base comparison group. If
#> this is not the desired base group, please convert cells$stage into a factor
#> and change the order of levels(cells$stage) so that the base group is at index
#> 1.
#> 
#> Number of cell type pairs: 1
#> Number of differentially localised cell type pairs: 
#> [1] 0

# Test all pairwise combinations of cell types without random effect of
# patient.
if (FALSE) { # \dontrun{
spicyTest <- spicy(diabetesData, condition = "stage", subject = "case")
} # }

# Test all pairwise combination of cell types with random effect of patient.
if (FALSE) { # \dontrun{
spicy(diabetesData, condition = "condition", subject = "subject")
} # }
```
