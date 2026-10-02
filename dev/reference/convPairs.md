# Converts colPairs object into an abundance matrix based on number of nearby interactions for every cell type.

Converts colPairs object into an abundance matrix based on number of
nearby interactions for every cell type.

## Usage

``` r
convPairs(cells, colPair, imageID = "imageID", cellType = "cellType")
```

## Arguments

- cells:

  A SingleCellExperiment that contains objects in the colPairs slot.

- colPair:

  The name of the object in the colPairs slot for which the dataframe is
  constructed from.

- imageID:

  The image ID if using SingleCellExperiment.

- cellType:

  The cell type if using SingleCellExperiment.

## Value

Matrix of abundances

## Examples

``` r
data("diabetesData")
images <- c("A09", "A11", "A16", "A17")
diabetesData <- diabetesData[
  , SummarizedExperiment::colData(diabetesData)$imageID %in% images
]

diabetesData_SPE <- SpatialExperiment::SpatialExperiment(diabetesData,
  colData = SummarizedExperiment::colData(diabetesData)
)
SpatialExperiment::spatialCoords(diabetesData_SPE) <- data.frame(
  SummarizedExperiment::colData(diabetesData_SPE)$x,
  SummarizedExperiment::colData(diabetesData_SPE)$y
) |>
  as.matrix()

SpatialExperiment::spatialCoordsNames(diabetesData_SPE) <- c("x", "y")

diabetesData_SPE <- imcRtools::buildSpatialGraph(diabetesData_SPE,
  img_id = "imageID",
  type = "knn",
  k = 20,
  coords = c("x", "y")
)
#> 'sample_id's are duplicated across 'SpatialExperiment' objects to cbind; appending sample indices.
#> The returned object is ordered by the 'imageID' entry.

pairAbundances <- convPairs(diabetesData_SPE,
  colPair = "knn_interaction_graph"
)
```
