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

diabetesData <- buildKnnGraph(diabetesData, k = 20)

pairAbundances <- convPairs(diabetesData,
  colPair = "knn_interaction_graph"
)
```
