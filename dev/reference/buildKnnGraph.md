# Build a k-nearest-neighbour graph between cells

Finds the `k` nearest neighbours of every cell within its own image and
stores them as a `colPair` of `cells`, ready for
[`convPairs()`](https://sydneybiox.github.io/spicyR/dev/reference/convPairs.md).
Distances are Euclidean, a cell is not its own neighbour, and ties are
broken as in
[`spatstat.geom::nnwhich()`](https://rdrr.io/pkg/spatstat.geom/man/nnwhich.html).
The graph is the one `imcRtools::buildSpatialGraph(type = "knn")`
builds, computed by the C++ core that
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
uses for `k`.

## Usage

``` r
buildKnnGraph(
  cells,
  k = 20,
  imageID = "imageID",
  spatialCoords = c("x", "y"),
  name = "knn_interaction_graph",
  cores = 1
)
```

## Arguments

- cells:

  A `SingleCellExperiment` or `SpatialExperiment`.

- k:

  The number of neighbours of each cell. Cells in images with at most
  `k` cells get no neighbours.

- imageID:

  The `colData` column holding the image of each cell.

- spatialCoords:

  The `colData` columns holding the coordinates. A `SpatialExperiment`'s
  own `spatialCoords()` are used unless these name `colData` columns.

- name:

  The name of the `colPair` the graph is stored in.

- cores:

  The number of threads; images are processed in parallel.

## Value

`cells`, with the graph in `colPair(cells, name)`: one edge from each
cell to each of its neighbours.

## Examples

``` r
data("diabetesData")
diabetesData <- buildKnnGraph(diabetesData, k = 10)
head(SingleCellExperiment::colPair(diabetesData, "knn_interaction_graph"))
#> SelfHits object with 6 hits and 0 metadata columns:
#>            from        to
#>       <integer> <integer>
#>   [1]         1         6
#>   [2]         1        11
#>   [3]         1        12
#>   [4]         1        20
#>   [5]         1        59
#>   [6]         1        66
#>   -------
#>   nnode: 253777
```
