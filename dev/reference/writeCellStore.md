# Write cells as a SpatialData Zarr store

Writes a `data.frame` of cells as the cell table of a SpatialData Zarr
store (an AnnData table under `tables/`, written with anndataR), which
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
reads a block of rows at a time and spatialdataR's `readSpatialData()`
reads as a `SpatialData` object. The coordinates become the `obsm` array
`"spatial"`, the other columns the table's `obs` (factors as categorical
columns), and the image column the table's region key. A table with one
row per image can be added as a second table.

## Usage

``` r
writeCellStore(
  cells,
  path,
  imageID = "imageID",
  spatialCoords = c("x", "y"),
  tableKey = "table",
  imageInfo = NULL,
  infoKey = "images"
)
```

## Arguments

- cells:

  A `data.frame` with one row per cell.

- path:

  The directory of the store (replaced if it exists).

- imageID, spatialCoords:

  The column of the image and the two columns of the coordinates.

- tableKey:

  The name of the cell table.

- imageInfo:

  Optional: a `data.frame` with one row per image (row names: the image
  IDs), written as the table `infoKey`.

- infoKey:

  The name of the per-image table.

## Value

`path`, invisibly.

## Examples

``` r
cells <- data.frame(imageID = factor(rep(c("a", "b"), each = 50)), cellType = factor(rep(c("T", "B"), 50)),
                    x = runif(100, 0, 100), y = runif(100, 0, 100))
if (requireNamespace("anndataR", quietly = TRUE) && requireNamespace("Rarr", quietly = TRUE)) {
  path <- writeCellStore(cells, tempfile(fileext = ".zarr"))
  list.files(file.path(path, "tables"))
}
#> [1] "table"     "zarr.json"
```
