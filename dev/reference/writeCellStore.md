# Write cells as a SpatialData Zarr store

Writes a `SpatialExperiment`, `SingleCellExperiment` or `data.frame` of
cells as a SpatialData Zarr store (Zarr v3), which
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
reads a block of rows at a time, spatialdataR's `readSpatialData()`
reads as a `SpatialData` object, and Python's `spatialdata` reads too.
The cells become a table (AnnData, under `tables/`):

## Usage

``` r
writeCellStore(
  cells,
  path,
  imageID = "imageID",
  spatialCoords = c("x", "y"),
  tableKey = "table",
  imageInfo = NULL,
  infoKey = "images",
  points = TRUE
)
```

## Arguments

- cells:

  A `SpatialExperiment`, `SingleCellExperiment` or `data.frame` with one
  row (column) per cell. The `colData`, `rowData` and `data.frame`
  columns are factors, character, integer, double or logical vectors,
  and every assay has a name.

- path:

  The directory of the store (replaced if it exists).

- imageID, spatialCoords:

  The column of the image and the two columns of the coordinates. A
  `SpatialExperiment`'s own `spatialCoords()` are used unless
  `spatialCoords` names two `colData` columns.

- tableKey:

  The name of the cell table.

- imageInfo:

  Optional: a `data.frame` with one row per image, written as the table
  `infoKey`. The image IDs are its `imageID` column, or else its row
  names.

- infoKey:

  The name of the per-image table.

- points:

  Whether to write the cells as points too, one points element per image
  (needs duckdb).

## Value

`path`, invisibly.

## Details

- `obs`: the `colData` (or the columns of a `data.frame`) other than the
  coordinates, factors as categorical columns. The image column is the
  table's region key (always categorical, its categories the levels of a
  factor, else the images in order of first appearance), and the cells
  are numbered 0, 1, ... in `"instance_id"`.

- `X` and `layers`: the first assay and the others (a sparse assay in
  compressed sparse row form).

- `var`: the `rowData`, with the row names as its index.

- `obsm`: the coordinates as `"spatial"`, and the `reducedDims`.

The cells of each image also become a points element named after the
image, with their coordinates and their `"instance_id"`, so that the
table annotates elements of the store (image IDs must then be valid
element names: letters, digits, `_`, `-` and `.`). The images of a
`SpatialExperiment` (`imgData()`) become image elements named
`<sample_id>_<image_id>`, each in the coordinate system of its sample
and scaled to the coordinates of the cells. A table with one row per
image can be added. Every array is written in chunks (of 2^20 rows for a
column) compressed with zstd, one chunk at a time, so writing needs
little memory beyond that of `cells`, however many cells it has.
`metadata()`, `altExps()` and column or row pairings are not written.

## Examples

``` r
cells <- data.frame(imageID = factor(rep(c("a", "b"), each = 50)), cellType = factor(rep(c("T", "B"), 50)),
                    x = runif(100, 0, 100), y = runif(100, 0, 100))
if (requireNamespace("Rarr", quietly = TRUE)) {
  path <- writeCellStore(cells, tempfile(fileext = ".zarr"))
  list.files(file.path(path, "tables"))
}
#> [1] "table"     "zarr.json"
```
