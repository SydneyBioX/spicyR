# Observed/expected proportion ratio for a fixed-k neighbourhood

Computes a per-image, per-cell-type-pair spatial localisation statistic
as a drop-in alternative to the fixed-radius L-function statistic
returned by
[`getPairwise`](https://sydneybiox.github.io/spicyR/dev/reference/getPairwise.md).
Instead of a search radius `r`, each reference cell's neighbourhood is
its `k` nearest neighbours (Euclidean distance, any cell type, the cell
itself excluded).

## Usage

``` r
getPairwiseProp(
  cells,
  imageID = "imageID",
  cellType = "cellType",
  spatialCoords = c("x", "y"),
  k = 15,
  from = NULL,
  to = NULL,
  cores = 1,
  includeZeroCells = FALSE,
  BPPARAM = NULL
)
```

## Arguments

- cells:

  A SingleCellExperiment, SpatialExperiment or data.frame.

- imageID:

  The name of the imageID column if a column name is provided.

- cellType:

  The name of the cellType column if a column name is provided.

- spatialCoords:

  The names of the spatial coordinates columns if provided.

- k:

  Number of nearest neighbours to search per reference cell.

- from:

  The reference cell types. Defaults to all cell types.

- to:

  The target cell types. Defaults to all cell types.

- cores:

  Number of threads for processing images in parallel (a
  BiocParallelParam object is also accepted, for backward
  compatibility).

- includeZeroCells:

  If FALSE (default), image-pairs where the reference or target cell
  type is absent are returned as NA. If TRUE, image-pairs whose
  reference type is absent (but whose target type is present, so
  \\p\_{0,j} \> 0\\) are floored at \\\hat\pi_j = 0\\, i.e. \\v_j = 0\\.
  Pairs whose target type is entirely absent (\\p\_{0,j} = 0\\) stay NA
  regardless, since the ratio is undefined.

- BPPARAM:

  A BiocParallelParam object; its number of workers overrides `cores`.
  Kept for backward compatibility.

## Value

A matrix with one row per image and one column per `from__to` cell-type
pair, containing the observed/expected proportion ratio.

## Details

For an ordered pair (`from` = reference type A, `to` = target type B),
the statistic for image \\j\\ is the plain observed-over-expected ratio
\$\$v_j = \hat\pi_j / p\_{0,j},\$\$ where \\\hat\pi_j\\ is the mean over
reference cells of the fraction of their `k` nearest neighbours that are
type B, and \\p\_{0,j}\\ is the image's background proportion of type-B
cells (B count / total cell count).

\\v_j = 1\\ is the null (complete spatial randomness); \\v_j \> 1\\
indicates attraction/enrichment and \\0 \le v_j \< 1\\ depletion. No
variance-stabilising transform (arcsine, log) is applied: \\v_j\\ is
bounded below at 0 and right-skewed, a rougher fit to the linear model's
Gaussian-noise assumption than a transformed statistic, but chosen for
direct interpretability ("observed is 1.4x expected").

The returned matrix is shaped identically to
[`getPairwise()`](https://sydneybiox.github.io/spicyR/dev/reference/getPairwise.md)'s
output (images in rows, `from__to` cell-type pairs in columns, same
order), so it can be passed straight to
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
as `alternateResult` to run the standard spicyR weighting and linear
(mixed) model pipeline on the ratio.

**Weights.**
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)'s
variance-weighting model (`weights = TRUE`, the default) is calibrated
to the L-function's numeric scale. It works for this ratio when the
dataset has enough cross-pair variance heterogeneity (e.g. several cell
types, some rare) to anchor the fit, but on a dataset with few cell
types or weak spatial structure it can degenerate and return an all-`NA`
p-value table. If that happens, refit with `weights = FALSE`.

No edge correction is applied: \\\hat\pi_j\\ estimates neighbourhood
composition rather than a count within an area, so a boundary cell's
truncated neighbourhood carries the same null label distribution as an
interior cell's. `k` controls only the variance of \\\hat\pi_j\\, not
the location of its null value.

## Examples

``` r
data("diabetesData")
propAssoc <- getPairwiseProp(diabetesData, k = 15)

# Use as the spicyR response, reusing the standard model:
# \donttest{
spicy(diabetesData,
  condition = "stage", subject = "case",
  alternateResult = propAssoc
)
#> `alternateResult` / `sigma` given: using method = "image" (the original spicyR test).
#> Dropping unused levels. Using stage = Non-diabetic as base comparison group. If this is not the desired base group, please convert cells$stage into a factor and change the order of levels(cells$stage) so that the base group is at index 1.
#> spicyR (image-level test): 256 pairs
#> Pairs with BH-adjusted p < 0.05:
#>         conditionOnset conditionLong-duration 
#>                      0                     13 
# }
```
