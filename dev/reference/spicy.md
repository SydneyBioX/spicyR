# Test for changes in the co-localisation of cell types between conditions

`spicy()` tests, for every ordered pair of cell types `from` → `to`,
whether the co-localisation of the two types differs between conditions,
or is associated with survival. A pair asks whether `to` cells are
placed near `from` cells more than other cells are.

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
  from = NULL,
  to = NULL,
  method = c("cell", "image"),
  effect = c("allocation", "count"),
  k = NULL,
  combine = c("maxT", "cauchy"),
  adjustAbundance = FALSE,
  variance = c("cr2", "hartung_knapp", "auto"),
  frailty = TRUE,
  labelClustering = FALSE,
  ref = NULL,
  cores = 1,
  ...
)
```

## Arguments

- cells:

  A `data.frame`, `SingleCellExperiment` or `SpatialExperiment` with one
  row (column) per cell.

- condition:

  The column of the image-level condition: two groups, or a
  [`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html) column
  for a survival outcome.

- subject:

  The column of the patient (unit) of each image. Images of one patient
  are combined; if omitted, every image is its own patient.

- covariates:

  Image- or patient-level columns to adjust for (cell method: added to
  the design of the excess, and the effect of each is reported as
  `<column>_effect` and `<column>_p_value`; survival: added to the null
  Cox model).

- imageID, cellType, spatialCoords:

  Column names of the image, cell type and coordinates.

- r:

  Radius (or radii) of the neighbourhood, in the units of the
  coordinates. Cell method: one radius (default 50), or several to be
  combined by `combine`. Image method: the radii of the L function
  (default 20, 50 and 100).

- from, to:

  Cell types to test (all ordered pairs by default).

- method:

  `"cell"` (spicyR Cell, the default) or `"image"` (the original spicyR
  test).

- effect:

  Cell method: `"allocation"` (the default), the extra fraction of `to`
  cells with at least one `from` cell within `r`; or `"count"`, the
  number of extra `from` cells within `r` of each `to` cell. With `k`,
  "within `r`" means among the cell's `k` nearest neighbours.

- k:

  Cell method: use the `k` nearest neighbours instead of a radius.

- combine:

  Cell method with several radii: `"maxT"` (max-T with the sandwich
  correlation across radii) or `"cauchy"` (Cauchy combination).

- adjustAbundance:

  Cell method: adjust the test for the log share of the `from` type in
  each image (default `FALSE`). Its effect is reported as
  `abundance_effect`. It does not separate more `from` cells from more
  densely packed ones, and it removes real effects when the share tracks
  the condition.

- variance:

  Cell method: `"cr2"` (the default: CR2 on Satterthwaite df),
  `"hartung_knapp"` (for very few patients: the model-based variance
  floored at CR2, on m - 2 df) or `"auto"` (`"hartung_knapp"` when a
  condition has at most 5 patients, `"cr2"` otherwise). The variance
  used is in `$variance`.

- frailty, labelClustering:

  Cell method: the patient frailty (on by default) and the
  label-clustering inflation of the within-image variance (off by
  default). The inflation of a `to` type is estimated from every counted
  type, so a pair's result does not depend on which other pairs are
  requested.

- ref:

  Cell method: the reference level of `condition`.

- cores:

  Number of threads (cell method) or cores (image method).

- ...:

  Arguments of the image method: `sigma`, `alternateResult`,
  `minLambda`, `weights`, `weightsByPair`, `weightFactor`,
  `weightZThreshold`, `window`, `window.length`, `edgeCorrect`,
  `includeZeroCells`, `verbose`, `BPPARAM`. Supplying `alternateResult`
  selects the image method.

## Value

A `SpicyResults` object.
[`topPairs()`](https://sydneybiox.github.io/spicyR/dev/reference/topPairs.md),
[`signifPlot()`](https://sydneybiox.github.io/spicyR/dev/reference/signifPlot.md),
[`spicyBoxPlot()`](https://sydneybiox.github.io/spicyR/dev/reference/spicyBoxPlot.md)
and
[`bind()`](https://sydneybiox.github.io/spicyR/dev/reference/bind.md)
work for both methods. For the cell method, `$cellResults` holds the
full table: for the allocation effect, the pair's `side` (`"attract"` or
`"avoid"`, which sets the scale of the effect), the effect in each
condition (`excess_ref`, `excess_comp`, at the average covariates), the
difference (`excess_difference`), its standard error, df, p-value and
BH-adjusted p-value, the frailty variance and, when the test was
adjusted, what for (`adjusted_for`), the effect and p-value of each
adjustment, and the unadjusted test (`unadjusted_difference`,
`unadjusted_p_value`, `unadjusted_p_adj`). `$effect` records which
effect was estimated.

## Details

**`method = "cell"` (the default).** In each image, the share of `to`
cells with at least one `from` cell within radius `r` is compared with
its exact expectation q if the `to` cells were a random choice among the
cells that are not `from` cells. The effect (`effect = "allocation"`,
the default) is the **fraction of `to` cells placed next to (or kept
away from) `from` cells**. For a pair that attracts (more `to` cells
next to `from` cells than q over all images together) it is (observed
share - q) / (1 - q): if a fraction f of the `to` cells were moved next
to `from` cells, the effect is f. For a pair that avoids it is (observed
share - q) / q: if a fraction f of the `to` cells that would have a
`from` cell nearby were moved away, the effect is -f. Either way it does
not depend on how many `from` cells there are or how densely they are
packed. The side is chosen once per pair from all images, without the
conditions, and is reported in the `side` column. `effect = "count"`
gives the number of extra `from` cells within `r` of each `to` cell
instead; it also reflects how many `from` cells surround a `to` cell
(depth of infiltration), but it grows with how densely the `from` cells
are packed, so a change in packing alone can appear as a change in
co-localisation. Images are combined within patients and patients within
conditions by a frailty GEE, and the difference between conditions is
tested with a CR2 cluster-robust variance on Satterthwaite degrees of
freedom, with **patients (`subject`) as the units**. The difference is
adjusted for any `covariates`, and, with `adjustAbundance = TRUE`, for
the log share of the `from` type in each image; the unadjusted test is
then reported alongside (`unadjusted_*` columns).

When nearly every cell has a `from` cell within `r` (q close to 1),
there is little room for attraction and an attracting pair's images
carry little information; a smaller `r` is more informative.

**`method = "image"` (the original spicyR test).** A per-image
L-function summary of each pair is compared between conditions with a
weighted linear model, or a mixed model when `subject` is given (Canete
et al. 2022).

## References

Canete NP et al. (2022). spicyR: spatial analysis of in situ cytometry
data in R. Bioinformatics 38(11), 3099-3105.
[doi:10.1093/bioinformatics/btac268](https://doi.org/10.1093/bioinformatics/btac268)

Bell RM, McCaffrey DF (2002). Bias reduction in standard errors for
linear regression with multi-stage samples. Survey Methodology 28(2),
169-181.

Pustejovsky JE, Tipton E (2018). Small-sample methods for cluster-robust
variance estimation and hypothesis testing in fixed effects models.
Journal of Business & Economic Statistics 36(4), 672-683.
[doi:10.1080/07350015.2016.1247004](https://doi.org/10.1080/07350015.2016.1247004)

Paule RC, Mandel J (1982). Consensus values and weighting factors.
Journal of Research of the National Bureau of Standards 87(5), 377-385.
[doi:10.6028/jres.087.022](https://doi.org/10.6028/jres.087.022)

Hartung J, Knapp G (2001). A refined method for the meta-analysis of
controlled clinical trials with binary outcome. Statistics in Medicine
20(24), 3875-3889.
[doi:10.1002/sim.1009](https://doi.org/10.1002/sim.1009)

Liu Y, Xie J (2020). Cauchy combination test: a powerful test with
analytic p-value calculation under arbitrary dependency structures.
Journal of the American Statistical Association 115(529), 393-402.
[doi:10.1080/01621459.2018.1554485](https://doi.org/10.1080/01621459.2018.1554485)

## Examples

``` r
data("diabetesData")
# spicyR Cell: patients ("case") are the units
# the extra fraction of Th cells, and of beta cells, with a Tc cell within 50 units
res <- spicy(diabetesData, condition = "stage", subject = "case", r = 50,
             from = "Tc", to = c("Th", "beta"))
topPairs(res)
#>        intercept coefficient   p.value adj.pvalue from to
#> Tc__Th 0.2573739   0.1042285 0.3792698  0.3792698   Tc Th
res$cellResults
#>                       from to         level    side excess_ref
#> Tc__Th__Onset           Tc Th         Onset attract  0.2573739
#> Tc__Th__Long-duration   Tc Th Long-duration attract  0.2573739
#>                       excess_difference         se       df   p_value     p_adj
#> Tc__Th__Onset                0.10422855 0.10860300 5.251571 0.3792698 0.3792698
#> Tc__Th__Long-duration        0.05161495 0.07179114 5.290259 0.5026656 0.5026656
#>                              tau2
#> Tc__Th__Onset         0.009930696
#> Tc__Th__Long-duration 0.009930696

# the extra number of Tc cells within 50 units of each Th cell
resCount <- spicy(diabetesData, condition = "stage", subject = "case", r = 50,
                  from = "Tc", to = "Th", effect = "count")

# the original image-level test
resImage <- spicy(diabetesData, condition = "stage", subject = "case",
                  from = "Tc", to = "Th", method = "image")
#> Dropping unused levels. Using stage = Non-diabetic as base comparison group. If this is not the desired base group, please convert cells$stage into a factor and change the order of levels(cells$stage) so that the base group is at index 1.
topPairs(resImage)
#>        intercept coefficient   p.value adj.pvalue from to
#> Tc__Th  1.622672    5.961811 0.6122509  0.6122509   Tc Th
```
