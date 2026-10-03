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
  k = NULL,
  combine = c("maxT", "cauchy"),
  adjustAbundance = TRUE,
  variance = c("cr2", "hartung_knapp"),
  frailty = TRUE,
  labelClustering = TRUE,
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

- k:

  Cell method: use the `k` nearest neighbours instead of a radius.

- combine:

  Cell method with several radii: `"maxT"` (max-T with the sandwich
  correlation across radii) or `"cauchy"` (Cauchy combination).

- adjustAbundance:

  Cell method: adjust the test for the log share of the `from` type in
  each image (default `TRUE`). Its effect is reported as
  `abundance_effect`. `FALSE` gives the test without it.

- variance:

  Cell method: `"cr2"` (CR2 on Satterthwaite df, the default) or
  `"hartung_knapp"` (for very few patients: the model-based variance
  floored at CR2, on m - 2 df).

- frailty, labelClustering:

  Cell method: the patient frailty and the label-clustering inflation of
  the within-image variance (both on by default). The inflation of a
  `to` type is estimated from every counted type, so a pair's result
  does not depend on which other pairs are requested.

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
full table: the excess in each condition (at the average abundance and
covariates), the difference, its standard error, df, p-value and
BH-adjusted p-value, the frailty variance, what the test was adjusted
for (`adjusted_for`), the effect and p-value of each adjustment, and the
unadjusted test (`unadjusted_difference`, `unadjusted_p_value`,
`unadjusted_p_adj`).

## Details

**`method = "cell"` (the default).** For each `to` cell, the number of
`from` cells within radius `r` is compared with its exact expectation if
the `to` cells were a random choice among the cells that are not `from`
cells in the same image. The effect is the **excess**: the number of
extra `from` cells within `r` of each `to` cell. Images are combined
within patients and patients within conditions by a frailty GEE, and the
difference between conditions is tested with a CR2 cluster-robust
variance on Satterthwaite degrees of freedom, with **patients
(`subject`) as the units**. By default the difference is adjusted for
how common the `from` type is in each image (the log of its share of all
cells) and for any `covariates`, so that a change in abundance alone
does not appear as a change in co-localisation. The unadjusted test is
reported alongside (`unadjusted_*` columns).

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
# extra Tc cells within 50 units of each Th cell and of each beta cell
res <- spicy(diabetesData, condition = "stage", subject = "case", r = 50,
             from = "Tc", to = c("Th", "beta"))
topPairs(res)
#>        intercept coefficient  p.value adj.pvalue from to
#> Tc__Th 0.6375863 -0.09109796 0.679257   0.679257   Tc Th
res$cellResults
#>                       from to         level excess_ref excess_difference
#> Tc__Th__Onset           Tc Th         Onset  0.6375863      -0.091097963
#> Tc__Th__Long-duration   Tc Th Long-duration  0.6375863      -0.008162035
#>                              se       df   p_value     p_adj       tau2
#> Tc__Th__Onset         0.2089707 5.573440 0.6792570 0.6792570 0.06581695
#> Tc__Th__Long-duration 0.2297366 5.839256 0.9728421 0.9728421 0.06581695
#>                       adjusted_for abundance_effect abundance_p_value
#> Tc__Th__Onset            abundance         0.476337       0.002122868
#> Tc__Th__Long-duration    abundance         0.476337       0.002122868
#>                       unadjusted_difference unadjusted_se unadjusted_df
#> Tc__Th__Onset                     0.3860306     0.2935793      5.614328
#> Tc__Th__Long-duration             0.1746664     0.0847683      5.775973
#>                       unadjusted_p_value unadjusted_p_adj
#> Tc__Th__Onset                 0.23969003       0.23969003
#> Tc__Th__Long-duration         0.08681395       0.08681395

# the original image-level test
resImage <- spicy(diabetesData, condition = "stage", subject = "case",
                  from = "Tc", to = "Th", method = "image")
#> Dropping unused levels. Using stage = Non-diabetic as base comparison group. If this is not the desired base group, please convert cells$stage into a factor and change the order of levels(cells$stage) so that the base group is at index 1.
topPairs(resImage)
#>        intercept coefficient   p.value adj.pvalue from to
#> Tc__Th  1.622671    5.961812 0.6122508  0.6122508   Tc Th
```
