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
  combine = c("maxT", "davies", "cauchy"),
  adjustAbundance = FALSE,
  variance = c("cr2", "hartung_knapp"),
  ref = NULL,
  cores = 1,
  tableKey = NULL,
  imageInfo = NULL,
  cellWeights = NULL,
  ...
)
```

## Arguments

- cells:

  A `data.frame`, `SingleCellExperiment` or `SpatialExperiment` with one
  row (column) per cell, or (cell method) the path of a SpatialData or
  AnnData Zarr store (read with the Rarr package). From a store, the
  columns are read a block of rows at a time and the coordinates a batch
  of whole images at a time, so the table is never held in memory; the
  results are those of the same cells given as a `data.frame`. Column
  names then refer to the table's `obs`; `imageID` defaults to the
  table's region key and `spatialCoords` to its `obsm` array
  `"spatial"`. A `SpatialData` object (spatialdataR) read from a store
  is read from that store in the same way; one built in memory is
  analysed from its table.

- condition:

  The column of the image-level condition: two groups, or a
  [`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html) column
  for a survival outcome (from a store, a column of `imageInfo`).

- subject:

  The column of the patient (unit) of each image. Images of one patient
  are combined; if omitted, every image is its own patient.

- covariates:

  Image- or patient-level columns to adjust for (cell method: added to
  the design of the excess, and the effect of each is reported as
  `<column>_effect` and `<column>_p_value`; survival: added to the null
  Cox model).

- imageID, cellType, spatialCoords:

  Column names of the image, cell type and coordinates. For a Zarr
  store, `spatialCoords` is either two `obs` columns or the name of one
  `obsm` array (its first two columns).

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
  correlation across radii; up to 6 radii), `"davies"` (Davies' bound
  for the same maximum, for many radii) or `"cauchy"` (Cauchy
  combination). The max-T p-value is the multivariate normal probability
  of the largest \|z\|, computed by the quasi-Monte Carlo method of Genz
  and Bretz with an absolute error of about 1e-4 (by quadrature for two
  radii); the same inputs always give the same value.

- adjustAbundance:

  Cell method: adjust the test for the log share of the `from` type in
  each image (default `FALSE`). Its effect is reported as
  `abundance_effect`. It does not separate more `from` cells from more
  densely packed ones, and it removes real effects when the share tracks
  the condition.

- variance:

  Cell method: `"cr2"` (the default: CR2 on Satterthwaite df, with the
  dispersion shared by the conditions, see Details) or `"hartung_knapp"`
  (the model-based variance, scaled by the Pearson statistic when that
  exceeds its df and never below CR2, on m - 2 df).

- ref:

  Cell method: the reference level of `condition`.

- cores:

  Number of threads (cell method) or cores (image method).

- tableKey:

  Cell method with a SpatialData store or object: the name of the cell
  table in `tables/` (by default the only table, or the one called
  `"table"`).

- imageInfo:

  Cell method: a table with one row per image, for `condition`,
  `subject`, `covariates` or a survival outcome that are not columns of
  the cell table: a `data.frame` whose `imageID` column, or else its row
  names, gives the image IDs, or the name of a table of the SpatialData
  store or object. Every image of the cells needs a row, and every row
  an image. A
  [`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html) column
  given as `condition` makes a survival analysis, also from a store.

- cellWeights:

  Optional cell weights in \[0, 1\], in the order of `cells`: a numeric
  vector with one weight per cell for its own cell type (for example a
  quality score), or a matrix (or sparse `Matrix`) of soft cell type
  membership with one row per cell and one named column per label (for
  example a classifier's confidence). With a matrix, `from` and `to`
  name its columns, and labels may overlap (a group and its subtypes).
  Cell method: a `to` cell counts with its weight, and its score is the
  chance that it has a `from` neighbour within `r` if each cell were a
  `from` cell with probability its `from` weight (allocation), or the
  summed `from` weight of its neighbours (count). The null shuffles the
  `to` weights over the cells whose `from` weight is below 1 (so a cell
  of the `from` type with weight 0.6 is a candidate with the remaining
  0.4), and its exact mean and variance replace those of random
  labelling; weights of 1 give the results of the cell types. The
  abundance share is the `from` label's summed weight over the image's
  cells. One radius and cells in memory; pairs of a label with itself
  are not tested. Image method: see
  [`getPairwise()`](https://sydneybiox.github.io/spicyR/dev/reference/getPairwise.md).

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
BH-adjusted p-value, the frailty variance, the between-image variance
(`tau2_image`) and, when the test was adjusted, what for
(`adjusted_for`), the effect and p-value of each adjustment, and the
unadjusted test (`unadjusted_difference`, `unadjusted_p_value`,
`unadjusted_p_adj`). `$effect` records which effect was estimated.

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
packed. The side is chosen once per pair from all images with a
condition, without the conditions, and is reported in the `side` column.
`effect = "count"` gives the number of extra `from` cells within `r` of
each `to` cell instead; it also reflects how many `from` cells surround
a `to` cell (depth of infiltration), but it grows with how densely the
`from` cells are packed, so a change in packing alone can appear as a
change in co-localisation. Images are combined within patients and
patients within conditions by a frailty GEE with a second random effect
for the images of a patient (each image is weighted by 1 / (its sampling
variance + the between-image variance), so one large image cannot
outweigh the patient's other images), and the difference between
conditions is tested with a CR2 cluster-robust variance on Satterthwaite
degrees of freedom, with **patients (`subject`) as the units**. The
difference is adjusted for any `covariates`, and, with
`adjustAbundance = TRUE`, for the log share of the `from` type in each
image; the unadjusted test is then reported alongside (`unadjusted_*`
columns).

When nearly every cell has a `from` cell within `r` (q close to 1),
there is little room for attraction and an attracting pair's images
carry little information; a smaller `r` is more informative.

The allocation effect is bounded: an image of an avoiding pair with no
`to` cell next to a `from` cell has the effect -1, and an image of an
attracting pair with every `to` cell next to one has +1. A patient is
informative at a radius when at least one of its images is not at -1 or
+1. A radius where one condition has no informative patient is not
tested: that condition's patients all have the same effect, so its
spread is exactly 0 and the standard error would rest on the other
condition alone. Its p-value is `NA` (with adjustment, `adjusted_for` is
`"none (no_information)"`), and it is left out of the combination over
radii. With more than two conditions, a pair where any condition has no
informative patient is not tested.

The conditions share one dispersion. Each condition's CR2 variance is
divided by its model-based variance, the ratios are averaged over the
conditions with their numbers of patients as weights, and the
model-based variance of the difference is multiplied by this average.
The degrees of freedom are those of CR2. A small condition's own spread
does not set the standard error alone, which matters when the conditions
are of unequal size and the patients' effects are skewed. With
adjustment, the ratios are those of the coefficient of each condition in
the adjusted design, and every coefficient's variance is scaled the same
way. `variance = "hartung_knapp"` keeps its own variance.

A patient's weight within its condition grows with the information of
its images, so one or two patients can carry nearly all of a condition's
weight. A radius where a condition's effective patients (Kish's (sum
w)^2 / sum w^2 of the patients' weights w) fall below the smaller of 3
and half the number of its patients that contribute at that radius is
not tested: its p-value is `NA` (with adjustment, `adjusted_for` is
`"none (few_effective_units)"`), and it is left out of the combination
over radii. A condition of 3 patients with equal weights can therefore
be tested. With more than two conditions, a pair where any condition
falls below this is not tested.

Images with a missing condition are left out of every step: the choice
of the side, the between-image variance and the tests.

**`method = "image"` (the original spicyR test).** A per-image
L-function summary of each pair is compared between conditions with a
weighted linear model, or a mixed model when `subject` is given (Canete
et al. 2022). With `cellWeights` (see
[`getPairwise()`](https://sydneybiox.github.io/spicyR/dev/reference/getPairwise.md)),
cells enter each label with a weight, for example the confidence of a
classifier, instead of through a single cell type.

## References

Canete NP et al. (2022). spicyR: spatial analysis of in situ cytometry
data in R. Bioinformatics 38(11), 3099-3105.
[doi:10.1093/bioinformatics/btac268](https://doi.org/10.1093/bioinformatics/btac268)

Bell RM, McCaffrey DF (2002). Bias reduction in standard errors for
linear regression with multi-stage samples. Survey Methodology 28(2),
169-181.

DerSimonian R, Laird N (1986). Meta-analysis in clinical trials.
Controlled Clinical Trials 7(3), 177-188.
[doi:10.1016/0197-2456(86)90046-2](https://doi.org/10.1016/0197-2456%2886%2990046-2)

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
#> Tc__Th  0.272768   0.1023755 0.4114593  0.4114593   Tc Th
res$cellResults
#>                       from to         level    side excess_ref
#> Tc__Th__Onset           Tc Th         Onset attract   0.272768
#> Tc__Th__Long-duration   Tc Th Long-duration attract   0.272768
#>                       excess_difference        se       df   p_value     p_adj
#> Tc__Th__Onset                 0.1023755 0.1146014 5.153215 0.4114593 0.4114593
#> Tc__Th__Long-duration         0.0398330 0.1117464 5.231229 0.7354259 0.7354259
#>                              tau2 tau2_image
#> Tc__Th__Onset         0.001709649 0.05857091
#> Tc__Th__Long-duration 0.001709649 0.05857091

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
