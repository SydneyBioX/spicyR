# Check that spicy() is calibrated on your data by shuffling the condition

Shuffles the condition across patients (or across images, without
`subject`), which removes any real difference between the groups while
keeping everything else about the study, and runs
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
on each shuffle. If the test is calibrated, about 5% of pairs have p \<
0.05 and about 1% p \< 0.01 in each shuffle. Pairs share cells, so the
percentage moves by a few points from one shuffle to the next; look at
the average over shuffles.

## Usage

``` r
calibrationCheck(
  cells,
  condition,
  subject = NULL,
  nShuffles = 10,
  thresholds = c(0.05, 0.01),
  ...
)

# S3 method for class 'spicyCalibration'
summary(object, ...)

# S3 method for class 'spicyCalibration'
plot(x, ...)
```

## Arguments

- cells:

  The cells, as for
  [`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md):
  a `SpatialExperiment`, `SingleCellExperiment` or `data.frame` (not a
  store).

- condition:

  The column of the condition to shuffle (a factor, character or `Surv`
  column).

- subject:

  The column of the patient. Each patient keeps one condition, moved as
  a whole; without it, the condition is shuffled across images.

- nShuffles:

  The number of shuffles.

- thresholds:

  The p-value thresholds to report.

- ...:

  Other arguments of
  [`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md),
  such as `r`, `imageID`, `cellType` or `covariates`.

- object, x:

  The result of `calibrationCheck()`.

## Value

A `data.frame` with one row per shuffle and threshold: `shuffle`,
`threshold` and `percent`, the percentage of pairs whose p-value is
below the threshold. [`summary()`](https://rdrr.io/r/base/summary.html)
gives the average percentage over the shuffles at each threshold, with
the expected one, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows each
shuffle against the expected percentage.

## Examples

``` r
data("diabetesData")
set.seed(1)
shuffles <- calibrationCheck(diabetesData, condition = "stage", subject = "case", r = 50, nShuffles = 2)
summary(shuffles)
#>   threshold expected   percent
#> 1      0.05        5 0.7989084
#> 2      0.01        1 0.0000000
plot(shuffles)
```
