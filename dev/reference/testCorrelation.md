# The correlation between the tests of a spicy() result

The tests of pairs that share a cell type, or that are nested (a cell
type and a group containing it), are correlated through the patients
they share. `testCorrelation()` estimates the correlation between every
two pairs' tests from the one fit, for a correction for multiple testing
that uses it in place of independence.

## Usage

``` r
testCorrelation(
  x,
  method = c("cr2", "permutation"),
  nPermutations = 1000,
  cores = 1
)
```

## Arguments

- x:

  The output of
  [`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
  (or
  [`kontextualEngine()`](https://sydneybiox.github.io/spicyR/dev/reference/kontextualEngine.md))
  with two conditions, one radius (or `k`), no adjustment (no
  `covariates`, `adjustAbundance = FALSE`) and `variance = "cr2"`.

- method:

  `"cr2"` (the default) or `"permutation"`.

- nPermutations:

  The number of permutations for `method = "permutation"`.

- cores:

  The number of threads for `method = "permutation"`.

## Value

A symmetric matrix with one row and one column per row of
`x$cellResults` (the same row names): the estimated correlation between
the pairs' tests. The diagonal is 1, and rows and columns of pairs not
tested (or, by permutation, whose z is the same in every permutation)
are `NA`.

## Details

Each pair's difference between conditions is a weighted mean over
patients: a patient's summary is the weighted mean of its images'
effects (`x$pairwiseAssoc`) with the images' weights (`x$imageWeights`),
and its weight is the sum of its images' weights times its condition's
information (`x$conditionInformation`).

`method = "cr2"`: each patient's CR2 influence on the difference, psi =
h (summary - its condition's mean) / sqrt(1 - h) with h the patient's
share of its condition's weight, negated in the reference condition. Its
squares sum to the test's squared standard error, and the cross-products
of two pairs' influences, divided by the two standard errors, estimate
the correlation of their differences. It is noisy when a few patients
carry most of a pair's weight and one condition is small.

`method = "permutation"`: the conditions are permuted over the patients
`nPermutations` times; for each permutation every test is recomputed in
closed form from the patients' summaries and weights (each condition's
weighted mean, its CR2 variance and Satterthwaite df, as
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
computes them; the between-patient variance is held at the fit's, and of
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)'s
rules for not testing a pair only the need for two patients in each
condition is applied), and its p-value is turned into a z with the sign
of the difference (0 when a condition has fewer than two of the pair's
patients). The correlation of the tests' z over the permutations is
returned. The permutations come from a fixed seed, so a result always
gives the same matrix, in R and in Python.

## See also

[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)

## Examples

``` r
data("diabetesData")
cells <- diabetesData[, diabetesData$stage != "Long-duration"]
res <- spicy(cells, condition = "stage", subject = "case", r = 50, from = c("Tc", "Th"), to = c("beta", "alpha"))
testCorrelation(res)
#>            Tc__beta Tc__alpha  Th__beta Th__alpha
#> Tc__beta  1.0000000 0.9336402 0.7045028 0.7501011
#> Tc__alpha 0.9336402 1.0000000 0.8968579 0.8537341
#> Th__beta  0.7045028 0.8968579 1.0000000 0.8805660
#> Th__alpha 0.7501011 0.8537341 0.8805660 1.0000000
testCorrelation(res, method = "permutation", nPermutations = 200)
#>            Tc__beta Tc__alpha  Th__beta Th__alpha
#> Tc__beta  1.0000000 0.9214592 0.7267367 0.7895594
#> Tc__alpha 0.9214592 1.0000000 0.8768567 0.8923144
#> Th__beta  0.7267367 0.8768567 1.0000000 0.8700752
#> Th__alpha 0.7895594 0.8923144 0.8700752 1.0000000
```
