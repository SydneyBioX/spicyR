# The most significant pairs, or the results of chosen pairs

Lists the pairs tested by
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md),
most significant first. Give `from` and/or `to` to see chosen pairs, and
`n = NULL` to see every pair.

## Usage

``` r
topPairs(
  x,
  coef = NULL,
  n = 10,
  adj = "fdr",
  cutoff = NULL,
  figures = NULL,
  from = NULL,
  to = NULL,
  full = FALSE
)
```

## Arguments

- x:

  The output from
  [`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md).

- coef:

  Which coefficient to list: a column of `x$p.value`, or with the cell
  method and more than two conditions, a level of the condition.

- n:

  The number of pairs to list. `NULL` lists every pair, as does leaving
  `n` out when `cutoff` is given.

- adj:

  Which p-value adjustment method to use, argument for p.adjust().

- cutoff:

  List only the pairs whose adjusted p-value is at most `cutoff`.

- figures:

  Round to `figures` significant figures.

- from, to:

  Only the pairs with these `from` and `to` cell types.

- full:

  Cell method: every column of the test, not only the main ones.

## Value

A data.frame, one row per pair, with the pair (`"from__to"`) as row
names.

## Details

For the cell method (the default of
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md))
the table has one row per pair (and per level with more than two
conditions), with the columns of the test: `from`, `to`, `excess_ref`
and `excess_comp` (the average effect in the reference and the
comparison group), `excess_difference`, `p_value` and `p_adj`; with
several radii, `r`, the radius with the strongest evidence; with
covariates, the effect and p-value of each (`<covariate>_effect`,
`<covariate>_p_value`); for survival, `hazard_ratio_sd`. `full = TRUE`
adds the scale of the effect (`side`), the standard error, degrees of
freedom, the variance between patients (`tau2`), what the test was
adjusted for and the unadjusted test. For the image method the columns
are `intercept`, `coefficient`, `p.value`, `adj.pvalue`, `from` and
`to`, as in spicyR 1.x.

## Examples

``` r

data(spicyTest)
topPairs(spicyTest)
#>                          intercept coefficient      p.value adj.pvalue
#> beta__delta           6.081500e+01  -15.433247 0.0006616394 0.09171144
#> delta__beta           6.090722e+01  -15.276362 0.0007164956 0.09171144
#> B__Th                -4.440892e-16   10.480124 0.0127535323 0.42316680
#> delta__delta          7.021912e+01  -16.358091 0.0155005358 0.42316680
#> Th__B                 2.220446e-15   10.003453 0.0173028313 0.42316680
#> B__unknown            7.771561e-16    4.584274 0.0182158621 0.42316680
#> otherimmune__naiveTc -3.179087e+00   11.944639 0.0199646993 0.42316680
#> unknown__macrophage   4.339587e+00   -5.274674 0.0222619394 0.42316680
#> unknown__B           -1.443290e-15    4.680750 0.0244004277 0.42316680
#> macrophage__unknown   4.305429e+00   -4.886438 0.0249693477 0.42316680
#>                             from         to
#> beta__delta                 beta      delta
#> delta__beta                delta       beta
#> B__Th                          B         Th
#> delta__delta               delta      delta
#> Th__B                         Th          B
#> B__unknown                     B    unknown
#> otherimmune__naiveTc otherimmune    naiveTc
#> unknown__macrophage      unknown macrophage
#> unknown__B               unknown          B
#> macrophage__unknown   macrophage    unknown

data("diabetesData")
res <- spicy(diabetesData, condition = "stage", subject = "case", r = 50, from = "Tc", to = c("Th", "beta"))
topPairs(res, from = "Tc", to = "Th")
#>                       from to         level excess_ref excess_difference
#> Tc__Th__Onset           Tc Th         Onset  0.3199387       0.072364141
#> Tc__Th__Long-duration   Tc Th Long-duration  0.3199387       0.001750948
#>                         p_value     p_adj
#> Tc__Th__Onset         0.5909207 0.5909207
#> Tc__Th__Long-duration 0.9852766 0.9852766
```
