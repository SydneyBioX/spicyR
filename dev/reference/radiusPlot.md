# The evidence for pairs at each radius

For a
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
result with several radii, plots the p-value of each chosen pair at
every radius, to show the scale at which the evidence is strongest. The
combined p-value over the radii is in
[`topPairs()`](https://sydneybiox.github.io/spicyR/dev/reference/topPairs.md).

## Usage

``` r
radiusPlot(results, from, to)
```

## Arguments

- results:

  The `SpicyResults` object returned by
  [`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
  with several radii (`r`).

- from, to:

  The `from` and `to` cell types of the pairs to show: every pair
  between them.

## Value

A ggplot.

## Examples

``` r
data("diabetesData")
res <- spicy(diabetesData, condition = "stage", subject = "case", r = c(25, 50, 75),
             from = c("Tc", "Th"), to = "beta")
#> Error: several radii are supported for two conditions; give one `r`.
radiusPlot(res, from = c("Tc", "Th"), to = "beta")
#> Error: object 'res' not found
```
