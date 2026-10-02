# Results of the image-level test on diabetesData

The result of
`spicy(diabetesData, condition = "stage", subject = "case", method = "image")`:
the original image-level test of spicyR, for all pairs of cell types,
comparing the onset and long-duration stages with non-diabetic donors.
It is used in examples, so that they run quickly.

## Usage

``` r
data("spicyTest")
```

## Format

A `SpicyResults` object (see
[`spicy`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)).

## Source

Made by `inst/scripts/make-spicyTest.R` from
[`diabetesData`](https://sydneybiox.github.io/spicyR/dev/reference/diabetesData.md).

## Examples

``` r
data("spicyTest")
topPairs(spicyTest, n = 5)
#>                  intercept coefficient      p.value adj.pvalue  from    to
#> beta__delta   6.081500e+01   -15.43325 0.0006616394 0.09171144  beta delta
#> delta__beta   6.090722e+01   -15.27636 0.0007164956 0.09171144 delta  beta
#> B__Th        -4.440892e-16    10.48012 0.0127535323 0.42316680     B    Th
#> delta__delta  7.021912e+01   -16.35809 0.0155005358 0.42316680 delta delta
#> Th__B         2.220446e-15    10.00345 0.0173028313 0.42316680    Th     B
```
