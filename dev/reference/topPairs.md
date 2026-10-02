# A table of the significant results from spicy tests

A table of the significant results from spicy tests

## Usage

``` r
topPairs(x, coef = NULL, n = 10, adj = "fdr", cutoff = NULL, figures = NULL)
```

## Arguments

- x:

  The output from spicy.

- coef:

  Which coefficient to list.

- n:

  Extract the top n most significant pairs.

- adj:

  Which p-value adjustment method to use, argument for p.adjust().

- cutoff:

  A p-value threshold to extract significant pairs.

- figures:

  Round to `figures` significant figures.

## Value

A data.frame

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
```
