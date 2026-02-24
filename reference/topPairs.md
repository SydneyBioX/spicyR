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

  Round to \`figures\` significant figures.

## Value

A data.frame

## Examples

``` r
data(spicyTest)
topPairs(spicyTest)
#>                          intercept coefficient      p.value adj.pvalue
#> beta__delta           1.815458e+02  -48.735693 0.0005033247 0.07169649
#> delta__beta           1.817943e+02  -48.166076 0.0005601288 0.07169649
#> B__unknown           -3.606296e-15   11.770938 0.0052338392 0.42051606
#> delta__delta          2.089550e+02  -52.061196 0.0125129422 0.42051606
#> unknown__macrophage   1.007337e+01  -15.826919 0.0207410908 0.42051606
#> unknown__B           -5.161932e-15   12.142848 0.0225855419 0.42051606
#> macrophage__unknown   1.004424e+01  -14.471666 0.0244668075 0.42051606
#> B__Th                 1.571479e-15   26.847934 0.0245039854 0.42051606
#> otherimmune__naiveTc -9.292508e+00   33.584755 0.0255812944 0.42051606
#> ductal__ductal        1.481580e+01   -8.632569 0.0266935703 0.42051606
#>                             from         to
#> beta__delta                 beta      delta
#> delta__beta                delta       beta
#> B__unknown                     B    unknown
#> delta__delta               delta      delta
#> unknown__macrophage      unknown macrophage
#> unknown__B               unknown          B
#> macrophage__unknown   macrophage    unknown
#> B__Th                          B         Th
#> otherimmune__naiveTc otherimmune    naiveTc
#> ductal__ductal            ductal     ductal
```
