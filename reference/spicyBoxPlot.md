# Plots boxplot for a specified cell-cell relationship

Plots boxplot for a specified cell-cell relationship

## Usage

``` r
spicyBoxPlot(results, from = NULL, to = NULL, rank = NULL)
```

## Arguments

- results:

  Data frame obtained from spicy.

- from:

  Cell type which you would like to compare to the to cell type.

- to:

  Cell type which you would like to compare to the from cell type.

- rank:

  Ranking of cell type in terms of p-value, the smaller the p-value the
  higher the rank.

## Value

a ggplot2 boxplot

## Examples

``` r
data(spicyTest)

spicyBoxPlot(spicyTest,
             rank = 1)
#> Warning: Removed 53 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).

```
