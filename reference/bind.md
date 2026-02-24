# Produces a dataframe showing L-function metric for each imageID entry.

Produces a dataframe showing L-function metric for each imageID entry.

## Usage

``` r
bind(results, pairName = NULL)
```

## Arguments

- results:

  Spicy test result obtained from spicy.

- pairName:

  A string specifying the pairwise interaction of interest. If NULL, all
  pairwise interactions are shown.

## Value

A data.frame containing the colData related to the results.

## Examples

``` r
data(spicyTest)
df <- bind(spicyTest)
```
