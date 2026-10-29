# The per-image values of every pair, or of one pair, as a data frame

Returns, for each image, its condition and the value of each pair that
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
tested: for the cell method its effect (by default the extra fraction of
`to` cells with a `from` cell within `r`; with `effect = "count"`, the
extra `from` cells per `to` cell), or the L-function summary for the
image method. Use it for your own plots or models, or with `from` and
`to` to look up the images behind one pair.

## Usage

``` r
bind(results, pairName = NULL, from = NULL, to = NULL, parent = NULL)
```

## Arguments

- results:

  The `SpicyResults` object returned by
  [`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md).

- pairName:

  A pair, as `"from__to"`. If NULL, all pairs are returned.

- from, to:

  One pair: a table of that pair alone, with its value in a column
  `effect` (`value` for the image method) and, for the cell method, each
  image's `weight` in the test.

- parent:

  For a Kontextual result, the parent population of the pair given by
  `from` and `to`.

## Value

A data.frame with one row per image, named by image: `imageID`,
`condition`, `subject` (if given) and one column per pair, named
`"from__to"`, or with `from` and `to` the columns `effect` (or `value`)
and `weight`. `NA` where the pair could not be measured in an image.

## Examples

``` r

data(spicyTest)
df <- bind(spicyTest)

data("diabetesData")
res <- spicy(diabetesData, condition = "stage", subject = "case", r = 50, from = "Tc", to = "Th")
head(bind(res, from = "Tc", to = "Th"))
#>     imageID condition subject    effect weight
#> A09     A09     Onset    6362        NA     NA
#> A11     A11     Onset    6362        NA     NA
#> A16     A16     Onset    6362 1.0000000     NA
#> A17     A17     Onset    6362 0.2318507     NA
#> A25     A25     Onset    6362 0.1641414     NA
#> A39     A39     Onset    6362        NA     NA
```
