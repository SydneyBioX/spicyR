# The per-image values of every pair, as a data frame

Returns, for each image, its condition and the value of each pair that
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
tested: the excess (extra `from` cells per `to` cell beyond chance) for
the cell method, or the L-function summary for the image method. Use it
for your own plots or models.

## Usage

``` r
bind(results, pairName = NULL)
```

## Arguments

- results:

  The `SpicyResults` object returned by
  [`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md).

- pairName:

  A pair, as `"from__to"`. If NULL, all pairs are returned.

## Value

A data.frame with one row per image: `imageID`, `condition`, `subject`
(if given) and one column per pair, named `"from__to"`. `NA` where the
pair could not be measured in an image.

## Examples

``` r

data(spicyTest)
df <- bind(spicyTest)
```
