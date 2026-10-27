# Plots result of signifPlot.

Plots result of signifPlot.

## Usage

``` r
signifPlot(
  results,
  fdr = FALSE,
  type = "bubble",
  breaks = NULL,
  comparisonGroup = NULL,
  colours = c("#4575B4", "white", "#D73027"),
  marksToPlot = NULL,
  cutoff = 0.05,
  contextColours = NULL,
  contextLabels = waiver()
)
```

## Arguments

- results:

  A spicy results object

- fdr:

  TRUE if FDR correction is used.

- type:

  Whether to make a bubble plot (`"bubble"`) or a heatmap (any other
  value). Note: For survival results a bubble plot will be used.

- breaks:

  Vector of 3 numbers giving breaks used in legend. The first number is
  the minimum, the second is the maximum, the third is the number of
  breaks.

- comparisonGroup:

  A string specifying the name of the outcome group to compare with the
  base group.

- colours:

  Vector of colours to use to colour legend.

- marksToPlot:

  Vector of marks to include in plot.

- cutoff:

  significance threshold for circles in bubble plot.

- contextColours:

  Used for `Kontextual` results. A named list specifying the colours for
  each context. By default the Tableau colour palette is used.

- contextLabels:

  Used for `Kontextual` results. A named list to change the default
  labels for each context.

## Value

a ggplot object

## Details

The heatmap has a tile for each pair, `from` in the rows and `to` in the
columns. A tile's colour is the pair's signed log10 p-value for the
condition (of the BH-adjusted p-value with `fdr = TRUE`): -log10(p) when
the coefficient is positive and log10(p) when it is negative, so the
most significant positive differences are at the top of the scale and
the most significant negative ones at the bottom. `breaks` sets the
scale; values beyond it take the end colours, and pairs without a
p-value are grey. The values behind either plot are available from
[`topPairs()`](https://sydneybiox.github.io/spicyR/dev/reference/topPairs.md)
(`n = NULL` for every pair).

## Examples

``` r
data(spicyTest)

p <- signifPlot(spicyTest, breaks = c(-3, 3, 0.5))
p

signifPlot(spicyTest, type = "heatmap")

```
