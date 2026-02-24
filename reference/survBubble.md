# Plots survival results from spicy.

Plots survival results from spicy.

## Usage

``` r
survBubble(
  result,
  fdr = FALSE,
  cutoff = 0.05,
  colourGradient = c("#4575B4", "white", "#D73027"),
  marksToPlot = NULL,
  contextColours = NULL,
  contextLabels = waiver()
)
```

## Arguments

- result:

  A spicyResults object that contains survival results.

- fdr:

  TRUE if FDR correction is used.

- cutoff:

  Significance threshold for circles in bubble plot.

- colourGradient:

  A vector of colours, used to define the low, medium, and high values
  for the colour scale.

- marksToPlot:

  Vector of marks to include in bubble plot.

- contextColours:

  Used for `Kontextual` results. A named list specifying the colours for
  each context. By default the Tableau colour palette is used.

- contextLabels:

  Used for `Kontextual` results. A named list to change the default
  labels for each context.

## Value

A ggplot object.
