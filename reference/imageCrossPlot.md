# Pairwise cross-plot of spatial associations

Builds a pairwise scatterplot of marker-marker associations for a single
image using the output of \`getPairwise()\`.

## Usage

``` r
imageCrossPlot(
  result,
  image = NULL,
  colourGradient = c("#4575B4", "white", "#D73027"),
  marksToPlot = NULL,
  limits = NULL
)
```

## Arguments

- result:

  A matrix or data frame produced by \`getPairwise()\` where rownames
  are image IDs and colnames are concatenated marker pairs in the form
  \`"from\_\_to"\`. Must contain at least one row with name matching
  \`image\`.

- image:

  Character scalar. Which image (row in \`result\`) to plot. If
  \`NULL\`, defaults to the first rowname of \`result\`.

- colourGradient:

  Character vector of length 3 giving the low, mid, and high colours for
  a diverging palette used by \`scale_colour_gradient2()\`. Default is
  \`c("#4575B4", "white", "#D73027")\`.

- marksToPlot:

  Optional character vector of marker names. If supplied, the plot is
  restricted to rows and columns where both \`from\` and \`to\` are in
  this set.

- limits:

  Numeric length-2 vector giving the lower and upper caps applied to
  association values for colour mapping. Use \`NULL\` to avoid clamping.
  Point sizes are clamped to \`c(0, max(abs(limits)))\`.

## Value

A ggplot2 object.

## Details

The function expects \`colnames(result)\` to contain \`"\_\_"\`
separating marker names; otherwise an error is thrown. Values are
transformed into two aesthetics:

- `value`: the (possibly clamped) signed association used for colour.

- `size`: the (possibly clamped) absolute association used for point
  size.

The colourbar shows three ticks at `lower`, `0`, and `upper` with labels
`"<lower"`, `"0"`, and `">upper"`.

## See also

[`getPairwise`](https://sydneybiox.github.io/spicyR/reference/getPairwise.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal example (toy data)
set.seed(1)
mks <- c("A","B","C")
pairs <- as.vector(outer(mks, mks, paste, sep="__"))
res <- matrix(rnorm(length(pairs)*2), nrow = 2,
              dimnames = list(c("img1","img2"), pairs))

p <- imageCrossPlot(res, image = "img1",
                       colourGradient = c("#4575B4","white","#D73027"),
                       marksToPlot = c("A","B","C"),
                       limits = c(-3, 3))
print(p)
} # }
```
