# Get proportions from a SummarizedExperiment.

Get proportions from a SummarizedExperiment.

## Usage

``` r
getProp(cells, feature = "cellType", imageID = "imageID")
```

## Arguments

- cells:

  A SingleCellExperiment, SpatialExperiment or data.frame.

- feature:

  The feature of interest

- imageID:

  The imageID's

## Value

Proportions

## Examples

``` r
data("diabetesData")
prop <- getProp(diabetesData)
```
