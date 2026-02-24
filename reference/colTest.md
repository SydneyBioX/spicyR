# Perform a simple wilcoxon-rank-sum test or t-test on the columns of a data frame

Perform a simple wilcoxon-rank-sum test or t-test on the columns of a
data frame

## Usage

``` r
colTest(df, condition, type = NULL, feature = NULL, imageID = "imageID")
```

## Arguments

- df:

  A data.frame or SingleCellExperiment, SpatialExperiment

- condition:

  The condition of interest

- type:

  The type of test, "wilcox", "ttest" or "survival".

- feature:

  Can be used to calculate the proportions of this feature for each
  image

- imageID:

  The imageID's if presenting a SingleCellExperiment

## Value

Proportions

## Examples

``` r
# Test for an association with long-duration diabetes
# This is clearly ignoring the repeated measures...
data("diabetesData")
diabetesData <- spicyR:::.format_data(
  diabetesData, "imageID", "cellType", c("x", "y"), FALSE
)
props <- getProp(diabetesData)
condition <- spicyR:::getImagePheno(diabetesData)$stage
names(condition) <- spicyR:::getImagePheno(diabetesData)$imageID
condition <- condition[condition %in% c("Long-duration", "Onset")]
test <- colTest(props[names(condition), ], condition)
```
