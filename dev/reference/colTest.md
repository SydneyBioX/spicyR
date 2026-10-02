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

A data.frame with one row per column of `df` (for example per cell
type), ordered by p-value: the mean in each condition, the test
statistic, the p-value (`pval`), the Benjamini-Hochberg adjusted p-value
(`adjPval`) and the column name (`cluster`).

## Examples

``` r
# Test for a difference in cell-type proportions between onset and long-duration diabetes
# (images are treated as independent here, which ignores that each donor has several)
data("diabetesData")
props <- getProp(diabetesData)
images <- unique(as.data.frame(SummarizedExperiment::colData(diabetesData))[, c("imageID", "stage")])
condition <- setNames(as.character(images$stage), images$imageID)
condition <- condition[condition %in% c("Long-duration", "Onset")]
test <- colTest(props[names(condition), ], condition)
head(test)
#>         mean in group Long-duration mean in group Onset tval.t    pval adjPval
#> acinar                      5.9e-01             0.49000    6.6 5.2e-09 8.3e-08
#> stromal                     1.4e-02             0.03300   -5.8 2.1e-07 1.7e-06
#> beta                        2.5e-05             0.02900   -6.2 3.1e-07 1.7e-06
#> ductal                      2.1e-01             0.26000   -5.3 1.1e-06 4.4e-06
#> naiveTc                     3.2e-04             0.00099   -3.3 1.7e-03 5.1e-03
#> Th                          3.1e-03             0.00620   -3.3 1.9e-03 5.1e-03
#>         cluster
#> acinar   acinar
#> stromal stromal
#> beta       beta
#> ductal   ductal
#> naiveTc naiveTc
#> Th           Th
```
