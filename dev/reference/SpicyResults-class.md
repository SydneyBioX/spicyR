# The results of spicy()

[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
returns a `SpicyResults` object: a list with the elements below. Use
[`topPairs()`](https://sydneybiox.github.io/spicyR/dev/reference/topPairs.md),
[`signifPlot()`](https://sydneybiox.github.io/spicyR/dev/reference/signifPlot.md),
[`spicyBoxPlot()`](https://sydneybiox.github.io/spicyR/dev/reference/spicyBoxPlot.md)
and
[`bind()`](https://sydneybiox.github.io/spicyR/dev/reference/bind.md) on
it, or take elements with `$`.

## Details

Elements for the cell method (`method = "cell"`, the default):

- `cellResults`:

  the full table, one row per pair (and per level with more than two
  conditions): `excess_ref`, `excess_comp`, `excess_difference`, `se`,
  `df`, `p_value`, `p_adj`, `tau2`; with covariates or
  `adjustAbundance = TRUE`, what the test was adjusted for
  (`adjusted_for`), the effect and p-value of each adjustment
  (`<covariate>_effect`, `<covariate>_p_value`, `abundance_effect`) and
  the unadjusted test (`unadjusted_difference`, `unadjusted_se`,
  `unadjusted_df`, `unadjusted_p_value`, `unadjusted_p_adj`). For
  survival: `score_coefficient`, `p_value`, `p_adj`, `hazard_ratio_sd`,
  `hr_p_value` and, with `adjustAbundance = TRUE`, the p-value and
  hazard ratio without it. With several radii: `r`, the radius with the
  strongest evidence, and `p_value_best_radius`.

- `radiusResults`:

  with several radii, the table at every radius.

- `pairwiseAssoc`:

  a list with, for every pair, the effect in every image (`NA` where the
  pair could not be measured).

- `effect`:

  `"allocation"` (the extra fraction of `to` cells with a `from` cell
  within `r`), `"count"` (the extra `from` cells within `r` of each `to`
  cell), or `"kontextual"` (Statial's Kontextual test).

- `imageWeights`:

  a list with, for every pair, each image's weight in the test: its
  share of its condition's information (sums to 1 within each
  condition).

- `coefficient`, `p.value`, `se`, `statistic`, `df`:

  matrices with one row per pair, as for the image method: the reference
  excess and the differences, and their tests.

- `condition`, `imageID`, `subject`, `nCells`, `r`, `k`:

  the condition and patient of every image, the cell counts, and the
  settings of the analysis.

For the image method (`method = "image"`) the object holds the matrices
`coefficient`, `p.value`, `se`, `statistic` and `df`, and
`pairwiseAssoc`, as in spicyR 1.x.

## See also

[`spicy`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md),
[`topPairs`](https://sydneybiox.github.io/spicyR/dev/reference/topPairs.md),
[`bind`](https://sydneybiox.github.io/spicyR/dev/reference/bind.md)
