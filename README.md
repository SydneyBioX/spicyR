# spicyR <img src="inst/spicyR.png" align="right" width="140" alt="spicyR hex sticker" />

**Test whether cell types co-localise differently between groups of patients.**

Do T cells gather around tumour cells more in one group of patients than in another? spicyR answers questions
like this for every pair of cell types in imaging and spatial transcriptomics data, comparing groups of patients
or relating co-localisation to survival. For each image it counts how many cells of one type lie within a radius of
each cell of another, and compares that count with what random labelling of the same cells would give, so holes,
folds and dense regions in the tissue are not mistaken for biology. Patients, not images or cells, are the units of
the test.

![T cells around proliferating tumour cells in breast cancer, compared between ER-negative and ER-positive patients](man/figures/spicyR_overview.png)

## Quick start

```r
library(spicyR)

res <- spicy(cells, condition = "response", subject = "patient", r = 25)
topPairs(res)                 # the most significant pairs
signifPlot(res)               # every pair at a glance
spicyBoxPlot(res, from = "CD8 T cells", to = "Tumour")
```

`cells` can be a `SpatialExperiment`, a `SingleCellExperiment` or a `data.frame` with one row per cell, holding
the image, cell type and coordinates of each cell.

## What you get

- A table with one row per pair of cell types: the number of extra neighbours per cell in each group, the
  difference, a p-value and an FDR-adjusted p-value.
- A plot of every pair at once, and the per-patient values behind any pair.
- The same test with covariates, several radii, more than two groups, or a survival outcome.

## Installation

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("spicyR")
```

The development version: `remotes::install_github("SydneyBioX/spicyR", ref = "spicyR2")`.

## Learn more

- [Introduction to spicyR](https://sydneybiox.github.io/spicyR/dev/articles/spicyR.html): a full analysis of
  breast cancer imaging mass cytometry data.
- [The original image-level test](https://sydneybiox.github.io/spicyR/dev/articles/image_method.html), the
  default before spicyR 2.0, still available with `method = "image"`.
- [spicyr](https://github.com/SydneyBioX/spicyr-py), the same analysis in Python, for SpatialData and AnnData
  objects.

## Citation

Canete NP, Iyengar SS, Ormerod JT, Baharlou H, Harman AN, Patrick E (2022). spicyR: spatial analysis of in situ
cytometry data in R. *Bioinformatics* 38(11), 3099–3105.
[doi:10.1093/bioinformatics/btac268](https://doi.org/10.1093/bioinformatics/btac268)

## Contact

Questions and bug reports: [GitHub issues](https://github.com/SydneyBioX/spicyR/issues) or
[ellis.patrick@sydney.edu.au](mailto:ellis.patrick@sydney.edu.au).

spicyR 2.0 is in development on the `spicyR2` branch; a paper describing its test is in preparation. Authors:
Nicolas Canete, Ellis Patrick, Sadiq Dohadwalla, Elijah Willie, Nicholas Robertson, Alex Qin, Farhan Ameen and
Shreya Rao.
