# spicyR

<img src="https://raw.githubusercontent.com/SydneyBioX/spicyR/devel/inst/spicyR.png" align="right" width="200" alt="spicyR hex sticker" />

**Calibrated tests for changes in cell-type co-localisation between groups of patients.**

spicyR asks, for every pair of cell types, whether one sits closer to the other in one group of patients than in
another, or whether their co-localisation is associated with outcome. It works with any cell-resolution spatial
omics data: imaging mass cytometry, CODEX, MIBI, Xenium, CosMx and others.

> **spicyR 2.0 (in development on the `spicyR2` branch).** The default method is now **spicyR Cell**. The original
> image-level test is still available with `method = "image"`.

## What spicyR Cell does

For a pair *from* → *to*, spicyR Cell counts the *to* cells within radius *r* of each *from* cell.

- **The null is random labelling of the observed cells.** Each image's count is compared with its exact
  expectation if cell-type labels were shuffled among the cells that are there. Tissue shape, holes and density
  are conditioned on, not modelled. The effect is the **excess**: extra *to* cells per *from* cell.
- **Patients are the units.** Images are combined within patients by a frailty GEE, and the difference between
  conditions is tested with a small-sample cluster-robust (CR2) variance on Satterthwaite degrees of freedom.
- **Abundance is kept apart from attraction.** Every result also gives the difference at equal availability of
  the counted type.
- **And:** several conditions, covariates, k nearest neighbours, several radii (max-T or Cauchy), survival
  outcomes, and a fast C++ core: every pair of 22 cell types from 400,000 cells in seconds.

## Quick start

```r
library(spicyR)
res <- spicy(cells, condition = "condition", subject = "patient", r = 25)
res                    # summary
topPairs(res)          # the most significant pairs
res$cellResults        # the full table, including the availability-adjusted difference
signifPlot(res)        # every pair at a glance
spicyBoxPlot(res, from = "T cells", to = "Tumour")
```

`cells` can be a `data.frame`, a `SingleCellExperiment` or a `SpatialExperiment`.

## Installation

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("spicyR")                                   # release
remotes::install_github("SydneyBioX/spicyR", ref = "spicyR2")    # spicyR 2.0 (development)
```

## Python

[`spicyr`](https://github.com/SydneyBioX/spicyr-py) is the Python twin. It uses the same C++ core and gives the same
results, and it works with SpatialData, AnnData and pandas objects.

## Issues and questions

- Bugs and feature requests: [GitHub issues](https://github.com/SydneyBioX/spicyR/issues).
- Questions: [ellis.patrick@sydney.edu.au](mailto:ellis.patrick@sydney.edu.au).

## Authors

Nicolas Canete, Ellis Patrick (maintainer), Sadiq Dohadwalla, Elijah Willie, Nicholas Robertson, Alex Qin,
Farhan Ameen, Shreya Rao.

## Citation

Canete NP, Iyengar SS, Ormerod JT, Baharlou H, Harman AN, Patrick E (2022). spicyR: spatial analysis of in situ
cytometry data in R. *Bioinformatics* 38(11), 3099–3105.
[doi:10.1093/bioinformatics/btac268](https://doi.org/10.1093/bioinformatics/btac268).
The spicyR Cell method: manuscript in preparation.
