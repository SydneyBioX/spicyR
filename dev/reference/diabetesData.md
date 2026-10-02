# Imaging mass cytometry of human pancreas in type 1 diabetes (Damond et al. 2019)

Cell types and positions from imaging mass cytometry of pancreas
sections from 12 donors at three stages of type 1 diabetes (Damond et
al. 2019): 4 non-diabetic donors, 4 at onset and 4 with long-duration
disease. The object holds 253,777 cells from 120 images, 10 images per
donor. It has no assays: only the cell data used by spicyR.

## Usage

``` r
data("diabetesData")
```

## Format

A `SingleCellExperiment` with one column per cell and these `colData`
columns:

- imageID:

  the image (character)

- cellID, imageCellID:

  cell identifiers, in the whole data set and within its image

- x, y:

  the cell's coordinates in the image, in micrometres

- cellType:

  the cell type assigned by the authors (factor)

- case:

  the donor (integer)

- slide, part:

  the slide, and the part of the pancreas (head, body or tail)

- group, stage:

  the stage of type 1 diabetes, as a code and as a factor with levels
  `"Non-diabetic"`, `"Onset"` and `"Long-duration"`

## Source

Damond N et al. (2019). Mendeley Data,
[doi:10.17632/cydmwsfztj.1](https://doi.org/10.17632/cydmwsfztj.1) ,
under the CC BY 4.0 licence. How the subset was made is described in
`inst/scripts/make-diabetesData.R`.

## References

Damond N, Engler S, Zanotelli VRT, et al. (2019). A map of human type 1
diabetes progression by imaging mass cytometry. *Cell Metabolism* 29(3),
755-768.
[doi:10.1016/j.cmet.2018.11.014](https://doi.org/10.1016/j.cmet.2018.11.014)

## Examples

``` r
data("diabetesData")
table(unique(as.data.frame(SummarizedExperiment::colData(diabetesData))[, c("case", "stage")])$stage)
#> 
#>  Non-diabetic         Onset Long-duration 
#>             4             4             4 
```
