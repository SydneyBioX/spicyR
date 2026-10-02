## Provenance of data/diabetesData.rda.
##
## Source: Damond N et al. (2019), "A map of human type 1 diabetes progression by imaging mass cytometry",
## Cell Metabolism 29(3), 755-768, doi:10.1016/j.cmet.2018.11.014. Data: Mendeley Data, doi:10.17632/cydmwsfztj.1
## (version 1, 31 January 2019), licence CC BY 4.0.
##
## The subset was made for spicyR 1.x in 2020 from the authors' single-cell table (cell positions, the cell types
## they assigned, donor and image metadata): all 12 donors, 10 images per donor (120 images, 253,777 cells), with no
## marker intensities. It was later converted from spicyR's former SegmentedCells class to a SingleCellExperiment
## (2024), and unused factor levels were dropped. The script that chose the 10 images per donor was not kept, so
## this file records what the object contains rather than rebuilding it. The check below confirms that the object
## shipped with the package is the one described.
data("diabetesData", package = "spicyR")
cd <- as.data.frame(SummarizedExperiment::colData(diabetesData))
stopifnot(
  ncol(diabetesData) == 253777,
  length(unique(cd$imageID)) == 120,
  length(unique(cd$case)) == 12,
  all(table(unique(cd[, c("case", "imageID")])$case) == 10),
  identical(levels(cd$stage), c("Non-diabetic", "Onset", "Long-duration"))
)
