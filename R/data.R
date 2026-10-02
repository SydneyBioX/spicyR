#' Imaging mass cytometry of human pancreas in type 1 diabetes (Damond et al. 2019)
#'
#' Cell types and positions from imaging mass cytometry of pancreas sections from 12 donors at three stages of type
#' 1 diabetes (Damond et al. 2019): 4 non-diabetic donors, 4 at onset and 4 with long-duration disease. The object
#' holds 253,777 cells from 120 images, 10 images per donor. It has no assays: only the cell data used by spicyR.
#'
#' @format A \code{SingleCellExperiment} with one column per cell and these \code{colData} columns:
#' \describe{
#'   \item{imageID}{the image (character)}
#'   \item{cellID, imageCellID}{cell identifiers, in the whole data set and within its image}
#'   \item{x, y}{the cell's coordinates in the image, in micrometres}
#'   \item{cellType}{the cell type assigned by the authors (factor)}
#'   \item{case}{the donor (integer)}
#'   \item{slide, part}{the slide, and the part of the pancreas (head, body or tail)}
#'   \item{group, stage}{the stage of type 1 diabetes, as a code and as a factor with levels
#'     \code{"Non-diabetic"}, \code{"Onset"} and \code{"Long-duration"}}
#' }
#' @source Damond N et al. (2019). Mendeley Data, \doi{10.17632/cydmwsfztj.1}, under the CC BY 4.0 licence. How the
#'   subset was made is described in \code{inst/scripts/make-diabetesData.R}.
#' @references Damond N, Engler S, Zanotelli VRT, et al. (2019). A map of human type 1 diabetes progression by
#'   imaging mass cytometry. \emph{Cell Metabolism} 29(3), 755-768. \doi{10.1016/j.cmet.2018.11.014}
#' @usage data("diabetesData")
#' @examples
#' data("diabetesData")
#' table(unique(as.data.frame(SummarizedExperiment::colData(diabetesData))[, c("case", "stage")])$stage)
#' @aliases diabetesData
"diabetesData"


#' Results of the image-level test on diabetesData
#'
#' The result of \code{spicy(diabetesData, condition = "stage", subject = "case", method = "image")}: the original
#' image-level test of spicyR, for all pairs of cell types, comparing the onset and long-duration stages with
#' non-diabetic donors. It is used in examples, so that they run quickly.
#'
#' @format A \code{SpicyResults} object (see \code{\link{spicy}}).
#' @source Made by \code{inst/scripts/make-spicyTest.R} from \code{\link{diabetesData}}.
#' @usage data("spicyTest")
#' @examples
#' data("spicyTest")
#' topPairs(spicyTest, n = 5)
#' @aliases spicyTest
"spicyTest"
