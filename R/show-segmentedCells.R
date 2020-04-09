#' Show SegmentedCells
#'
#' This outputs critical information about aSegmentedCells.
#'
#' @section usage:
#' `show(object)`
#'
#' @param object A SegmentedCells.
#'
#' @return Information of the number of images, cells, intenisties, morphologies and phenotypes.
#'
#' @examples
#' ### Something that resembles cellProfiler data
#'
#' set.seed(51773)
#'
#' n = 10
#'
#' cells <- data.frame(row.names = seq_len(n))
#' cells$ObjectNumber <- seq_len(n)
#' cells$ImageNumber <- rep(1:2,c(n/2,n/2))
#' cells$AreaShape_Center_X <- runif(n)
#' cells$AreaShape_Center_Y <- runif(n)
#' cells$AreaShape_round <- rexp(n)
#' cells$AreaShape_diameter <- rexp(n, 2)
#' cells$Intensity_Mean_CD8 <- rexp(n, 10)
#' cells$Intensity_Mean_CD4 <- rexp(n, 10)
#'
#' cellExp <- SegmentedCells(cells, cellProfiler = TRUE)
#'
#' ### Cluster cell types
#' intensities <- intensity(cellExp)
#' kM <- kmeans(intensities,2)
#' cellType(cellExp) <- paste('cluster',kM$cluster, sep = '')
#'
#' cellExp
#'
#' @name show-SegmentedCells
#' @rdname show-SegmentedCells
#' @aliases show
NULL


.SegmentedCells_show <- function(object) {
  cat('A SegmentedCells object with... \n')
  cat("Number of images:", nrow(object), '\n')
  cat("Number of cells:", length(cellID(object)), '\n')
  
  uniqueCellTypes <- unique(cellType(object))
  .showCat("Number of cell types: ", as.character(uniqueCellTypes))
  
  .showCat("Number of intensities: ", colnames(object[1, 'cellMarks'][[1]]))
  
  .showCat("Number of morphologies: ", colnames(object[1, 'cellMorph'][[1]]))
  
  .showCat("Number of image phenotypes: ", colnames(object[1, 'imagePheno'][[1]]))
}

if (!isGeneric("show"))
  setGeneric("show", function(object)
    standardGeneric("show"))

setMethod("show", signature(object = "SegmentedCells"), function(object) {
  .SegmentedCells_show(object)
})


.showCat <- function (fmt,
                      vals = character(),
                      exdent = 2,
                      ...)
{
  vals <- ifelse(nzchar(vals), vals, "''")
  lbls <- paste('[', paste(.selectSome(vals), collapse = ", ") , ']')
  txt <- paste(fmt, length(vals), lbls)
  cat(strwrap(txt, exdent = exdent), sep = "\n")
}

.selectSome <- function (obj, maxToShow = 4)
{
  len <- length(obj)
  if (maxToShow < 3)
    maxToShow <- 3
  if (len > maxToShow) {
    maxToShow <- maxToShow - 1
    bot <- ceiling(maxToShow / 2)
    top <- len - (maxToShow - bot - 1)
    nms <- obj[c(seq_len(bot), seq(top, len, 1))]
    c(as.character(nms[seq_len(bot)]), "...", as.character(nms[-seq_len(bot)]))
  }
  else if (is.factor(obj))
    as.character(obj)
  else
    obj
}
