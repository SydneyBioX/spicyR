#' Show SegmentedCellExperiment object
#' 
#' This outputs critical information about aSegmentedCellExperiment object.
#' 
#' @section usage:
#' `show(object)`
#'
#' @param object A SegmentedCellExperiment object.
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
#' cellExp <- SegmentedCellExperiment(cells, cellProfiler = TRUE)
#' 
#' ### Cluster cell types
#' intensities <- intensity(cellExp)
#' kM <- kmeans(intensities,2)
#' cellType(cellExp) <- paste('cluster',kM$cluster, sep = '')
#' 
#' cellExp
#' 
#' @name show-SegmentedCellExperiment
#' @rdname show-SegmentedCellExperiment
#' @aliases show
NULL


.SegmentedCellExperiment_show <- function(object) {
  cat('A SegmentedCellExperiment with... \n')
  cat("Number of images:", nrow(object),'\n')
  cat("Number of cells:", length(cellID(object)),'\n')
  
  uniqueCellTypes <- unique(cellType(object))
  .showCat("Number of cell types: ", uniqueCellTypes)
  
  .showCat("Number of intensities: ", colnames(object[1,'intensity'][[1]]))
  
  .showCat("Number of morphologies: ", colnames(object[1,'morphology'][[1]]))
  
  .showCat("Number of image phenotypes: ", colnames(object[1,'phenotype'][[1]]))
}

if (!isGeneric("show")) setGeneric("show", function(object) standardGeneric("show"))

setMethod("show", signature(object = "SegmentedCellExperiment"), function(object) {
  .SegmentedCellExperiment_show(object)
})


.showCat <- function (fmt, vals = character(), exdent = 2, ...) 
{
  vals <- ifelse(nzchar(vals), vals, "''")
  lbls <- paste('[', paste(.selectSome(vals), collapse = ", ") ,']')
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
    bot <- ceiling(maxToShow/2)
    top <- len - (maxToShow - bot - 1)
    nms <- obj[c(1:bot, top:len)]
    c(as.character(nms[1:bot]), "...", as.character(nms[-c(1:bot)]))
  }
  else if (is.factor(obj)) 
    as.character(obj)
  else obj
}

