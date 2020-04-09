################################################################################
#
# Generics for spicy
#
################################################################################


#' A table of the significant results from spicy tests
#'
#' @param x The output from spicy.
#' @param coef Which coefficient to list.
#' @param n Extract the top n most significant pairs.
#' @param adj Which p-value adjustment method to use, argument for p.adjust().
#' @param cutoff A p-value threshold to extract significant pairs.
#'
#' @return A data.frame
#'
#' @examples
#' 
#' data(spicyTest)
#' topPairs(spicyTest)
#' 
#' @aliases 
#' topPairs,spicy-method
#' topPairs
#' @rdname topPairs
#' @export
setGeneric("topPairs", function(x,
                           coef = 1,
                           n = 10,
                           adj = 'fdr',
                           cutoff = NULL)
    standardGeneric("topPairs"))
setMethod("topPairs", "SpicyResults", function(x,
                                   coef = 1,
                                   n = 10,
                                   adj = 'fdr',
                                   cutoff = NULL) {
    useCondition <- grep('condition', colnames(x$p.value))[coef]
    pval <- as.numeric(x$p.value[[useCondition]])
    adj.pvalue <- p.adjust(pval, adj)
    
    comp <- x$comparison
    
    results <-
        data.frame(
            intercept = x$coefficient[, "(Intercept)"],
            coefficient = x$coefficient[, useCondition],
            p.value = pval,
            adj.pvalue = adj.pvalue,
            from = comp$from,
            to = comp$to
        )
    rownames(results) <- rownames(x$coefficient)
    if (length(results$p.value) > 0) {
        results <- results[order(results$p.value), ]
        
        if (is.null(cutoff) &
            !is.null(n))
            return(results[seq_len(pmin(n, nrow(results))), ])
        if (is.null(n) &
            !is.null(cutoff))
            return(results[results$adj.pvalue < 0.05, ])
        if (is.null(n) &
            !is.null(cutoff))
            return(results[intersect(Which(results$adj.pvalue < 0.05), seq_len(pmin(n, nrow(results)))), ])
    }
    
})





################################################################################
#
# Generics for SegmentedCells
#
################################################################################



#' Accessors for SegmentedCells
#'
#' Methods to access various components of the `SegmentedCells` object.
#'
#' @usage cellSummary(x, imageID = NULL, bind = TRUE)
#' @usage cellSummary(x, imageID = NULL) <- value
#' @usage cellMarks(x, imageID = NULL, bind = TRUE)
#' @usage cellMarks(x, imageID = NULL) <- value
#' @usage cellMorph(x, imageID = NULL, bind = TRUE)
#' @usage cellMorph(x, imageID = NULL) <- value
#' @usage imagePheno(x, imageID = NULL, bind = TRUE, expand = FALSE)
#' @usage imagePheno(x, imageID = NULL) <- value
#' @usage imageID(x, imageID = NULL)
#' @usage cellID(x, imageID = NULL)
#' @usage cellID(x) <- value
#' @usage imageCellID(x, imageID = NULL)
#' @usage imageCellID(x) <- value
#' @usage cellType(x, imageID = NULL)
#' @usage cellType(x, imageID = NULL) <- value
#' @usage filterCells(x, select)
#'
#' @param x A `SegmentedCells` object.
#' @param imageID A vector of imageIDs to specifically extract.
#' @param bind When false outputs a list of DataFrames split by imageID
#' @param expand Used to expand the phenotype information from per image to per cell.
#' @param value The relevant information used to replace.
#' @param select A logical vector of the cells to be kept.
#'
#' @section Descriptions:
#' \describe{
#' \item{`cellSummary`:}{
#' Retrieves the DataFrame containing `x` and `y` coordinates of each cell as well as `cellID`, `imageID` and `cellType`.
#' imageID can be used to select specific images and bind=FALSE outputs the information as a list split by imageID.
#' }
#'
#' \item{`cellMorph`:}{
#' Retrieves the DataFrame containing morphology information.
#' }
#'
#' \item{`cellMarks`:}{
#' Retrieves the DataFrame containing intensity of gene or protein markers.
#' }
#'
#' \item{`imagePheno`:}{
#' Retrieves the DataFrame containing the phenotype information for each image. Using expand = TRUE will produce a DataFrame with the number of rows equal to the number of cells.
#' }
#' }
#'
#' @return DataFrame or a list of DataFrames
#' @name Accessors
#'
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
#' intensities <- cellMarks(cellExp)
#' kM <- kmeans(intensities,2)
#' cellType(cellExp) <- paste('cluster',kM$cluster, sep = '')
#'
#' cellSummary(cellExp, imageID = 1)
#'
#' @aliases
#' cellSummary,SegmentedCells-method
#' cellSummary<-,SegmentedCells-method
#' cellMarks,SegmentedCells-method
#' cellMarks<-,SegmentedCells-method
#' cellMorph,SegmentedCells-method
#' cellMorph<-,SegmentedCells-method
#' imagePheno,SegmentedCells-method
#' imagePheno<-,SegmentedCells-method
#' imageID,SegmentedCells-method
#' cellType,SegmentedCells-method
#' cellType<-,SegmentedCells-method
#' imageCellID,SegmentedCells-method
#' imageCellID<-,SegmentedCells-method
#' cellID,SegmentedCells-method
#' cellID<-,SegmentedCells-method
#' filterCells,SegmentedCells-method
#' cellSummary
#' cellSummary<-
#' cellMarks
#' cellMarks<-
#' cellMorph
#' cellMorph<-
#' imagePheno
#' imagePheno<-
#' imageID
#' cellType
#' cellType<-
#' imageCellID
#' imageCellID<-
#' cellID
#' cellID<-
#' filterCells


### Get cellSummary information for each cell.

#' @export
#' @importFrom BiocGenerics do.call rbind
setGeneric("cellSummary", function(x, imageID = NULL, bind = TRUE)
    standardGeneric("cellSummary"))
setMethod("cellSummary", "SegmentedCells", function(x, imageID = NULL, bind = TRUE) {
    if (!is.null(imageID)) {
        x <- x[imageID,]
    }
    if (bind == FALSE) {
        return(x$cellSummary)
    }
    if (bind == TRUE) {
        return(cbind(
            imageID = rep(rownames(x), unlist(lapply(x[, "cellSummary"], nrow))),
            BiocGenerics::do.call("rbind", x$cellSummary)
        ))
    }
    
})

#' @export
#' @importFrom S4Vectors split
setGeneric("cellSummary<-", function(x, imageID = NULL, value)
    standardGeneric("cellSummary<-"))
setMethod("cellSummary<-", "SegmentedCells", function(x, imageID = NULL, value) {
    if (is.null(imageID))
        imageID <- rownames(x)
    if (nrow(value) == length(imageID)) {
        x[imageID,]@listData$cellSummary <- value
        return(x)
    }
    
    if (nrow(value) == length(imageID(x, imageID))) {
        value <- value[, c("cellID", "imageCellID", "x", "y", "cellType")]
        by <-
            rep(imageID, unlist(lapply(x[imageID, "cellSummary"], nrow)))
        by <- factor(by, levels = unique(by))
        x[imageID,]@listData$cellSummary <- S4Vectors::split(value, by)
        return(x)
    }
})




### Get imageIDs for each cell, not sure if this should also report rownames(df)

#' @export
setGeneric("imageID", function(x, imageID = NULL)
    standardGeneric("imageID"))
setMethod("imageID", "SegmentedCells", function(x, imageID = NULL) {
    if (!is.null(imageID)) {
        x <- x[imageID,]
    }
    rep(rownames(x), unlist(lapply(x$cellSummary, nrow)))
})




### Get cellIDs

#' @export
#' @importFrom BiocGenerics do.call rbind
setGeneric("cellID", function(x, imageID = NULL)
    standardGeneric("cellID"))
setMethod("cellID", "SegmentedCells", function(x, imageID = NULL) {
    if (!is.null(imageID)) {
        x <- x[imageID,]
    }
    BiocGenerics::do.call("rbind", x$cellSummary)$cellID
})

#' @export
setGeneric("cellID<-", function(x, value)
    standardGeneric("cellID<-"))
setMethod("cellID<-", "SegmentedCells", function(x, value) {
    loc <- cellSummary(x)
    
    if (nrow(loc) != length(value)) {
        stop("There is not enough or too many cellIDs")
    }
    
    loc$cellID <- value
    cellSummary(x) <- loc
})




### Get imageCellID

#' @export
setGeneric("imageCellID", function(x, imageID = NULL)
    standardGeneric("imageCellID"))
setMethod("imageCellID", "SegmentedCells", function(x, imageID = NULL) {
    if (!is.null(imageID)) {
        x <- x[imageID,]
    }
    BiocGenerics::do.call("rbind", x$cellSummary)$imageCellID
})

#' @export
setGeneric("imageCellID<-", function(x, value)
    standardGeneric("imageCellID<-"))
setMethod("imageCellID<-", "SegmentedCells", function(x, value) {
    loc <- cellSummary(x)
    
    if (nrow(loc) != length(value)) {
        stop("There is not enough or too many imageCellIDs")
    }
    
    loc$imageCellID <- value
    
    cellSummary(x) <- loc
})




### Get marker intensity information

#' @export
setGeneric("cellMarks", function(x, imageID = NULL, bind = TRUE)
    standardGeneric("cellMarks"))
setMethod("cellMarks", "SegmentedCells", function(x, imageID = NULL, bind = TRUE) {
    if (!is.null(imageID)) {
        x <- x[imageID,]
    }
    if (bind == FALSE) {
        return(x$cellMarks)
    }
    if (bind == TRUE) {
        return(BiocGenerics::do.call("rbind", x$cellMarks))
    }
    
})

#' @export
setGeneric("cellMarks<-", function(x, imageID = NULL, value)
    standardGeneric("cellMarks<-"))
setMethod("cellMarks<-", "SegmentedCells", function(x, imageID = NULL, value) {
    if (is.null(imageID))
        imageID <- rownames(x)
    if (nrow(value) == length(imageID)) {
        x[imageID,]@listData$cellMarks <- value
        return(x)
    }
    
    if (nrow(value) == length(imageID(x))) {
        by <- rep(rownames(x), unlist(lapply(x$cellMarks, nrow)))
        by <- factor(by, levels = unique(by))
        x[imageID,]@listData$cellMarks <-
            S4Vectors::split(value, by)
        return(x)
    }
})




### Get morphology information

#' @export
setGeneric("cellMorph", function(x, imageID = NULL, bind = TRUE)
    standardGeneric("cellMorph"))
setMethod("cellMorph", "SegmentedCells", function(x, imageID = NULL, bind = TRUE) {
    if (!is.null(imageID)) {
        x <- x[imageID,]
    }
    if (bind == FALSE) {
        return(x$cellMorph)
    }
    if (bind == TRUE) {
        return(BiocGenerics::do.call("rbind", x$cellMorph))
    }
    
})

#' @export
setGeneric("cellMorph<-", function(x, imageID = NULL, value)
    standardGeneric("cellMorph<-"))
setMethod("cellMorph<-", "SegmentedCells", function(x, imageID = NULL,
                                                     value) {
    if (is.null(imageID))
        imageID <- rownames(x)
    if (nrow(value) == length(imageID)) {
        x[imageID,]@listData$cellMorph <- value
        return(x)
    }
    
    if (nrow(value) == length(imageID(x, imageID))) {
        by <- rep(rownames(x), unlist(lapply(x$cellMorph, nrow)))
        by <- factor(by, levels = unique(by))
        
        x[imageID,]@listData$cellMorph <-
            S4Vectors::split(value, by)
        return(x)
    }
})




### Get cell type information

#' @export
setGeneric("cellType", function(x, imageID = NULL)
    standardGeneric("cellType"))
setMethod("cellType", "SegmentedCells", function(x, imageID = NULL) {
    if (!is.null(imageID)) {
        x <- x[imageID,]
    }
    BiocGenerics::do.call("rbind", x$cellSummary)$cellType
})

#' @export
setGeneric("cellType<-", function(x, imageID = NULL, value)
    standardGeneric("cellType<-"))
setMethod("cellType<-", "SegmentedCells", function(x, imageID = NULL, value) {
    if (is.null(imageID))
        imageID <- rownames(x)
    loc <- cellSummary(x, imageID = imageID)
    
    if (nrow(loc) != length(value)) {
        stop("There is not enough or too many cellTypes")
    }
    
    loc$cellType <- value
    
    cellSummary(x, imageID = imageID) <- loc
    x
})





### Get and add image phenotype data to the object

#' @export
setGeneric("imagePheno", function(x,
                                 imageID = NULL,
                                 bind = TRUE,
                                 expand = FALSE)
    standardGeneric("imagePheno"))
setMethod("imagePheno", "SegmentedCells", function(x,
                                                  imageID = NULL,
                                                  bind = TRUE,
                                                  expand = FALSE) {
    if (!is.null(imageID)) {
        x <- x[imageID,]
    }
    if (expand) {
        pheno <- BiocGenerics::do.call("rbind", x$imagePheno)
        rownames(pheno) <- pheno$imageID
        return(pheno[imageID(x),])
    } else {
        pheno <- BiocGenerics::do.call("rbind", x$imagePheno)
        rownames(pheno) <- pheno$imageID
        return(pheno[rownames(x),])
    }
})


#' @export
setGeneric("imagePheno<-", function(x, imageID = NULL, value)
    standardGeneric("imagePheno<-"))
setMethod("imagePheno<-", "SegmentedCells", function(x, imageID = NULL, value) {
    if (is.null(imageID))
        imageID <- rownames(x)
    use <- intersect(value$imageID, imageID)
    rownames(value) <- value$imageID
    x[use,]@listData$imagePheno <-
        S4Vectors::split(value[use,], use)
    x[unique(use),]
})




as.data.frame.SegmentedCells <- function(x, ...) {
    loc <- cellSummary(x)
    int <- cellMarks(x)
    morph <- cellMorph(x)
    pheno <- imagePheno(x, expand = TRUE)
    pheno <- pheno[!colnames(pheno) %in% "imageID"]
    
    if (length(colnames(int)) > 0)
        colnames(int) <- paste("intensity", colnames(int), sep = "_")
    if (length(colnames(morph)) > 0)
        colnames(morph) <- paste("morphology", colnames(morph), sep = "_")
    if (length(colnames(pheno)) > 0)
        colnames(pheno) <- paste("phenotype", colnames(pheno), sep = "_")
    
    if (prod(dim(loc)) > 0 &
        prod(dim(int)) > 0 & prod(dim(morph)) > 0 & prod(dim(pheno)) > 0) {
        return(as.data.frame(cbind(loc, int, morph, pheno)))
    }
    
    if (prod(dim(loc)) > 0 & prod(dim(int)) > 0 &
        prod(dim(morph)) > 0) {
        return(as.data.frame(cbind(loc, int, morph)))
    }
    
    if (prod(dim(loc)) > 0 & prod(dim(int)) > 0 &
        prod(dim(pheno)) > 0) {
        return(as.data.frame(cbind(loc, int, pheno)))
    }
    
    if (prod(dim(loc)) > 0 &
        prod(dim(morph)) > 0 & prod(dim(pheno)) > 0) {
        return(as.data.frame(cbind(loc, morph, pheno)))
    }
    
    if (prod(dim(loc)) > 0 & prod(dim(int)) > 0) {
        return(as.data.frame(cbind(loc, int)))
    }
    
    if (prod(dim(loc)) > 0 & prod(dim(morph)) > 0) {
        return(as.data.frame(cbind(loc, morph)))
    }
    
    if (prod(dim(loc)) > 0 & prod(dim(pheno)) > 0) {
        return(as.data.frame(cbind(loc, pheno)))
    }
    
    NA
}


if (!isGeneric("as.data.frame"))
    setGeneric("as.data.frame", function(x, ...)
        standardGeneric("as.data.frame"))
setMethod("as.data.frame",
          "SegmentedCells",
          as.data.frame.SegmentedCells)


#' @export
setGeneric("filterCells", function(x, select)
    standardGeneric("filterCells"))
setMethod("filterCells", "SegmentedCells", function(x, select) {
    df <- as.data.frame(x)
    if (length(select) != nrow(df))
        stop("length of select must equal nrow of SegmentedCells")
    SegmentedCells(df[select, ])
})


################################################################################
#
# Generics for lisa
#
################################################################################


#' Adding regions to SegmentedCells
#'
#' A method for setting and getting regions for a SegmentedCells object.
#'
#' @usage cellSummary(x, imageID = NULL, bind = TRUE)
#' @usage cellSummary(x, imageID = NULL) <- value
#' @usage cellMarks(x, imageID = NULL, bind = TRUE)
#' @usage cellMarks(x, imageID = NULL) <- value
#' @usage cellMorph(x, imageID = NULL, bind = TRUE)
#' @usage cellMorph(x, imageID = NULL) <- value
#' @usage imagePheno(x, imageID = NULL, bind = TRUE, expand = FALSE)
#' @usage imagePheno(x, imageID = NULL) <- value
#' @usage imageID(x, imageID = NULL)
#' @usage cellID(x, imageID = NULL)
#' @usage cellID(x) <- value
#' @usage imageCellID(x, imageID = NULL)
#' @usage imageCellID(x) <- value
#' @usage cellType(x, imageID = NULL)
#' @usage cellType(x, imageID = NULL) <- value
#'
#' @param x A `SegmentedCells` object.
#' @param imageID A vector of imageIDs to specifically extract.
#' @param bind When false outputs a list of DataFrames split by imageID
#' @param expand Used to expand the phenotype information from per image to per cell.
#' @param value The relevant information used to replace.
#'
#' @section Descriptions:
#' \describe{
#' \item{`cellSummary`:}{
#' Retrieves the DataFrame containing `x` and `y` coordinates of each cell as well as `cellID`, `imageID` and `cellType`.
#' imageID can be used to select specific images and bind=FALSE outputs the information as a list split by imageID.
#' }
#'
#' \item{`cellMorph`:}{
#' Retrieves the DataFrame containing morphology information.
#' }
#'
#' \item{`cellMarks`:}{
#' Retrieves the DataFrame containing intensity of gene or protein markers.
#' }
#'
#' \item{`imagePheno`:}{
#' Retrieves the DataFrame containing the image phenotype information. Using expand = TRUE will produce a DataFrame with the number of rows equal to the number of cells.
#' }
#' }
#'
#' @return DataFrame or a list of DataFrames
#' @name Accessors
#'
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
#' intensities <- cellMarks(cellExp)
#' kM <- kmeans(intensities,2)
#' cellType(cellExp) <- paste('cluster',kM$cluster, sep = '')
#'
#' cellSummary(cellExp, imageID = 1)
#'
#' @aliases
#' region,SegmentedCells-method
#' region<-,SegmentedCells-method
#' region
#' region<-





### Get regions information

#' @export
#' @importFrom BiocGenerics do.call rbind
setGeneric("region", function(x, imageID = NULL, annot = FALSE)
    standardGeneric("region"))
setMethod("region", "SegmentedCells", function(x, imageID = NULL, annot = TRUE) {
    if (!is.null(imageID)) {
        x <- x[imageID,]
    }
    if (is.null(x$region))
        stop("There is no region information in your SegmentedCells yet")
    if (annot)
        return(data.frame(cellSummary(x), region = BiocGenerics::do.call("rbind", x$region)))
    
    BiocGenerics::do.call("rbind", x$region)
})




#' @export
#' @importFrom S4Vectors DataFrame split
setGeneric("region<-", function(x, imageID = NULL, value)
    standardGeneric("region<-"))
setMethod("region<-", "SegmentedCells", function(x, imageID = NULL, value) {
    if (is.null(imageID))
        imageID <- rownames(x)
    if (length(value) == length(imageID(x, imageID))) {
        if (is.null(x$region)) {
            x <- DataFrame(x)
            by <- rep(rownames(x), unlist(lapply(x$cellSummary, nrow)))
            by <- factor(by, levels = unique(by))
            x$region <-
                S4Vectors::split(DataFrame(region = value), by)
            x <- new("SegmentedCells", x)
        }
        if (!is.null(x$region))
            by <- rep(rownames(x), unlist(lapply(x$cellSummary, nrow)))
        by <- factor(by, levels = unique(by))
        x[imageID,]@listData$region <-
            S4Vectors::split(DataFrame(region = value), by)
    }
    x
})
