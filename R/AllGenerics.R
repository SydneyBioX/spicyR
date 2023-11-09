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
#' topPairs,SpicyResults-method
#' topPairs
#' @rdname topPairs
#' @export
setGeneric("topPairs", function(x,
                                coef = NULL,
                                n = 10,
                                adj = 'fdr',
                                cutoff = NULL)
    standardGeneric("topPairs"))
setMethod("topPairs", "SpicyResults", function(x,
                                               coef = NULL,
                                               n = 10,
                                               adj = 'fdr',
                                               cutoff = NULL) {
    if(!methods::is(x,"SpicyResults")) stop("x are not results from spicy")
    
    if(is.null(coef)) coef <- grep('condition', colnames(x$p.value))[1]
    
    if(methods::is(coef,"character")&!coef%in%colnames(x$p.value)) stop("coef not a column name")
    if(methods::is(coef,"numeric")&!coef%in%seq_len(ncol(x$p.value))) stop("coef not a column name")
    if(length(coef)>1) warning("coef needs to be length 1, taking first entry.")
    useCondition <- coef[1]
    pval <- x$p.value[[useCondition]]
    adj.pvalue <- stats::p.adjust(pval, adj)
    
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
            return(results[which(results$adj.pvalue <= cutoff), ])
        if ((!is.null(n)) &
            (!is.null(cutoff)))
            return(results[which(results$adj.pvalue <= cutoff & seq_len(nrow(results))<=n), ])
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
#' @usage cellAnnotation(x, variable, imageID = NULL)
#' @usage cellAnnotation(x, variable, imageID = NULL) <- value

#'
#' @param x A `SegmentedCells` object.
#' @param imageID A vector of imageIDs to specifically extract.
#' @param bind When false outputs a list of DataFrames split by imageID
#' @param expand Used to expand the phenotype information from per image to per cell.
#' @param value The relevant information used to replace.
#' @param select A logical vector of the cells to be kept.
#' @param variable A variable to add or retrieve from cellSummary.
#'
#' @section Descriptions:
#' \describe{
#' \item{`cellSummary`:}{
#' Retrieves the DataFrame containing `x` and `y` coordinates of each cell as 
#' well as `cellID`, `imageID` and `cellType`. imageID can be used to select 
#' specific images and bind=FALSE outputs the information as a list split by imageID.
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
#' Retrieves the DataFrame containing the phenotype information for each image.
#'  Using expand = TRUE will produce a DataFrame with the number of rows equal 
#'  to the number of cells.
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
#' cellAnnotation,SegmentedCells-method
#' cellAnnotation<-,SegmentedCells-method
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
#' cellAnnotation
#' cellAnnotation<-
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
            imageID = factor(rep(rownames(x), unlist(lapply(x[, "cellSummary"], nrow))),
                             rownames(x)),
            BiocGenerics::do.call("rbind", x$cellSummary)
        ))
    }
    
})


#' @importFrom methods slot slot<-
.putData <- function(object,variable, value, image = NULL){
    methods::slot(object, "listData")[[variable]] <- value
    object
}


#' @export
#' @importFrom S4Vectors split
setGeneric("cellSummary<-", function(x, imageID = NULL, value)
    standardGeneric("cellSummary<-"))
setReplaceMethod("cellSummary", "SegmentedCells", function(x, imageID = NULL, value) {

    if (is.null(imageID))
        imageID <- rownames(x)
    if (nrow(value) == length(imageID)) {
        if(any(!colnames(x[1,1][[1]])%in%colnames(value[[1]]))){
            stop("There are colnames of value that aren't in cellSummary")}
        x <- .putData(x, "cellSummary", value, imageID)
        return(x)
    }
    
    if (nrow(value) == length(imageID(x, imageID))) {
        colNames <- colnames(x[1,1][[1]])
        if(any(!colNames%in%colnames(value))){
            stop("There are colnames of value that aren't in cellSummary")}
        value <- value[, colNames]
        by <-
            rep(imageID, unlist(lapply(x[imageID, "cellSummary"], nrow)))
        by <- factor(by, levels = unique(by))
        value <- S4Vectors::split(value, by)
        x <- .putData(x, "cellSummary", value, imageID)
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
    
    factor(rep(rownames(x), unlist(lapply(x$cellSummary, nrow))),
           rownames(x))
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
setReplaceMethod("cellID", "SegmentedCells", function(x, value) {
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
setReplaceMethod("imageCellID", "SegmentedCells", function(x, value) {
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
setReplaceMethod("cellMarks", "SegmentedCells", function(x, imageID = NULL, value) {
    if (is.null(imageID))
        imageID <- rownames(x)
    if (nrow(value) == length(imageID)) {
        x <- .putData(x, "cellMarks", value, imageID)
        return(x)
    }
    
    if (nrow(value) == length(imageID(x))) {
        by <- rep(rownames(x), unlist(lapply(x$cellMarks, nrow)))
        by <- factor(by, levels = unique(by))
        value <- S4Vectors::split(value, by)
        x <- .putData(x, "cellMarks", value, imageID)
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
setReplaceMethod("cellMorph", "SegmentedCells", function(x, imageID = NULL,
                                                     value) {
    if (is.null(imageID))
        imageID <- rownames(x)
    if (nrow(value) == length(imageID)) {
        x <- .putData(x, "cellMorph", value, imageID)
        return(x)
    }
    
    if (nrow(value) == length(imageID(x, imageID))) {
        by <- rep(rownames(x), unlist(lapply(x$cellMorph, nrow)))
        by <- factor(by, levels = unique(by))
        value <- S4Vectors::split(value, by)
        x <- .putData(x, "cellMorph", value, imageID)
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
setReplaceMethod("cellType", "SegmentedCells", function(x, imageID = NULL, value) {
    if (is.null(imageID))
        imageID <- rownames(x)
    loc <- cellSummary(x, imageID = imageID)
    
    if (nrow(loc) != length(value)) {
        stop("There is not enough or too many cellTypes")
    }
    if(!methods::is(value,"factor"))value = factor(value)
    loc$cellType <- value
    
    cellSummary(x, imageID = imageID) <- loc
    x
})



### Get cell type information

#' @export
setGeneric("cellAnnotation", function(x, variable, imageID = NULL)
    standardGeneric("cellAnnotation"))
setMethod("cellAnnotation", "SegmentedCells", function(x, variable, imageID = NULL) {
    if (!is.null(imageID)) {
        x <- x[imageID,]
    }
    cS <- cellSummary(x, bind = TRUE)
    if(!variable%in%colnames(cS))stop("variable not in cellSummary")
    cS[,variable]
})

#' @export
setGeneric("cellAnnotation<-", function(x, variable, imageID = NULL, value)
    standardGeneric("cellAnnotation<-"))
setReplaceMethod("cellAnnotation", "SegmentedCells", function(x, variable, imageID = NULL, value) {
    if (is.null(imageID))
        imageID <- rownames(x)
    loc <- cellSummary(x, bind = TRUE)
    if(length(variable)!=1)stop("Sorry, I can only add one variable at a time currently")
    if(!variable%in%colnames(loc)){
        message(c("Creating variable ", variable))
        loc[,variable] = NA
    }
    
    
    if (sum(loc$imageID%in%imageID) != length(value)) {
        stop("You are trying to put too much or too little into ", variable)
    }
    loc[loc$imageID%in%imageID, variable] <- value
    by <- factor(loc$imageID, rownames(x))
    loc <- loc[, colnames(loc)!="imageID"]
    loc <- S4Vectors::split(loc, by)
    x <- .putData(x, "cellSummary", loc, imageID)
    return(x)
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
        if(dim(pheno)[1]>0){
            return(pheno[imageID(x),])
            }else{
            return(pheno)
            }
        
    } else {
        return(BiocGenerics::do.call("rbind", x$imagePheno))
    }
})


#' @export
setGeneric("imagePheno<-", function(x, imageID = NULL, value)
    standardGeneric("imagePheno<-"))
setReplaceMethod("imagePheno", "SegmentedCells", function(x, imageID = NULL, value) {
    if (is.null(imageID))
        imageID <- as.character(rownames(x))
    value$imageID <- as.character(value$imageID)
    use <- intersect(imageID,value$imageID)
    rownames(value) <- value$imageID
    value <- S4Vectors::split(value[use,], factor(use, levels = use))
    x <- .putData(x, "imagePheno", value, use)
    x[unique(use),]
})



#' @export
setGeneric("filterCells", function(x, select)
    standardGeneric("filterCells"))
setMethod("filterCells", "SegmentedCells", function(x, select) {
    imageID <- imageID(x)
    df <- cellSummary(x, bind = TRUE)
    if (length(select) != nrow(df))
        stop("length of select must equal nrow of SegmentedCells")
    loc <- S4Vectors::split(df[select,], imageID[select])
    x <- .putData(x, "cellSummary", loc, imageID)

    df <- cellMarks(x, bind = TRUE)
    if(length(select) == nrow(df)){
    loc <- S4Vectors::split(df[select,], imageID[select])
    x <- .putData(x, "cellMarks", loc, imageID)
    }
    
    df <- cellMorph(x, bind = TRUE)
    if(length(select) == nrow(df)){
    loc <- S4Vectors::split(df[select,], imageID[select])
    x <- .putData(x, "cellMorph", loc, imageID)
    }
    
    x
})

