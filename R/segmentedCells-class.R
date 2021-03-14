#' The SegmentedCells class
#'
#' The SegmentedCells S4 class is for storing data from segmented
#' imaging cytometry and spatial omics data. It extends DataFrame and defines
#' methods that take advantage of DataFrame nesting to represent elements of
#' cell-based experiments with spatial orientation that are commonly
#' encountered. This object is able to store information on a cell's spatial
#' location, cellType, morphology, intensity of gene/protein markers as well as
#' image level phenotype information.
#'
#' @param cellData A data frame that contains at least the columns x and y giving 
#' the location coordinates of each cell.
#' @param cellProfiler A logical indicating that cellData is in a format similar 
#' to what cellProfiler outputs.
#' @param spatialCoords The column names corresponding to spatial coordinates. 
#' eg. x, y, z...
#' @param cellTypeString The name of the column that contains cell type calls.
#' @param intensityString A string which can be used to identify the columns 
#' which contain marker intensities. (This needs to be extended to take the 
#' column names themselves.)
#' @param morphologyString A string which can be used to identify the columns 
#' which contains morphology information.
#' @param cellIDString The column name for cellID.
#' @param cellAnnotations A vector of variables that provide additional 
#' annotation of a cell.
#' @param imageCellIDString The column name for imageCellID.
#' @param imageIDString The column name for imageIDString.
#' @param phenotypeString A string which can be used to identify the columns 
#' which contains phenotype information.
#'
#' @return A SegmentedCells object
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
#' cells$ImageNumber <- rep(seq_len(2),c(n/2,n/2))
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
#' cellSummary(cellExp)
#'
#' @aliases
#' SegmentedCells
#' SegmentedCells,SegmentedCells-method
#'
#' @export
#' @rdname SegmentedCells
#' @importFrom methods new
#' @importFrom S4Vectors DataFrame split
#' @importFrom IRanges SplitDataFrameList
SegmentedCells <-
    function(cellData,
             cellProfiler = FALSE,
             spatialCoords = c("x", "y"),
             cellTypeString = "cellType",
             intensityString = "intensity_",
             morphologyString = "morphology_",
             phenotypeString = "phenotype_",
             cellIDString = "cellID",
             cellAnnotations = NULL,
             imageCellIDString = "imageCellID",
             imageIDString = "imageID") {
        
        if (cellTypeString == "cellType" & !"cellType" %in% colnames(cellData)) {
            message("There is no cellType column, setting to NA")
            cellData$cellType = NA
        }
        if (!is.null(cellTypeString)) {
            if (!cellTypeString %in% colnames(cellData)) 
                stop("cellTypeString is not a column name of cellData")
            cellData$cellType <- cellData[, cellTypeString]
        }
        
        
        if(!is(cellData$cellType,"factor")){
            cellData$cellType = factor(cellData$cellType)
        }
        
        if (!is.null(cellAnnotations)) {
            if (any(!cellAnnotations %in% colnames(cellData))) 
                stop("At least one element of cellAnnotations is not a column name of cellData")
        }
        
        if (!any(grepl(intensityString, colnames(cellData))) & intensityString != 
            "intensity_") 
            stop("intensityString is not in column names of cellData")
        if (!any(grepl(morphologyString, colnames(cellData))) & morphologyString != 
            "morphology_") 
            stop("morphologyString is not in column names of cellData")
        if (!any(grepl(phenotypeString, colnames(cellData))) & phenotypeString != 
            "phenotype_") 
            stop("phenotypeString is not in column names of cellData")
        if (!cellProfiler) {
            if (any(!spatialCoords %in% colnames(cellData))) 
                stop("spatialCoords are not column names of cellData")
            cellData$x <- cellData[, spatialCoords[1]]
            cellData$y <- cellData[, spatialCoords[2]]
            spatialCoords <- c("x", "y")
            if (cellIDString == "cellID" & !"cellID" %in% colnames(cellData)) {
                message("There is no cellID. I'll create these", 
                        "\n")
                cellData$cellID <- paste("cell", seq_len(nrow(cellData)), 
                                         sep = "_")
            }
            if (!is.null(cellIDString)) {
                if (!cellIDString %in% colnames(cellData)) 
                    stop("cellIDString is not a column name of cellData")
                if (length(unique(cellData$cellID)) != length(cellData$cellID)) 
                    stop("Your cellIDs are not unique to each cell ")
                cellData$cellID <- cellData[, cellIDString]
            }
            if (imageCellIDString != "imageCellID") {
                if (!imageIDString %in% colnames(cellData)) 
                    stop("imageIDString is not a column name of cellData")
                cellData$imageCellID <- cellData[, imageIDString]
            }
            if (imageIDString != "imageID") {
                if (!imageIDString %in% colnames(cellData)) 
                    stop("imageIDString is not a column name of cellData")
                cellData$imageID <- cellData[, imageIDString]
            }
            if (is.null(cellData$imageCellID)) {
                message("There is no image specific imageCellID. I'll create these", 
                        "\n")
                cellData$imageCellID <- paste("cell", seq_len(nrow(cellData)), 
                                              sep = "_")
            }
            if (length(cellData$imageCellID) != nrow(cellData)) 
                stop("The number of rows in cells does not equal the number of imageCellIDs")
            if (is.null(cellData$imageID)) {
                message("There is no imageID. I'll assume this is only one image and create an arbitrary imageID", 
                        "\n")
                cellData$imageID <- "image1"
            }
            if(!is(cellData$imageID, "factor")){
                cellData$imageID <- factor(cellData$imageID, unique(cellData$imageID))
            }
            cellData$cellID <- as.character(cellData$cellID)
        }
        if (cellProfiler) {
            if (imageIDString == "imageID" & "ImageNumber" %in% colnames(cellData) & 
                !"imageID" %in% colnames(cellData)) 
                cellData$imageID <- as.factor(cellData$ImageNumber)
            if (imageIDString != "imageID") {
                if (!imageIDString %in% colnames(cellData)) 
                    stop("imageIDString is not a column name of cellData")
                cellData$imageID <- cellData[, imageIDString]
            }
            if (is.null(cellData$imageID)) {
                message("There is no imageID. I'll assume this is only one image and create an arbitrary imageID", 
                        "\n")
                cellData$imageID <- "image1"
            }
            cellData$cellID <- paste("cell", seq_len(nrow(cellData)), 
                                     sep = "_")
            if (!is.null(cellIDString)) {
                if (!cellIDString %in% colnames(cellData)) 
                    stop("cellIDString is not a column name of cellData")
                if (length(unique(cellData$cellID)) != length(cellData$cellID)) 
                    stop("Your cellIDs are not unique to each cell ")
                cellData$cellID <- cellData[, cellIDString]
            }
            cellData$imageID <- as.factor(cellData$imageID)
            cellData$cellID <- as.character(cellData$cellID)
            if (imageCellIDString == "imageCellID" & "ObjectNumber" %in% 
                colnames(cellData) & !"imageCellID" %in% colnames(cellData)) 
                cellData$imageCellID <- cellData$ObjectNumber
            cellData$imageCellID <- paste("cell", seq_len(nrow(cellData)), 
                                          sep = "_")
            if (imageCellIDString != "imageCellID") {
                if (!imageCellIDString %in% colnames(cellData)) 
                    stop("imageCellIDString is not a column name of cellData")
                cellData$imagCelleID <- cellData[, imageIDString]
            }
            if (all(spatialCoords %in% c("x", "y")) & all(c("AreaShape_Center_X", 
                                                            "AreaShape_Center_Y") %in% colnames(cellData)) & 
                !all(c("x", "y") %in% colnames(cellData))) {
                spatialCoords <- c("AreaShape_Center_X", "AreaShape_Center_Y")
            }
            if (any(!spatialCoords %in% colnames(cellData))) 
                stop("spatialCoords are not column names of cellData")
            cellData$x <- cellData[, spatialCoords[1]]
            cellData$y <- cellData[, spatialCoords[2]]
            if (intensityString == "intensity_" & !any(grepl("intensity_", 
                                                             colnames(cellData))) & any(grepl("Intensity_Mean_", 
                                                                                              colnames(cellData)))) {
                intensityString <- "Intensity_Mean_"
            }
            if (morphologyString == "morphology_" & !any(grepl("morphology_", 
                                                               colnames(cellData))) & any(grepl("AreaShape_", colnames(cellData)))) {
                morphologyString <- "AreaShape_"
            }
            if (any(grepl(morphologyString, spatialCoords))) {
                cn <- spatialCoords[grep(morphologyString, spatialCoords)]
                cellData <- cellData[, !colnames(cellData) %in% cn]
            }
            spatialCoords <- c("x", "y")
        }
        cellData$imageID <- droplevels(cellData$imageID)
        df <- DataFrame(row.names = levels(cellData$imageID))
            cellData$cellType <- cellData[, cellTypeString]
            cellSummaryCols <- c("cellID", "imageCellID", spatialCoords, "cellType", cellAnnotations)
            cellSummary <- S4Vectors::split(DataFrame(cellData[,cellSummaryCols]), 
                                            cellData$imageID)
        df$cellSummary <- cellSummary[rownames(df),]
        df$cellMarks <- S4Vectors::split(DataFrame(), cellData$imageID)
        df$cellMorph <- S4Vectors::split(DataFrame(), cellData$imageID)
        if (any(grepl(intensityString, colnames(cellData)))) {
            markers <- cellData[, grep(intensityString, colnames(cellData))]
            colnames(markers) <- gsub(intensityString, "", colnames(markers))
            df$cellMarks <- S4Vectors::split(DataFrame(markers), 
                                             cellData$imageID)[rownames(df)]
        }
        if (any(grepl(morphologyString, colnames(cellData)))) {
            morphology <- cellData[, grep(morphologyString, colnames(cellData))]
            colnames(morphology) <- gsub(morphologyString, "", colnames(morphology))
            df$cellMorph <- S4Vectors::split(DataFrame(morphology), 
                                             cellData$imageID)[rownames(df)]
        }
        df$imagePheno <- S4Vectors::split(DataFrame(), cellData$imageID)
        if (any(grepl(phenotypeString, colnames(cellData)))) {
            phenotype <- cellData[, grep(phenotypeString, colnames(cellData))]
            colnames(phenotype) <- gsub(phenotypeString, "", colnames(phenotype))
            phenotype <- cbind(imageID = cellData$imageID, phenotype)
            phenotype <- unique(phenotype)
            phenotype <- S4Vectors::split(DataFrame(phenotype), phenotype$imageID)
            df$imagePheno <- phenotype[rownames(df)]
        }
        
        df$images <- S4Vectors::split(DataFrame(), cellData$imageID)[rownames(df)]
        df$masks <- S4Vectors::split(DataFrame(), cellData$imageID)[rownames(df)]
        df <- new("SegmentedCells", df)
        df
    }



#' as.data.frame
#'
#' Function to coerce a SegmentedCells object to a data frame.
#'
#' @param x A SegmentedCells object.
#' @param ... Other arguments.
#' 
#' @return A data.frame
#' 
#' ## Generate toy data
#' set.seed(51773)
#' x <- round(c(runif(200),runif(200)+1,runif(200)+2,runif(200)+3,
#'              runif(200)+3,runif(200)+2,runif(200)+1,runif(200)),4)
#'              y <- round(c(runif(200),runif(200)+1,runif(200)+2,runif(200)+3,
#'                           runif(200),runif(200)+1,runif(200)+2,runif(200)+3),4)
#' cellType <- factor(paste('c',rep(rep(c(1:2),rep(200,2)),4),sep = ''))
#' imageID <- rep(c('s1', 's2'),c(800,800))
#' cells <- data.frame(x, y, cellType, imageID)
#' 
#' ## Store data in SegmentedCells object
#' cellExp <- SegmentedCells(cells, cellTypeString = 'cellType')
#' 
#' ## Generate LISA
#' cellsDF <- as.data.frame(cellExp)
#' 
#' 
#' @return \code{NULL}
#'
#' @rdname as.data.frame.SegmentedCells
#' @method as.data.frame SegmentedCells
#' @export
as.data.frame.SegmentedCells <- function(x, ...) {
    loc <- cellSummary(x)
    int <- cellMarks(x)
    morph <- cellMorph(x)
    pheno <- imagePheno(x, expand = TRUE)
    pheno <- pheno[!colnames(pheno) %in% "imageID"]
#    if("region"%in%names(x))loc <- lisaClust::region(x, annot = TRUE)
    
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
    
    if (prod(dim(loc)) > 0) {
        return(as.data.frame(loc))
    }
    
    NA
}


setAs("SegmentedCells", "data.frame", function(from) {
    from <- as.data.frame.SegmentedCells(from)
    new("SegmentedCells", from)
})