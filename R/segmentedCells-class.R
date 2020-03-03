#' The segmentedCells class
#' 
#' The segmentedCells S4 class is for storing data from segmented 
#' imaging cytometry and spatial omics data. It extends DataFrame and defines 
#' methods that take advantage of DataFrame nesting to represent elements of 
#' cell-based experiments with spatial orientation that are commonly 
#' encountered. This object is able to store information on a cell's spatial 
#' location, cellType, morphology, intensity of gene/protein marks as well as 
#' image level phenotype information.
#'
#' @param cellData A data frame that contains at least the columns x and y giving the location of each cell.
#' @param cellProfiler A logical indicating that cellData is in a format similar to what cellProfiler outputs.
#' @param spatialCoords The column names corresponding to spatial coordinates. eg. x, y, z...
#' @param cellTypeString The name of the column that contains cell type calls.
#' @param intensityString A string which can be used to identify the columns which contain marker intensities. (This needs to be extended to take the column names themselves.)
#' @param morphologyString A string which can be used to identify the columns which contains morphology information.
#' @param cellIDString The column name for cellID.
#' @param imageCellIDString The column name for imageCellID.
#' @param imageIDString The column name for imageIDString
#' 
#' @return A segmentedCells object
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
#' cellExp <- segmentedCells(cells, cellProfiler = TRUE)
#' 
#' ### Cluster cell types
#' intensities <- intensity(cellExp)
#' kM <- kmeans(intensities,2)
#' cellType(cellExp) <- paste('cluster',kM$cluster, sep = '')
#' location(cellExp)
#' 
#' @aliases 
#' segmentedCells
#' segmentedCells,segmentedCells-method
#' 
#' @export
#' @rdname segmentedCells
#' @importFrom methods new
#' @importFrom S4Vectors DataFrame split
segmentedCells <- function(cellData, cellProfiler = FALSE, spatialCoords = NULL, 
    cellTypeString = NULL, intensityString = NULL, morphologyString = NULL, phenotypeString = NULL, cellIDString = NULL, 
    imageCellIDString = NULL, imageIDString = NULL) {
    
    ### Check variable names
    
    if (!is.null(cellIDString)) {
        if (!cellIDString %in% colnames(cellData)) 
            stop("cellIDString is not a column name of cellData")
        cellData$cellID <- cellData[, cellIDString]
    }
    
    if (!is.null(imageIDString)) {
        if (!imageIDString %in% colnames(cellData)) 
            stop("imageIDString is not a column name of cellData")
        cellData$imageCellID <- cellData[, imageIDString]
    }
    
    if (!is.null(imageIDString)) {
        if (!imageIDString %in% colnames(cellData)) 
            stop("imageIDString is not a column name of cellData")
        
        cellData$imageID <- cellData[, imageIDString]
    }
    
    if (is.null(cellTypeString) & "cellType" %in% colnames(cellData)) {
        cellTypeString <- "cellType"
    }
    
    if (!is.null(cellTypeString)) {
        if (!cellTypeString %in% colnames(cellData)) 
            stop("cellTypeString is not a column name of cellData")
        
        cellData$cellType <- cellData[, cellTypeString]
    }
    
    if (!is.null(intensityString)) {
        if (length(intensityString) > 1) 
            stop("intensityString needs to be NULL or length 1")
        if (length(grep(intensityString, colnames(cellData))) == 0) 
            stop("intensityString is not in column names of cellData")
        
    }
    
    
    if (!is.null(morphologyString)) {
        if (length(morphologyString) > 1) 
            stop("morphologyString needs to be NULL or length 1")
        if (length(grep(morphologyString, colnames(cellData))) == 0) 
            stop("morphologyString is not in column names of cellData")
    }
    
    
    
    if (!is.null(phenotypeString)) {
        if (length(phenotypeString) > 1) 
            stop("phenotypeString needs to be NULL or length 1")
        if (length(grep(phenotypeString, colnames(cellData))) == 0) 
            stop("phenotypeString is not in column names of cellData")
    }
    
    
    if (is.null(phenotypeString)) {
        if (length(grep("phenotype_", colnames(cellData))) > 0) 
            phenotypeString = "phenotype_"
    }
    
    ### Format variable names if not from cellProfiler output
    
    if (!cellProfiler) {
        
        if (is.null(intensityString)) {
            if (length(grep("intensity_", colnames(cellData))) == 0) 
                intensityString = "intensity_"
        }
        
        
        if (is.null(morphologyString)) {
            if (length(grep("morphology_", colnames(cellData))) == 0) 
                intensityString = "morphology_"
        }
        
        
        
        if (is.null(cellData$cellID)) {
            cat("There is no cellID. I'll create these", "\n")
            cellData$cellID <- paste("cell", seq_len(nrow(cellData)), sep = "_")
        }
        
        if (!is.null(spatialCoords)) {
            cellData$x <- cellData[,spatialCoords[1]]
        }
        
        if (!is.null(spatialCoords)) {
            cellData$y <- cellData[,spatialCoords[2]]
        }
        
        spatialCoords <- c("x", "y")
        
        if (is.null(cellData$x)) {
            stop("You need to include a 'x' column in the data.frame")
        }
        
        if (is.null(cellData$y)) {
            stop("You need to include a 'y' column in the data.frame")
        }
        
        
        if (is.null(cellData$imageCellID)) {
            cat("There is no image specific imageCellID. I'll create these", "\n")
            cellData$imageCellID <- paste("cell", seq_len(nrow(cellData)), sep = "_")
        }
        if (length(cellData$imageCellID) != nrow(cellData)) 
            stop("The number of rows in cells does not equal the number of imageCellIDs")
        
        if (is.null(cellData$imageID)) {
            cat("There is no imageID. I'll assume this is only one image and create an arbitrary imageID", 
                "\n")
            cellData$imageID <- "image1"
        }
        cellData$imageID <- as.factor(cellData$imageID)
        cellData$cellID <- as.character(cellData$cellID)
    }
    
    ### Format variable names if cellProfiler output
    
    if (cellProfiler) {
        
        cellData$imageID <- as.factor(cellData$ImageNumber)
        cellData$cellID <- paste('cell',1:nrow(cellData),sep='_')
        cellData$imageCellID <- cellData$ObjectNumber
        
        if (is.null(spatialCoords)) {
            spatialCoords <- c("AreaShape_Center_X", "AreaShape_Center_Y")
        }
        
        if (!is.null(spatialCoords)) {
            cellData$x <- cellData[,spatialCoords[1]]
        }
        
        if (!is.null(spatialCoords)) {
            cellData$y <- cellData[,spatialCoords[2]]
        }
        

        if (is.null(intensityString) & any(grepl("Intensity_Mean_", colnames(cellData)))) {
            intensityString <- "Intensity_Mean_"
        }
        
        if (is.null(morphologyString) & any(grepl("AreaShape_", colnames(cellData)))) {
            morphologyString <- "AreaShape_"
        }
        
        if (any(grepl(morphologyString,spatialCoords))) {
            cn <- spatialCoords[grep(morphologyString,spatialCoords)]
            cellData <- cellData[,!colnames(cellData)%in%cn]
        }
        
        spatialCoords <- c("x", "y")
        
    }
    
    
    ### Create location information
    
    df <- DataFrame(row.names = unique(cellData$imageID))
    
    if (!is.null(cellTypeString)) {
        cellData$cellType <- cellData[, cellTypeString]
        location <- S4Vectors::split(DataFrame(cellData[, c("cellID", "imageCellID", 
            spatialCoords, "cellType")]), cellData$imageID)
    } else {
        cellData$cellType <- NA
        location <- S4Vectors::split(DataFrame(cellData[, c("cellID", "imageCellID", 
            spatialCoords, "cellType")]), cellData$imageID)
    }
    
    df$location <- location
    
    
    ### Create intensity and morphology information
    
    df$intensity <- S4Vectors::split(DataFrame(), cellData$imageID)
    df$morphology <- S4Vectors::split(DataFrame(), cellData$imageID)
    
    if (!is.null(intensityString)) {
        intensity <- cellData[, grep(intensityString, colnames(cellData))]
        colnames(intensity) <- gsub(intensityString, "", colnames(intensity))
        df$intensity <- S4Vectors::split(DataFrame(intensity), cellData$imageID)
    }
    
    if (!is.null(morphologyString)) {
        morphology <- cellData[, grep(morphologyString, colnames(cellData))]
        colnames(morphology) <- gsub(morphologyString, "", colnames(morphology))
        df$morphology <- S4Vectors::split(DataFrame(morphology), cellData$imageID)
    }
    
    
    ### Create columns in DataFrame for storing phenotype information and potentially
    ### images and masks.
    
    df$phenotype <- S4Vectors::split(DataFrame(), cellData$imageID)
    
    if(!is.null(phenotypeString)){
        phenotype <- cellData[, grep(phenotypeString, colnames(cellData))]
        colnames(phenotype) <- gsub(phenotypeString, "", colnames(phenotype))
        phenotype <- cbind(imageID = cellData$imageID, phenotype)
        phenotype <- unique(phenotype)
        phenotype <- S4Vectors::split(DataFrame(phenotype), phenotype$imageID)
        df$phenotype <- phenotype[rownames(df)]
    }
    
    df$images <- S4Vectors::split(DataFrame(), cellData$imageID)
    df$masks <- S4Vectors::split(DataFrame(), cellData$imageID)
    
    ### Create segmentedCells object.
    
    df <- new("segmentedCells", df)
    df
}


