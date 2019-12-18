### Convert something into a ImagingCytometryExperiment

as.ImagingCytometryExperiment <- function(x){
  if(any(grepl('Intensity_mean',colnames(cells)))){
    
    df = DataFrame(row.names = unique(x$imageID))
    location <- split(DataFrame(x[,c('cellID','cellType','x','y')]), x$imageID)
    df$location <- location
    intensity <- x[,grep('Intensity_mean',colnames(x))]
    colnames(intensity) <- gsub('Intensity_mean_', '', intensity)
    df$intensity = split(DataFrame(intensity), x$imageID)
    morphology <- x[,grep('Area',colnames(x))]
    colnames(morphology) <- gsub('Area_', '', morphology)
    df$morphology = split(DataFrame(morphology), x$imageID)
    df$clinical = split(DataFrame(), x$imageID)
    return(df)
  }
}


### Some functions to get and make location data

location <- function(x, image = NULL, bind = FALSE){
  if(!is.null(image)){
    x = x[image,]
  }
  if(bind==FALSE){
    return(x$location)
  }
  if(bind==TRUE){
    return(do.call('rbind',x$location))
  }
}


'location<-' <- function(x, value){
  if(nrow(value)==nrow(x)){
  x$location <- value
  return(x)
  }
  
  if(nrow(value)==length(imageID(x))){
    x$location <- split(value,rep(rownames(x),unlist(lapply(x$location,nrow))))
    return(x)
  }
    
}


### Get imageIDs for each cell, not sure if this should also report rownames(df)

imageID <- function(x, image = NULL){
  if(!is.null(image)){
    x = x[image,]
  }
  rep(rownames(x),unlist(lapply(x$location,nrow)))
}

### Get cellTypes for each cell

cellType <- function(x, image = NULL){
  if(!is.null(image)){
    x = x[image,]
  }
  do.call('rbind',x$location)$cellType
}

'cellType<-' <- function(x, value){
  
  loc <- location(x, bind = TRUE)
  
  if(nrow(loc)!=nrow(x)){
    stop('There is not enough or too many cell types')
  }
  
  loc$cellType <- value
  
  location(x) <- loc
  
}


### Get cellID

cellID <- function(x, image = NULL){
  if(!is.null(image)){
    x = x[image,]
  }
  do.call('rbind',x$location)$cellID
}


'cellID<-' <- function(x, value){
  
  loc <- location(x, bind = TRUE)
  
  if(nrow(loc)!=nrow(x)){
    stop('There is not enough or too many cellIDs')
  }
  
  loc$cellID <- value
  
  location(x) <- loc

}



### Get uniqueCellID

uniqueCellID <- function(x, image = NULL){
  if(!is.null(image)){
    x = x[image,]
  }
  do.call('rbind',x$location)$uniqueCellID
}


makeUniqueCellID <- function(x){
 loc <- location(x,bind=TRUE)
 loc$uniqueCellID <- paste('cell', seq_len(nrow(loc)),sep = '')
 location(x) <- loc
 x
}




#### I can access and add clinical data to the object

clinical <- function(x, image = NULL, bind = FALSE){
  if(!is.null(image)){
    x = x[image,]
  }
    do.call('rbind',x$clinical)
  }



'clinical<-' <- function(x, value, image = NULL){
  if(is.null(image)) image <- rownames(x)
  use <- intersect(value$imageID,image)
  x <- x[image,] 
  x$clinical <- split(value,image)
  x
}





# library(S4Vectors)
# 
# 
# 
# ### Something that resembles cellProfiler data
# 
# cells <- data.frame(row.names = 1:10)
# cells$cellID <- paste('cell',1:10,sep = '')
# cells$imageID <- paste('image',rep(1:2,c(4,6)),sep = '')
# cells$x <- 1:10
# cells$y <- 10:1
# cells$cellType <- paste('cellType',rep(1:2,5),sep = '')
# cells$Area_round <- 1
# cells$Area_diameter <- 2
# cells$Intensity_mean_CD8 <- 3
# cells$Intensity_mean_CD4 <- 11:20
# 
# 
# 
# ### Lets make a CellImagingExperiment... obviously we can make this a class
# 
# df <- as.ImagingCytometryExperiment(cells)
# 
# 
# 
# clinicalData <- DataFrame(sex = c('M',"F"), height = c(120,210), weight = c(30, 140), imageID = c('image1','image2'))
# clinical(df) <- clinicalData
# 
# ### I can get the locations and set them again
# z <- location(df)
# z
# 
# ### If I give it a list of DataFrames
# location(df) <- z
# 
# 
# z <- location(df, bind = TRUE)
# z
# ### If I give it a single DataFrame
# location(df) <- z
# 
# 
# ### I can pull out information for specific images
# df['image1',]
# location(df['image1',])
# location(df,'image1')
# location(df,c('image1','image2'),TRUE)
# 
# cellType(df)
# 
# 
# 
# ### We could also add an image(df) which contained stacks maybe even image(df,'image1', 'CD8')
