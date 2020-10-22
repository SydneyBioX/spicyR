#' Diabetes IMC data 
#'
#' This is a subset of the Damond et al 2019 imaging mass cytometry dataset. The data 
#' contains cells in the pancreatic islets of individuals with early onset 
#' diabetes and healthy controls. 
#' The object contains single-cell data of 160 images from 8 subjects, 
#' with 20 images per subject.
#'
#' @format diabetesDataSub a SegmentedCells object
#' @aliases 
#' diabetesDataSub
"diabetesDataSub"



#' Results from spicy for diabetesDataSub
#'
#' Results from the call:
#' spicyTest <- spicy(diabetesDataSub, 
#'                    condition = "condition", 
#'                    subject = "subject")
#'
#' @format spicyTest a spicy object
#' @aliases 
#' spicyTest
"spicyTest"



#' Results from spicy with bootstrap for diabetesDataSub
#'
#' Results from the call:
#' spicyTestBootstrap <- spicy(diabetesDataSub, 
#'                            condition = "condition", 
#'                            subject = "subject", 
#'                            nsim = 199)
#'
#' @format spicyTestBootstrap a spicy object
#' @aliases 
#' spicyTestBootstrap
"spicyTestBootstrap"
