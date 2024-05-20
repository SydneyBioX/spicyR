#' Diabetes IMC data in SCE format.
#'
#' This is a subset of the Damond et al 2019 imaging mass cytometry dataset. The data 
#' contains cells in the pancreatic islets of individuals with early onset 
#' diabetes and healthy controls. 
#' The object contains single-cell data of 160 images from 8 subjects, 
#' with 20 images per subject.
#'
#' Converted into a SingleCellExperiment format.
#'
#' @format diabetesData_SCE a SingleCellExperiment object
#' @usage data("diabetesData")
#' @aliases 
#' diabetesData
"diabetesData"


#' Results from spicy for diabetesData
#'
#' Results from the call:
#' spicyTest <- spicy(diabetesData, 
#'                    condition = "condition", 
#'                    subject = "subject")
#'
#' @format spicyTest a spicy object
#' @usage data("spicyTest")
#' @aliases 
#' spicyTest
"spicyTest"
