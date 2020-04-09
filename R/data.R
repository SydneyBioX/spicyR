#' Melanoma responders/non-responders dataset.
#'
#' A dataset containing cell-types of 27 melanoma patients that did and did not
#' respond to treatment. Each subject has 5 images. The data has already been 
#' summarised into 9 cell-types.
#'
#' @format melanomaResponders a SegmentedCells object
#' @aliases 
#' melanomaResponders
"melanomaResponders"



#' Results from spicy for melanomaResponders data
#'
#' Results from the call:
#' spicyTest <- spicy(melanomaResponders, 
#'                    condition = "condition", 
#'                    subject = "subject")
#'
#' @format spicyTest a spicy object
#' @aliases 
#' spicyTest
"spicyTest"



#' Results from spicy with bootstrap for melanomaResponders data
#'
#' Results from the call:
#' spicyTestBootstrap <- spicy(melanomaResponders, 
#'                            condition = "condition", 
#'                            subject = "subject", 
#'                            nsim = 199)
#'
#' @format spicyTestBootstrap a spicy object
#' @aliases 
#' spicyTestBootstrap
"spicyTestBootstrap"
