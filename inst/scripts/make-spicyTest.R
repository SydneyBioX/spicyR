## Makes data/spicyTest.rda: the image-level test of spicyR on diabetesData.
## Run from the package root, with this version of spicyR installed:
##   Rscript inst/scripts/make-spicyTest.R
library(spicyR)
data("diabetesData", package = "spicyR")
spicyTest <- spicy(diabetesData, condition = "stage", subject = "case", method = "image")
save(spicyTest, file = file.path("data", "spicyTest.rda"), compress = "xz")
