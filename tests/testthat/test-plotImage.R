test_that("plotImage() shows every cell type without a pair, and the pair with one", {
  data("diabetesData", package = "spicyR")
  p <- plotImage(diabetesData, "A09")
  expect_s3_class(p, "ggplot")
  expect_identical(sort(unique(p$data$cellType)),
                   sort(unique(as.character(diabetesData$cellType[diabetesData$imageID == "A09"]))))
  # a cell type has the same colour in every image
  pal <- function(im) { q <- plotImage(diabetesData, im); g <- ggplot2::ggplot_build(q)$data[[1]]
                        unique(data.frame(type = q$data$cellType, colour = g$colour)) }
  both <- merge(pal("A09"), pal(unique(diabetesData$imageID)[2]), by = "type")
  expect_gt(nrow(both), 1)
  expect_identical(both$colour.x, both$colour.y)
  expect_s3_class(plotImage(diabetesData, "A09", from = "acinar", to = "alpha", r = 50), "ggplot")
  expect_error(plotImage(diabetesData, "A09", from = "acinar"), "both")
  expect_error(plotImage(diabetesData, "A09", r = 50), "circles")
})

test_that("plotImage() draws several images as panels, in the order given", {
  data("diabetesData", package = "spicyR")
  ims <- c("A11", "A09")
  p <- plotImage(diabetesData, ims, from = "acinar", to = "alpha", r = 50)
  expect_s3_class(p$facet, "FacetWrap")
  expect_identical(levels(p$layers[[1]]$data$image), ims)
  expect_equal(nrow(p$layers[[1]]$data), sum(diabetesData$imageID %in% ims))
  expect_s3_class(plotImage(diabetesData, ims)$facet, "FacetWrap")
  expect_error(plotImage(diabetesData, c("A09", "nothere")), "nothere")
})
