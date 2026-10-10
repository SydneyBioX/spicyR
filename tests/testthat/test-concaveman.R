test_that("the C++ concave hull gives concaveman's polygons exactly", {
  # polygons from concaveman::concaveman() 1.2.0, for the points of a concave window around 150 cells
  ref <- readRDS(test_path("concaveman_reference.rds"))
  for (case in names(ref)) {
    r <- ref[[case]]
    expect_identical(.concaveHull(r$points[, 1], r$points[, 2], r$concavity, r$length_threshold), r$polygon,
                     label = case)
  }
})
