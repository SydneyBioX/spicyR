test_that("Data formatting works of different formats.", {
  skip_if_not_installed("SpatialExperiment")
  diabetesData_SPE <- SpatialExperiment::SpatialExperiment(
    diabetesData,
    colData = SummarizedExperiment::colData(diabetesData)
  )

  SpatialExperiment::spatialCoords(diabetesData_SPE) <- data.frame(
    SummarizedExperiment::colData(diabetesData_SPE)$x,
    SummarizedExperiment::colData(diabetesData_SPE)$y
  ) |> as.matrix()
  SpatialExperiment::spatialCoordsNames(diabetesData_SPE) <- c("x", "y")

  expect_no_error(
    .format_data(diabetesData_SPE, "imageID", "cellType", c("x", "y"), FALSE)
  )
  expect_no_error(
    .format_data(
      SummarizedExperiment::colData(diabetesData),
      "imageID", "cellType", c("x", "y"), FALSE
    )
  )
  expect_no_error(
    .format_data(
      subset(as.data.frame(SummarizedExperiment::colData(diabetesData)),
             select = -c(cellID, imageCellID)),
      "imageID", "cellType", c("x", "y"), FALSE
    )
  )
})

test_that(
  "Test .format_data input validation",
  {
    expect_error(
      .format_data(
        diabetesData, "imageIDs", "cellType", c("x", "y"), FALSE
      )
    )
    expect_error(
      .format_data(
        diabetesData, "imageID", "cellType", c("x", "b"), FALSE
      )
    )
  }
)

test_that(".par_lapply() and .par_mapply() give the serial results and leave BLAS threads as they were", {
    skip_on_os("windows")
    X <- setNames(as.list(1:20), letters[1:20])
    blasThreads <- RhpcBLASctl::blas_get_num_procs()
    expect_identical(.par_lapply(X, function(x, y) x * y, y = 3, cores = 2), lapply(X, function(x, y) x * y, y = 3))
    expect_identical(.par_lapply(X, sqrt, cores = 2), lapply(X, sqrt))
    expect_identical(
        .par_mapply(function(a, b, k) a + b + k, a = X, b = rev(X), MoreArgs = list(k = 1), cores = 2),
        mapply(function(a, b, k) a + b + k, a = X, b = rev(X), MoreArgs = list(k = 1), SIMPLIFY = FALSE)
    )
    expect_identical(RhpcBLASctl::blas_get_num_procs(), blasThreads)
    expect_error(.par_lapply(1:4, function(x) if (x == 3) stop("bad task") else x, cores = 2), "bad task")
})
