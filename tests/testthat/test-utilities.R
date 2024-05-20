test_that("Data formatting works of different formats.", {
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
      SummarizedExperiment::colData(diabetesData) %>%
        as.data.frame() %>%
        select(-cellID, -imageCellID),
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
