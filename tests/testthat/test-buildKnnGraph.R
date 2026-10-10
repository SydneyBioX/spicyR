test_that("buildKnnGraph() gives each cell its k nearest neighbours in its image, as nnwhich()", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("spatstat.geom")
  cells <- diabetesData[, diabetesData$imageID %in% c("A09", "A11", "A16")]
  # interleave the images, so the graph must be mapped back to the original cell order
  cells <- cells[, order(seq_len(ncol(cells)) %% 7)]
  g <- SingleCellExperiment::colPair(buildKnnGraph(cells, k = 5, cores = 2), "knn_interaction_graph")

  expected <- do.call(rbind, lapply(split(seq_len(ncol(cells)), cells$imageID), function(i) {
    nn <- spatstat.geom::nnwhich(cells$x[i], cells$y[i], k = 1:5)
    cbind(rep(i, each = 5), i[t(nn)])
  }))
  expected <- expected[order(expected[, 1], expected[, 2]), ]
  expect_equal(cbind(S4Vectors::from(g), S4Vectors::to(g)), unname(expected))
})

test_that("cells in images with at most k cells get no neighbours", {
  skip_if_not_installed("SingleCellExperiment")
  a11 <- which(diabetesData$imageID == "A11")[1:5]
  cells <- diabetesData[, c(which(diabetesData$imageID == "A09"), a11)]
  g <- SingleCellExperiment::colPair(buildKnnGraph(cells, k = 5), "knn_interaction_graph")
  expect_false(any(cells$imageID[S4Vectors::from(g)] == "A11"))
  expect_equal(length(g), 5 * sum(cells$imageID == "A09"))
})
