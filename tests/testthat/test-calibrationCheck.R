test_that("calibrationCheck() shuffles the condition across patients and reports the share of small p-values", {
  data("diabetesData", package = "spicyR")
  d <- diabetesData[, diabetesData$stage %in% c("Non-diabetic", "Onset")]
  set.seed(1)
  out <- suppressMessages(calibrationCheck(d, condition = "stage", subject = "case", r = 50, nShuffles = 2))
  expect_identical(names(out), c("shuffle", "threshold", "percent"))
  expect_identical(out$shuffle, c(1L, 1L, 2L, 2L))
  expect_identical(out$threshold, c(0.05, 0.01, 0.05, 0.01))
  # the same as shuffling by hand: each patient keeps one condition, moved as a whole
  set.seed(1)
  cd <- SummarizedExperiment::colData(d)
  first <- !duplicated(cd$case)
  lab <- cd$stage[first][sample.int(sum(first))]
  d$shuffledCondition <- lab[match(cd$case, cd$case[first])]
  p <- topPairs(suppressMessages(spicy(d, condition = "shuffledCondition", subject = "case", r = 50)), n = NULL)$p.value
  expect_equal(out$percent[1:2], 100 * c(mean(p < 0.05, na.rm = TRUE), mean(p < 0.01, na.rm = TRUE)))
})

test_that("calibrationCheck() leaves a missing condition missing", {
  data("diabetesData", package = "spicyR")
  d <- diabetesData[, diabetesData$stage %in% c("Non-diabetic", "Onset")]
  d$stage <- as.character(d$stage)
  d$stage[d$case == unique(d$case)[1]] <- NA   # one patient with no condition
  set.seed(2)
  expect_no_error(out <- suppressMessages(calibrationCheck(d, condition = "stage", subject = "case", r = 50,
                                                           nShuffles = 1)))
  expect_error(calibrationCheck(d, condition = "nothere", subject = "case", r = 50), "columns not found")
  expect_error(calibrationCheck(tempdir(), condition = "stage"), "in memory")
})
