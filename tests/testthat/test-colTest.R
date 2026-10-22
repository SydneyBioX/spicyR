test_that("colTest(subject =) pools each patient's images, as imageID = <patient column> did", {
  data("diabetesData", package = "spicyR")
  d <- diabetesData[, diabetesData$stage %in% c("Long-duration", "Onset")]
  bySubject <- colTest(d, condition = "stage", feature = "cellType", subject = "case", type = "wilcox")
  byImageArg <- colTest(d, condition = "stage", feature = "cellType", imageID = "case", type = "wilcox")
  expect_identical(bySubject, byImageArg)
  # each patient is one unit: proportions of the patient's cells, all images pooled
  props <- getProp(d, feature = "cellType", imageID = "case")
  expect_identical(nrow(props), length(unique(d$case)))
})

test_that("colTest() reports the median in each group with a Wilcoxon test", {
  data("diabetesData", package = "spicyR")
  d <- diabetesData[, diabetesData$stage %in% c("Long-duration", "Onset")]
  res <- colTest(d, condition = "stage", feature = "cellType", type = "wilcox")
  props <- getProp(d, feature = "cellType")
  stage <- unique(as.data.frame(SummarizedExperiment::colData(d))[, c("imageID", "stage")])
  g <- factor(setNames(as.character(stage$stage), stage$imageID)[rownames(props)])
  expect_true(all(paste("median in group", levels(g)) %in% names(res)))
  ct <- res$cluster[1]
  expect_equal(res[1, paste("median in group", levels(g)[1])],
               signif(median(props[g == levels(g)[1], ct]), 2))
})

test_that("colTest() matches the rows of a matrix to a named condition", {
  data("diabetesData", package = "spicyR")
  props <- getProp(diabetesData)
  images <- unique(as.data.frame(SummarizedExperiment::colData(diabetesData))[, c("imageID", "stage")])
  condition <- setNames(as.character(images$stage), images$imageID)
  condition <- condition[condition %in% c("Long-duration", "Onset")]
  shuffled <- condition[sample(length(condition))]
  expect_identical(colTest(props, shuffled), colTest(props[names(condition), ], condition))
  # a named Surv outcome
  surv <- survival::Surv(seq_along(condition), rep(c(1, 0), length.out = length(condition)))
  names(surv) <- names(condition)
  expect_identical(colTest(props, surv[rev(seq_along(surv))], type = "survival"),
                   colTest(props[names(condition), ], surv, type = "survival"))
  # names that are not row names of props: taken in row order, as before
  expect_identical(colTest(props[names(condition), ], setNames(condition, paste0("x", seq_along(condition)))),
                   colTest(props[names(condition), ], unname(condition)))
})

test_that("colTest() takes a data.frame outcome whose row names are not image IDs in row order", {
  data("diabetesData", package = "spicyR")
  props <- getProp(diabetesData)
  outcome <- data.frame(time = seq_len(nrow(props)), event = rep(c(1, 0), length.out = nrow(props)))
  expect_identical(colTest(props, outcome, type = "survival"),
                   colTest(props, as.matrix(outcome), type = "survival"))
  rownames(outcome) <- rev(rownames(props))   # named rows are matched by name
  expect_identical(colTest(props, outcome, type = "survival"),
                   colTest(props, outcome[rownames(props), ], type = "survival"))
})
