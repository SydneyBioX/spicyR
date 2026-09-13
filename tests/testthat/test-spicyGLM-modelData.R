## Tests for `storeModelData` on spicyGLM(): opt-in attachment of the raw
## per-image, per-reference-cell table (modelDataGen()/getPairwiseAssoc()
## output) that the GLM fit was built from, for both poisson and binomial.

skip_if_not_installed("SpatialExperiment")
skip_if_not_installed("brglm2")

suppressWarnings(suppressPackageStartupMessages(library(SpatialExperiment)))

## --- synthetic data (mirrors test-spicyGLM-binomial.R's make_spe()) --------

make_spe <- function(seed = 42, n_per_cond = 6, cells_per_img = 400,
                     enrich_cond = "B") {
  set.seed(seed)
  conds <- c("A", "B")
  blocks <- list()
  for (cond in conds) {
    for (s in seq_len(n_per_cond)) {
      sid <- paste0(cond, "_s", s)
      nc <- cells_per_img
      ct <- sample(c("Tcell", "Tumour", "Bcell"), nc, replace = TRUE,
                   prob = c(0.3, 0.4, 0.3))
      x <- runif(nc, 0, 1000); y <- runif(nc, 0, 1000)
      if (cond == enrich_cond) {
        tum <- ct == "Tumour"
        cx <- mean(x[tum]); cy <- mean(y[tum])
        tt <- which(ct == "Tcell")
        mv <- sample(tt, length(tt) %/% 2)
        x[mv] <- 0.5 * x[mv] + 0.5 * cx
        y[mv] <- 0.5 * y[mv] + 0.5 * cy
      }
      blocks[[length(blocks) + 1]] <- data.frame(
        imageID = sid, subject = sid, condition = cond,
        cellType = ct, x = x, y = y, stringsAsFactors = FALSE
      )
    }
  }
  df <- do.call(rbind, blocks)
  SpatialExperiment(
    assays = list(counts = matrix(0, 1, nrow(df))),
    colData = DataFrame(df[, c("imageID", "subject", "condition", "cellType")]),
    spatialCoords = as.matrix(df[, c("x", "y")])
  )
}

run_glm <- function(spe, ..., quiet = TRUE) {
  f <- function() spicyGLM(spe, condition = "condition", subject = "subject",
                           imageID = "imageID", cellType = "cellType",
                           spatialCoords = c("x", "y"), cores = 1, ...)
  if (quiet) suppressMessages(suppressWarnings(f())) else f()
}

poisson_cols <- c("from", "to", "condition", "subject", "imageID", "cellID", "n", "density")
binomial_cols <- c("from", "to", "condition", "subject", "imageID", "cellID", "n", "k", "prop", "p0")

## --- default: no modelData attached ----------------------------------------

test_that("storeModelData defaults to FALSE (no $modelData attached)", {
  spe <- make_spe()
  poisRes <- run_glm(spe, r = 50, from = "Tcell", to = "Tumour", family = "poisson")
  binRes  <- run_glm(spe, k = 10, from = "Tcell", to = "Tumour", family = "binomial")
  expect_null(poisRes$modelData)
  expect_null(binRes$modelData)
})

## --- single pair -------------------------------------------------------

test_that("storeModelData attaches the correct table for a single poisson pair", {
  spe <- make_spe()
  res <- run_glm(spe, r = 50, from = "Tcell", to = "Tumour", family = "poisson",
                 storeModelData = TRUE)
  md <- res$modelData
  expect_s3_class(md, "data.frame")
  expect_identical(names(md), poisson_cols)
  expect_true(all(md$from == "Tcell"))
  expect_true(all(md$to == "Tumour"))
  expect_true(all(md$n >= 0))
  expect_true(all(md$density > 0))
})

test_that("storeModelData attaches the correct table for a single binomial pair, with prop = n/k", {
  spe <- make_spe()
  res <- run_glm(spe, k = 10, from = "Tcell", to = "Tumour", family = "binomial",
                 storeModelData = TRUE)
  md <- res$modelData
  expect_s3_class(md, "data.frame")
  expect_identical(names(md), binomial_cols)
  expect_true(all(md$k == 10))
  expect_equal(md$prop, md$n / md$k)
})

## --- multi-pair ----------------------------------------------------------

test_that("storeModelData attaches a named list keyed 'from__to' for multi-pair runs", {
  spe <- make_spe()
  res <- run_glm(spe, k = 12, family = "binomial", storeModelData = TRUE)
  expect_type(res$modelData, "list")
  expect_true(all(grepl("__", names(res$modelData))))
  expect_true(all(vapply(names(res$modelData), function(nm) {
    df <- res$modelData[[nm]]
    if (nrow(df) == 0) return(TRUE)
    identical(names(df), binomial_cols)
  }, logical(1))))
})

test_that("storeModelData poisson multi-pair table matches getPairwiseAssoc() shape", {
  spe <- make_spe()
  res <- run_glm(spe, r = 50, family = "poisson", storeModelData = TRUE)
  expect_type(res$modelData, "list")
  nonEmpty <- Filter(function(df) nrow(df) > 0, res$modelData)
  expect_true(length(nonEmpty) > 0)
  for (df in nonEmpty) expect_identical(names(df), poisson_cols)
})
