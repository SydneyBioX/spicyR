## Regression tests for the reference/base condition group bug: modelDataGen() used to
## rebuild `condition` from the raw column per pair, and computeImage()'s per-image
## as.factor() on a scalar collapsed to a single-level factor -- so dplyr::bind_rows()'s
## level-union-by-first-appearance order could differ pair to pair, depending on which
## images survived that specific pair's cell-type filtering. Fixed by factoring
## `df$condition` once in modelDataGen(), on the full dataset, before any per-image split.

skip_if_not_installed("SpatialExperiment")
skip_if_not_installed("brglm2")

suppressWarnings(suppressPackageStartupMessages(library(SpatialExperiment)))

## --- synthetic data: engineered so different pairs see different "first surviving
## image", which is what exposes the bug. condition is plain character, not pre-factored. --

make_spe_asymmetric <- function(seed = 1, cells_per_img = 300) {
  set.seed(seed)
  imgs <- data.frame(
    imageID   = c("Img1", "Img2", "Img3", "Img4", "Img5", "Img6"),
    condition = c("A", "B", "A", "B", "A", "B"),
    hasRare   = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  blocks <- list()
  for (i in seq_len(nrow(imgs))) {
    types <- c("Tcell", "Tumour")
    probs <- c(0.5, 0.5)
    if (imgs$hasRare[i]) {
      types <- c(types, "Rare")
      probs <- c(0.45, 0.45, 0.10)
    }
    ct <- sample(types, cells_per_img, replace = TRUE, prob = probs)
    x <- runif(cells_per_img, 0, 1000); y <- runif(cells_per_img, 0, 1000)
    blocks[[i]] <- data.frame(
      imageID = imgs$imageID[i], subject = imgs$imageID[i],
      condition = imgs$condition[i], cellType = ct, x = x, y = y,
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, blocks)
  SpatialExperiment(
    assays = list(counts = matrix(0, 1, nrow(df))),
    colData = DataFrame(df[, c("imageID", "subject", "condition", "cellType")]),
    spatialCoords = as.matrix(df[, c("x", "y")])
  )
}

## Pair "Tcell__Tumour": present in Img1..Img6 -> first surviving image is Img1 (condition A).
## Pair "Rare__Tumour": Img1 has no Rare cells -> first surviving image is Img2 (condition B).
## Both conditions keep >= 2 qualifying patients for both pairs (CR2 needs at least two
## independent patients per group), so neither pair is skipped -- they just used to disagree
## on which group is the reference.

extract_base_group <- function(msgs) {
  hit <- grep("as base comparison group", msgs, value = TRUE)
  stopifnot(length(hit) >= 1)
  sub(".*Using condition = (\\S+) as base.*", "\\1", hit[1])
}

## --- symptom 1 + symptom 2: cross-pair consistency and agreement with the announcement ---

test_that("poisson: every pair agrees on conditionRef, matching the announced base group", {
  spe <- make_spe_asymmetric()
  msgs <- character(0)
  res <- withCallingHandlers(
    suppressWarnings(spicyGLM(spe, condition = "condition", subject = "subject",
                              imageID = "imageID", cellType = "cellType",
                              spatialCoords = c("x", "y"), cores = 1,
                              r = 50, from = c("Tcell", "Rare"), to = "Tumour",
                              family = "poisson")),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  )

  gr <- res$GLMresults
  target <- gr[paste0(gr$from, "__", gr$to) %in% c("Tcell__Tumour", "Rare__Tumour"), ]
  expect_equal(nrow(target), 2)
  expect_length(unique(target$conditionRef), 1)
  expect_equal(unique(target$conditionRef), extract_base_group(msgs))
})

test_that("binomial: every pair agrees on conditionRef, matching the announced base group", {
  spe <- make_spe_asymmetric()
  msgs <- character(0)
  res <- withCallingHandlers(
    suppressWarnings(spicyGLM(spe, condition = "condition", subject = "subject",
                              imageID = "imageID", cellType = "cellType",
                              spatialCoords = c("x", "y"), cores = 1,
                              k = 10, from = c("Tcell", "Rare"), to = "Tumour",
                              family = "binomial")),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  )

  gr <- res$GLMresults
  target <- gr[paste0(gr$from, "__", gr$to) %in% c("Tcell__Tumour", "Rare__Tumour"), ]
  expect_equal(nrow(target), 2)
  expect_length(unique(target$conditionRef), 1)
  expect_equal(unique(target$conditionRef), extract_base_group(msgs))
})

## --- explicit ref parameter -------------------------------------------------

test_that("ref lets the user override the base group, consistently across pairs", {
  spe <- make_spe_asymmetric()
  res <- suppressMessages(suppressWarnings(
    spicyGLM(spe, condition = "condition", subject = "subject",
             imageID = "imageID", cellType = "cellType",
             spatialCoords = c("x", "y"), cores = 1,
             r = 50, from = c("Tcell", "Rare"), to = "Tumour",
             family = "poisson", ref = "B")
  ))

  gr <- res$GLMresults
  target <- gr[paste0(gr$from, "__", gr$to) %in% c("Tcell__Tumour", "Rare__Tumour"), ]
  expect_equal(nrow(target), 2)
  expect_equal(unique(target$conditionRef), "B")
})

test_that("an invalid ref errors with the available levels", {
  spe <- make_spe_asymmetric()
  expect_error(
    suppressMessages(suppressWarnings(
      spicyGLM(spe, condition = "condition", subject = "subject",
               imageID = "imageID", cellType = "cellType",
               spatialCoords = c("x", "y"), cores = 1,
               r = 50, from = "Tcell", to = "Tumour",
               family = "poisson", ref = "NotALevel")
    )),
    "not a level of `condition`.*A.*B"
  )
})
