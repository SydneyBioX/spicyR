## spicy() (cell method) from a SpatialData / AnnData Zarr store: the same results as the same cells in a data.frame.

skip_if_not_installed("Rarr")

store_cells <- function(seed = 2) {
  d <- sim_cells(seed, n_pat = 4, n_cells = 250)
  d$imageID <- factor(d$imageID, levels = rev(unique(d$imageID)))   # categories not in sorted order
  d$cellType <- factor(d$cellType)
  d$patient <- factor(d$patient)
  img <- unique(as.character(d$imageID)); d$age <- (seq_along(img) * 7 %% 13)[match(as.character(d$imageID), img)]
  d
}

expect_same_spicy <- function(d, path, ...) {
  a <- suppressMessages(spicy(d, ...))
  b <- suppressMessages(spicy(path, ...))
  expect_identical(b, a)
}

test_that("a store gives the results of the same cells as a data.frame", {
  d <- store_cells()
  path <- write_cell_store(d, tempfile(fileext = ".zarr"), chunk = 700L)
  expect_same_spicy(d, path, condition = "condition", subject = "patient", imageID = "imageID", r = 30)
  expect_same_spicy(d, path, condition = "condition", subject = "patient", imageID = "imageID", r = 30, effect = "count")
  expect_same_spicy(d, path, condition = "condition", imageID = "imageID", k = 10)
  expect_same_spicy(d, path, condition = "condition", subject = "patient", imageID = "imageID", r = c(10, 30, 50))
  expect_same_spicy(d, path, condition = "condition", imageID = "imageID", r = c(10, 30), effect = "count",
                    combine = "davies")
  expect_same_spicy(d, path, condition = "condition", imageID = "imageID", r = 30, covariates = "age",
                    adjustAbundance = TRUE, ref = "B")
  expect_same_spicy(d, path, condition = "condition", subject = "patient", imageID = "imageID", r = 30,
                    from = "T", to = c("tumour", "B"), cores = 3)
  # defaults: the region key and obsm "spatial"
  a <- suppressMessages(spicy(d, condition = "condition", subject = "patient", r = 30))
  expect_identical(suppressMessages(spicy(path, condition = "condition", subject = "patient", r = 30)), a)
  unlink(path, recursive = TRUE)
})

test_that("blocks of rows, batches of images, scattered images, strings and obs coordinates give the same results", {
  d <- store_cells(3)
  set.seed(4); d <- d[sample(nrow(d)), ]; rownames(d) <- NULL   # each image scattered over the rows
  d$cellType[c(5, 77, 1500)] <- NA                               # cells without a type
  d$imageID <- as.character(d$imageID); d$patient <- as.character(d$patient)
  d$condition <- factor(d$condition, levels = c("B", "A", "C"))  # an unused level
  path <- write_cell_store(d, tempfile(fileext = ".zarr"), coords_in_obs = TRUE, chunk = 333L)
  old <- options(spicyR.store_block = 1000, spicyR.store_batch = 1500); on.exit(options(old))
  for (cores in c(1, 2))
    expect_same_spicy(d, path, condition = "condition", subject = "patient", imageID = "imageID",
                      spatialCoords = c("x", "y"), r = c(20, 40), cores = cores)
  expect_same_spicy(d, path, condition = "condition", imageID = "imageID", spatialCoords = c("x", "y"), k = 8)
  # in row order, with string cell types too
  d2 <- d[order(d$imageID), ]; rownames(d2) <- NULL; d2$cellType <- as.character(d2$cellType)
  path2 <- write_cell_store(d2, tempfile(fileext = ".zarr"), chunk = 500L)
  expect_same_spicy(d2, path2, condition = "condition", subject = "patient", imageID = "imageID", r = 30,
                    effect = "count")
  unlink(c(path, path2), recursive = TRUE)
})

test_that("three conditions and a Zarr v2 store", {
  d <- store_cells(5)
  pat <- levels(d$patient); d$condition <- factor(c("x", "y", "z")[(match(d$patient, pat) - 1) %% 3 + 1])
  path <- write_cell_store(d, tempfile(fileext = ".zarr"), zarr_version = 2L)
  expect_same_spicy(d, path, condition = "condition", subject = "patient", imageID = "imageID", r = 30, ref = "z")
  unlink(path, recursive = TRUE)
})

test_that("a store is checked as a data.frame is", {
  d <- store_cells()
  d$stage <- d$condition; d$stage[c(1, 400)] <- rev(levels(d$stage))[as.integer(d$stage[c(1, 400)])]
  path <- write_cell_store(d, tempfile(fileext = ".zarr"))
  expect_error(spicy(path, condition = "stage", imageID = "imageID"), "must be constant within each image")
  expect_error(spicy(path, condition = "nope", imageID = "imageID"), "not found")
  expect_error(spicy(path, condition = "condition", subject = "nope", imageID = "imageID"), "`subject` column not found")
  expect_error(spicy(path, condition = "condition", covariates = "nope"), "covariates not found")
  expect_error(spicy(path, condition = "condition", table = "other"), "table 'other' not found")
  expect_error(spicy(path, condition = "condition", from = "nope"), "cell type not found")
  expect_error(spicy(path, condition = "condition", method = "image"), "method = \"cell\"")
  unlink(path, recursive = TRUE)
})

test_that("a store written by the Python spatialdata package is read", {
  path <- Sys.getenv("SPICYR_SDATA_STORE")
  skip_if(!nzchar(path) || !dir.exists(path), "no SpatialData store from Python given")
  tab <- file.path(path, "tables", "table")
  rd <- function(n) {
    p <- file.path(tab, "obs", n)
    if (identical(Rarr::read_zarr_attributes(p)[["encoding-type"]], "categorical")) {
      lev <- as.vector(Rarr::read_zarr_array(file.path(p, "categories")))
      factor(lev[as.vector(Rarr::read_zarr_array(file.path(p, "codes"))) + 1L], levels = lev)
    } else as.vector(Rarr::read_zarr_array(p))
  }
  key <- Rarr::read_zarr_attributes(tab)$region_key
  xy <- Rarr::read_zarr_array(file.path(tab, "obsm", "spatial"))
  d <- data.frame(imageID = rd(key), cellType = rd("cellType"), condition = rd("condition"), subject = rd("subject"),
                  x = xy[, 1], y = xy[, 2])
  a <- suppressMessages(spicy(d, condition = "condition", subject = "subject", r = 25))
  expect_identical(suppressMessages(spicy(path, condition = "condition", subject = "subject", r = 25)), a)
})
