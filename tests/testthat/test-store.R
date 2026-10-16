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
  expect_error(spicy(path, condition = "condition", tableKey = "other"), "table 'other' not found")
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
  ct <- intersect(c("cellType", "cell_type"), unlist(Rarr::read_zarr_attributes(file.path(tab, "obs"))[["column-order"]]))[1]
  xy <- Rarr::read_zarr_array(file.path(tab, "obsm", "spatial"))
  d <- data.frame(imageID = rd(key), cellType = rd(ct), condition = rd("condition"), subject = rd("subject"),
                  x = xy[, 1], y = xy[, 2])
  rm(xy)
  a <- suppressMessages(spicy(d, condition = "condition", subject = "subject", r = 25))
  expect_identical(suppressMessages(spicy(path, condition = "condition", subject = "subject", cellType = ct, r = 25)), a)
})

test_that("a SpatialData object (spatialdataR) is read from its store, or from its table in memory", {
  pkg <- "spatialdataR"   # not on CRAN or Bioconductor: used only when installed
  skip_if_not_installed(pkg)
  skip_if_not_installed("SingleCellExperiment")
  path <- Sys.getenv("SPICYR_SDATA_STORE")
  skip_if(!nzchar(path) || !dir.exists(path), "no SpatialData store from Python given")
  x <- suppressWarnings(getExportedValue(pkg, "readSpatialData")(path, images = FALSE, labels = FALSE, points = FALSE,
                                                                 shapes = FALSE))
  ct <- intersect(c("cellType", "cell_type"), unlist(Rarr::read_zarr_attributes(file.path(path, "tables", "table", "obs"))[["column-order"]]))[1]
  a <- suppressMessages(spicy(path, condition = "condition", subject = "subject", cellType = ct, r = 25))
  expect_identical(suppressMessages(spicy(x, condition = "condition", subject = "subject", cellType = ct, r = 25)), a)
  # without a store behind it: the table in memory
  tab <- methods::slot(x, "tables")
  sce <- tab[["table"]]; m <- SingleCellExperiment::int_metadata(sce); m$source_path <- NULL
  SingleCellExperiment::int_metadata(sce) <- m; tab[["table"]] <- sce; methods::slot(x, "tables") <- tab
  b <- suppressMessages(spicy(x, condition = "condition", subject = "subject", cellType = ct, r = 25))
  expect_identical(b, a)
})

test_that("imageInfo gives the image-level columns that the cell table lacks", {
  d <- store_cells(6)
  first <- !duplicated(d$imageID)
  info <- d[first, c("imageID", "condition", "patient", "age")]
  rownames(info) <- as.character(info$imageID); info <- info[rev(seq_len(nrow(info))), ]   # any row order
  lean <- d[, setdiff(names(d), c("condition", "patient", "age"))]
  a <- suppressMessages(spicy(d, condition = "condition", subject = "patient", imageID = "imageID", r = 30,
                              covariates = "age"))
  # a data.frame, by its imageID column or by its row names
  expect_identical(suppressMessages(spicy(lean, condition = "condition", subject = "patient", imageID = "imageID",
                                          r = 30, covariates = "age", imageInfo = info)), a)
  expect_identical(suppressMessages(spicy(lean, condition = "condition", subject = "patient", imageID = "imageID",
                                          r = 30, covariates = "age", imageInfo = info[, -1])), a)
  # columns of the cell table come from the cells; imageInfo only adds the others
  expect_identical(suppressMessages(spicy(d, condition = "condition", subject = "patient", imageID = "imageID",
                                          r = 30, covariates = "age", imageInfo = info["age"])), a)
  # from a store: a data.frame, or a table of the store
  path <- write_cell_store(lean, tempfile(fileext = ".zarr"), chunk = 700L)
  write_info_table(info[, -1], path)
  for (ii in list(info, "images"))
    expect_identical(suppressMessages(spicy(path, condition = "condition", subject = "patient", imageID = "imageID",
                                            r = 30, covariates = "age", imageInfo = ii)), a)
  # unknown and missing images
  extra <- rbind(info, info[1, ]); extra$imageID <- as.character(extra$imageID); extra$imageID[nrow(extra)] <- "nope"
  expect_error(spicy(lean, condition = "condition", imageID = "imageID", imageInfo = extra), "not found among the cells: nope")
  expect_error(spicy(lean, condition = "condition", imageID = "imageID", imageInfo = info[-1, ]), "images missing from `imageInfo`")
  expect_error(spicy(path, condition = "condition", imageID = "imageID", imageInfo = info[-1, ]), "images missing from `imageInfo`")
  expect_error(spicy(lean, condition = "condition", imageID = "imageID", imageInfo = rbind(info, info[1, ])),
               "more than one row")
  unlink(path, recursive = TRUE)
})

test_that("a survival outcome from imageInfo, for a data.frame and for a store", {
  d <- sim_cells(7); d$imageID <- factor(d$imageID); d$cellType <- factor(d$cellType); d$patient <- factor(d$patient)
  pat <- levels(d$patient); set.seed(8)
  tm <- stats::rexp(length(pat), 0.1); ev <- stats::rbinom(length(pat), 1, 0.7)
  first <- !duplicated(d$imageID)
  info <- d[first, c("imageID", "patient")]; rownames(info) <- as.character(info$imageID)
  i <- match(info$patient, pat); info$os <- survival::Surv(tm[i], ev[i])
  d$os <- info$os[match(d$imageID, info$imageID)]
  a <- suppressMessages(spicy(d, condition = "os", subject = "patient", imageID = "imageID", r = 30, from = "tumour",
                              to = c("T", "B")))
  expect_true(all(c("hazard_ratio_sd", "p_value") %in% names(a$cellResults)))
  lean <- d[, setdiff(names(d), c("os", "condition"))]
  path <- write_cell_store(lean, tempfile(fileext = ".zarr"), chunk = 700L)
  path2 <- write_cell_store(lean[, setdiff(names(lean), "patient")], tempfile(fileext = ".zarr"), chunk = 700L)
  expect_identical(suppressMessages(spicy(lean, condition = "os", subject = "patient", imageID = "imageID", r = 30,
                                          from = "tumour", to = c("T", "B"), imageInfo = info["os"])), a)
  expect_identical(suppressMessages(spicy(path, condition = "os", subject = "patient", imageID = "imageID", r = 30,
                                          from = "tumour", to = c("T", "B"), imageInfo = info["os"])), a)
  # the subject from imageInfo too
  expect_identical(suppressMessages(spicy(path2, condition = "os", subject = "patient", imageID = "imageID", r = 30,
                                          from = "tumour", to = c("T", "B"), imageInfo = info)), a)
  unlink(c(path, path2), recursive = TRUE)
})

test_that("writeCellStore() writes a store that spicy() and spatialdataR read", {
  old <- options(spicyR.store_write_chunk = 700L); on.exit(options(old))   # several chunks, the last one partial
  d <- store_cells(9)
  info <- unique(d[, c("imageID", "condition", "patient")]); rownames(info) <- as.character(info$imageID)
  info$imageID <- NULL
  path <- writeCellStore(d[, c("imageID", "cellType", "x", "y", "age")], tempfile(fileext = ".zarr"), imageInfo = info)
  a <- suppressMessages(spicy(d, condition = "condition", subject = "patient", r = 30, covariates = "age"))
  expect_identical(suppressMessages(spicy(path, condition = "condition", subject = "patient", r = 30,
                                          covariates = "age", imageInfo = "images")), a)
  skip_if_not_installed("spatialdataR")
  x <- suppressWarnings(spatialdataR::readSpatialData(path))
  expect_identical(suppressMessages(spicy(x, condition = "condition", subject = "patient", r = 30, covariates = "age",
                                          imageInfo = "images")), a)
  unlink(path, recursive = TRUE)
})

test_that("a store with character image IDs and missing values is read by Python's spatialdata", {
  py <- "/dskh/nobackup/ellisp/team-patrick/scratch/venv_spicyr_dev/bin/python"
  skip_if_not(file.exists(py), "the spicyr-py development Python is not available")
  d <- store_cells(4); d$imageID <- as.character(d$imageID)
  d$score <- ifelse(seq_len(nrow(d)) %% 7 == 0, NA, seq_len(nrow(d)))   # a nullable integer column
  path <- writeCellStore(d[, c("imageID", "cellType", "x", "y", "score")], tempfile(fileext = ".zarr"))
  # an error if spatialdata refuses the table, or anndata warns of an old format
  code <- paste0("import warnings, anndata, spatialdata as sd; ",
                 "warnings.simplefilter('error', anndata._warnings.OldFormatWarning); ",
                 "t = sd.read_zarr('", path, "').tables['table']; print(t.obs['imageID'].dtype, t.n_obs)")
  out <- suppressWarnings(system2(py, c("-c", shQuote(code)), stdout = TRUE, stderr = TRUE))
  expect_null(attr(out, "status"))
  expect_true(any(grepl(paste("category", nrow(d)), out)))
  unlink(path, recursive = TRUE)
})
