## writeCellStore() of a whole SpatialExperiment: assays, rowData, colData, reducedDims, coordinates and images come
## back from spatialdataR and from Python's spatialdata.

skip_if_not_installed("Rarr")
skip_if_not_installed("SpatialExperiment")
skip_if_not_installed("Matrix")

full_spe <- function(n = 300, m = 7) {
  set.seed(1)
  counts <- Matrix::rsparsematrix(m, n, density = 0.3, rand.x = function(k) as.double(rpois(k, 5) + 1))
  dimnames(counts) <- list(paste0("G", 1:m), paste0("cell", 1:n))
  logc <- as.matrix(log1p(counts))
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts, logcounts = logc),
    rowData = S4Vectors::DataFrame(symbol = paste0("g", 1:m), mean = unname(rowMeans(logc))),
    colData = S4Vectors::DataFrame(imageID = rep(paste0("i", 1:6), each = n / 6), group = rep(c("A", "B"), each = n / 2),
                                   cellType = factor(sample(c("T", "B", "tumour"), n, TRUE)),
                                   score = ifelse(seq_len(n) %% 11 == 0, NA, seq_len(n))),
    sample_id = rep(c("s1", "s2"), each = n / 2),
    spatialCoords = cbind(x = runif(n, 0, 500), y = runif(n, 0, 600)))
  SingleCellExperiment::reducedDim(spe, "PCA") <- matrix(rnorm(n * 3), n, 3)
  png <- list.files(system.file("extdata", "10xVisium", package = "SpatialExperiment"), "lowres", recursive = TRUE,
                    full.names = TRUE)
  for (k in 1:2)
    spe <- SpatialExperiment::addImg(spe, imageSource = png[k], scaleFactor = 0.05, sample_id = paste0("s", k),
                                     image_id = "lowres", load = TRUE)
  spe
}

rgb_of <- function(spe, s) {
  ras <- as.matrix(SpatialExperiment::imgRaster(SpatialExperiment::getImg(spe, s, "lowres")))
  array(t(grDevices::col2rgb(as.vector(ras))), c(dim(ras), 3L))   # y, x, channel
}

test_that("a SpatialExperiment is written whole, and spatialdataR reads it back", {
  skip_if_not_installed("spatialdataR")
  old <- options(spicyR.store_write_chunk = 64L); on.exit(options(old))   # several chunks
  spe <- full_spe()
  path <- writeCellStore(spe, tempfile(fileext = ".zarr"))
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  sd <- suppressWarnings(spatialdataR::readSpatialData(path))
  tab <- spatialdataR::table(sd)
  expect_equal(unname(as.matrix(SummarizedExperiment::assay(tab, "X"))),
               unname(as.matrix(SummarizedExperiment::assay(spe, "counts"))))
  expect_equal(unname(as.matrix(SummarizedExperiment::assay(tab, "logcounts"))),
               unname(SummarizedExperiment::assay(spe, "logcounts")))
  expect_identical(rownames(tab), rownames(spe))
  expect_identical(colnames(tab), colnames(spe))
  expect_identical(as.character(SummarizedExperiment::rowData(tab)$symbol), SummarizedExperiment::rowData(spe)$symbol)
  expect_equal(unname(SummarizedExperiment::rowData(tab)$mean), SummarizedExperiment::rowData(spe)$mean)
  expect_identical(as.character(tab$cellType), as.character(spe$cellType))
  expect_identical(as.numeric(tab$score), as.numeric(spe$score))
  expect_equal(unname(as.matrix(SingleCellExperiment::reducedDim(tab, "PCA"))),
               unname(SingleCellExperiment::reducedDim(spe, "PCA")))
  expect_equal(unname(as.matrix(SingleCellExperiment::reducedDim(tab, "spatial"))),
               unname(SpatialExperiment::spatialCoords(spe)))
  for (s in c("s1", "s2")) {
    a <- as.array(spatialdataR::data(spatialdataR::image(sd, paste0(s, "_lowres"))))   # channel, y, x
    expect_identical(as.vector(aperm(a, c(2, 3, 1))), as.vector(rgb_of(spe, s)))
  }
  # the cells of each image as points, in the coordinate system of their sample
  xy <- SpatialExperiment::spatialCoords(spe)
  for (im in unique(spe$imageID)) {
    pt <- as.data.frame(spatialdataR::point(sd, im))
    pt <- pt[order(pt$instance_id), ]
    rows <- which(spe$imageID == im)
    expect_identical(as.integer(pt$instance_id), rows - 1L)
    expect_equal(unname(sf::st_coordinates(pt$geometry)), unname(xy[rows, ]))
  }
  # spicy() reads the store as it reads the object
  a <- suppressMessages(spicy(spe, condition = "group", cellType = "cellType", r = 40))
  b <- suppressMessages(spicy(path, condition = "group", cellType = "cellType", r = 40))
  expect_identical(b, a)
})

test_that("Python's spatialdata reads a whole SpatialExperiment", {
  py <- "/dskh/nobackup/ellisp/team-patrick/scratch/venv_spicyr_dev/bin/python"
  skip_if_not(file.exists(py), "the spicyr-py development Python is not available")
  old <- options(spicyR.store_write_chunk = 64L); on.exit(options(old))
  spe <- full_spe()
  path <- writeCellStore(spe, tempfile(fileext = ".zarr"))
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  code <- paste0(
    "import json, warnings, numpy as np, anndata, spatialdata as sd\n",
    "from spatialdata.transformations import get_transformation\n",
    "warnings.simplefilter('error', UserWarning)\n",   # e.g. a table annotating an element that is not there

    "warnings.simplefilter('error', anndata._warnings.OldFormatWarning)\n",
    "s = sd.read_zarr('", path, "'); t = s.tables['table']\n",
    "out = dict(X=t.X.toarray().tolist(), logcounts=t.layers['logcounts'].tolist(),\n",
    "  spatial=t.obsm['spatial'].tolist(), PCA=t.obsm['PCA'].tolist(), obs=list(t.obs_names), var=list(t.var_names),\n",
    "  symbol=list(t.var['symbol']), na=int(t.obs['score'].isna().sum()))\n",
    "for k in ['s1', 's2']:\n",
    "    im = s.images[k + '_lowres']\n",
    "    out[k] = np.moveaxis(im.values, 0, -1).ravel(order='F').tolist()\n",
    "    out[k + '_scale'] = list(get_transformation(im, to_coordinate_system=k).scale)\n",
    "for k in sorted(s.points):\n",
    "    p = s.points[k].compute().sort_values('instance_id')\n",
    "    out['pt_' + k] = dict(x=p['x'].tolist(), id=p['instance_id'].tolist(),\n",
    "                         cs=list(get_transformation(s.points[k], get_all=True)))\n",
    "print(json.dumps(out))\n")
  f <- tempfile(fileext = ".py"); writeLines(code, f)
  res <- system2(py, c("-I", f), stdout = TRUE, stderr = FALSE)
  expect_null(attr(res, "status"))
  got <- jsonlite::fromJSON(res[length(res)])
  expect_equal(got$X, unname(t(as.matrix(SummarizedExperiment::assay(spe, "counts")))))
  expect_equal(got$logcounts, unname(t(SummarizedExperiment::assay(spe, "logcounts"))))
  expect_equal(got$spatial, unname(SpatialExperiment::spatialCoords(spe)))
  expect_equal(got$PCA, unname(SingleCellExperiment::reducedDim(spe, "PCA")))
  expect_identical(got$obs, colnames(spe))
  expect_identical(got$var, rownames(spe))
  expect_identical(got$symbol, SummarizedExperiment::rowData(spe)$symbol)
  expect_identical(got$na, sum(is.na(spe$score)))
  for (im in unique(spe$imageID)) {
    rows <- which(spe$imageID == im)
    expect_equal(got[[paste0("pt_", im)]]$x, unname(SpatialExperiment::spatialCoords(spe)[rows, 1]))
    expect_identical(as.integer(got[[paste0("pt_", im)]]$id), rows - 1L)
    expect_identical(got[[paste0("pt_", im)]]$cs, spe$sample_id[rows[1]])
  }
  for (s in c("s1", "s2")) {
    expect_identical(as.integer(got[[s]]), as.integer(as.vector(rgb_of(spe, s))))
    expect_equal(got[[paste0(s, "_scale")]], c(1, 20, 20))
  }
})

test_that("unsorted cells over many chunks give one points file per image", {
  skip_if_not_installed("spatialdataR")
  old <- options(spicyR.store_write_chunk = 50L); on.exit(options(old))
  set.seed(3)
  d <- data.frame(imageID = sample(c("a", "b", "c"), 1000, TRUE), x = runif(1000), y = runif(1000))
  path <- writeCellStore(d, tempfile(fileext = ".zarr"))
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  for (im in c("a", "b", "c")) {
    expect_identical(list.files(file.path(path, "points", im, "points.parquet")), "part.0.parquet")
    pt <- as.data.frame(spatialdataR::readPoint(file.path(path, "points", im)))
    pt <- pt[order(pt$instance_id), ]
    expect_identical(as.integer(pt$instance_id), which(d$imageID == im) - 1L)
    expect_equal(unname(sf::st_coordinates(pt$geometry)[, 2]), d$y[d$imageID == im])
  }
  expect_false(dir.exists(file.path(path, "points", ".tmp")))
  # image IDs that are not element names: no points, with a warning
  d$imageID <- sub("a", "a b", d$imageID)
  expect_warning(p2 <- writeCellStore(d, tempfile(fileext = ".zarr")), "valid SpatialData element names")
  expect_false(dir.exists(file.path(p2, "points")))
  unlink(p2, recursive = TRUE)
})
