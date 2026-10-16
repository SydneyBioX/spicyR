## ---- the cell method from a SpatialData (Zarr) store ------------------------------------------------------------
## A SpatialData store's cell table is an AnnData group (tables/<name>) whose columns are Zarr arrays, read here with
## Rarr. Nothing of the length of the table is held in R: one pass over the image, cell type, condition and subject
## columns, a block of rows at a time, gives the cells per image and type, the first row of each image (and of each
## image and type) and the image-level values; the coordinates and cell types are then read a batch of whole images
## at a time for the neighbour sums (cell_batch_totals). The context has the fields .cell_context() gives a
## data.frame of the same cells, so every later step is shared.

## The obs of a (small) table of the store as a data.frame, its index as the row names: a per-image table.
.store_obs_frame <- function(path, table) {
  # a table's own group (as a SpatialData object read from a store points to): the store is two levels up
  if (basename(dirname(path)) == "tables" &&
      identical(tryCatch(Rarr::read_zarr_attributes(path)[["encoding-type"]], error = function(e) NULL), "anndata"))
    path <- dirname(dirname(path))
  st <- .store_table(path, table)
  obs <- Rarr::read_zarr_attributes(file.path(st$path, "obs"))
  n <- st$n
  d <- lapply(stats::setNames(st$obs_columns, st$obs_columns), function(nm) .store_read(.store_column(st, nm), seq_len(n)))
  d <- as.data.frame(d, stringsAsFactors = FALSE, optional = TRUE)
  if (!is.null(obs[["_index"]])) rownames(d) <- as.character(as.vector(Rarr::read_zarr_array(file.path(st$path, "obs", obs[["_index"]]))))
  d
}

## Whether `x` names a Zarr store (a directory).
.is_store <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && dir.exists(x)

## Rows read per block in the pass over the columns, cells per thread in a batch of whole images, and at most this
## many cells in a batch of more than one image.
.store_block <- function() getOption("spicyR.store_block", 2^23)
.store_batch <- function() getOption("spicyR.store_batch", 2^24)
.store_batch_max <- function() getOption("spicyR.store_batch_max", 2^26)

## The AnnData group of the cell table: the store itself when it is an AnnData group, else tables/<table> (the only
## table, or "table", when `table` is NULL).
.store_table <- function(path, table = NULL) {
  .need("Rarr", "to read a SpatialData store")
  attrs <- function(p) tryCatch(Rarr::read_zarr_attributes(p), error = function(e) list())
  if (identical(attrs(path)[["encoding-type"]], "anndata") && is.null(table)) {
    tab <- path
  } else {
    tdir <- file.path(path, "tables")
    have <- if (dir.exists(tdir)) basename(list.dirs(tdir, recursive = FALSE)) else character(0)
    have <- have[vapply(file.path(tdir, have), function(p) identical(attrs(p)[["encoding-type"]], "anndata"), TRUE)]
    if (!length(have)) stop("no table found in the SpatialData store ", path, ".", call. = FALSE)
    if (is.null(table)) table <- if (length(have) == 1L) have else if ("table" %in% have) "table" else
      stop("the store has several tables (", paste(have, collapse = ", "), "); choose one with `table`.", call. = FALSE)
    if (!table %in% have) stop("table '", table, "' not found; the store has: ", paste(have, collapse = ", "), call. = FALSE)
    tab <- file.path(tdir, table)
  }
  obs <- attrs(file.path(tab, "obs"))
  idx <- obs[["_index"]]
  n <- if (!is.null(idx)) .store_dim(file.path(tab, "obs", idx))[1] else NA
  list(path = tab, attrs = attrs(tab), obs_columns = unlist(obs[["column-order"]]), n = as.numeric(n),
       obsm = if (dir.exists(file.path(tab, "obsm"))) basename(list.dirs(file.path(tab, "obsm"), recursive = FALSE)) else character(0))
}

.store_dim <- function(p) unlist(Rarr::zarr_overview(p, as_data_frame = TRUE)$dim)
.store_chunk <- function(p) unlist(Rarr::zarr_overview(p, as_data_frame = TRUE)$chunk_dim)

## How to read one obs column: a categorical (integer codes and their categories), a nullable array (values and a
## mask) or a plain array (numbers or strings).
.store_column <- function(st, name) {
  if (!name %in% st$obs_columns && !dir.exists(file.path(st$path, "obs", name)))
    stop("column '", name, "' not found in the table's obs; columns: ", paste(st$obs_columns, collapse = ", "), call. = FALSE)
  p <- file.path(st$path, "obs", name)
  a <- Rarr::read_zarr_attributes(p)
  type <- a[["encoding-type"]]
  if (identical(type, "categorical")) {
    cats <- as.vector(Rarr::read_zarr_array(file.path(p, "categories")))
    return(list(name = name, kind = "categorical", path = file.path(p, "codes"), levels = as.character(cats),
                ordered = isTRUE(a$ordered), n = .store_dim(file.path(p, "codes"))[1]))
  }
  if (!is.null(type) && startsWith(type, "nullable"))
    return(list(name = name, kind = "nullable", path = file.path(p, "values"), mask = file.path(p, "mask"),
                n = .store_dim(file.path(p, "values"))[1]))
  list(name = name, kind = "array", path = p, n = .store_dim(p)[1])
}

## Rows `rows` (increasing) of a column as an R vector: a categorical as a factor, a nullable array with NA where
## masked, strings as character.
.store_read <- function(col, rows) {
  v <- as.vector(Rarr::read_zarr_array(col$path, index = list(rows)))
  if (col$kind == "categorical") {
    v <- as.integer(v) + 1L; v[v < 1L] <- NA_integer_
    return(structure(v, levels = col$levels, class = c(if (col$ordered) "ordered", "factor")))
  }
  if (col$kind == "nullable") v[as.logical(Rarr::read_zarr_array(col$mask, index = list(rows)))] <- NA
  v
}

## Codes of rows a..b of a column: for a categorical, its own codes (0-based, -1 where missing; base 0); for any
## other column, 1-based codes against `lev`, a growing record of the values in order of first appearance (an
## environment; base 1).
.store_codes <- function(col, a, b, lev = NULL) {
  if (col$kind == "categorical") {
    v <- Rarr::read_zarr_array(col$path, index = list(a:b)); dim(v) <- NULL
    if (!is.integer(v)) v <- as.integer(v)
    return(structure(v, base = 0L))
  }
  v <- .store_read(col, a:b)
  u <- unique(v); u <- u[!is.na(u)]
  lev$values <- c(lev$values, u[!u %in% lev$values])
  structure(match(v, lev$values), base = 1L)
}

## Rows a..b of an image-level column as compared within images: a categorical's codes, else its values.
.store_values <- function(col, a, b) {
  if (col$kind != "categorical") return(.store_read(col, a:b))
  v <- Rarr::read_zarr_array(col$path, index = list(a:b)); dim(v) <- NULL
  if (!is.integer(v)) v <- as.integer(v)
  v
}

## Values from .store_values() as R values (a categorical's codes as a factor).
.store_as_values <- function(col, v) {
  if (col$kind != "categorical") return(v)
  v <- as.integer(v) + 1L; v[v < 1L] <- NA_integer_
  structure(v, levels = col$levels, class = c(if (col$ordered) "ordered", "factor"))
}

## Coordinates of rows a..b: an obsm array (first two columns) or two obs columns.
.store_coords <- function(st, coords, a, b) {
  if (coords$obsm) {
    if (coords$both) {   # both columns in each chunk: one read
      m <- Rarr::read_zarr_array(coords$path, index = list(a:b, 1:2))
      return(list(x = m[, 1], y = m[, 2]))
    }
    col <- function(j) { v <- Rarr::read_zarr_array(coords$path, index = list(a:b, j)); dim(v) <- NULL; v }
    return(list(x = col(1L), y = col(2L)))
  }
  list(x = .store_read(coords$x, a:b), y = .store_read(coords$y, a:b))
}

## The coordinates: `spatialCoords` names two obs columns or one obsm array; by default the obsm array "spatial".
.store_coord_spec <- function(st, spatialCoords) {
  if (length(spatialCoords) == 1L) {
    if (!spatialCoords %in% st$obsm) stop("obsm '", spatialCoords, "' not found in the table; obsm: ",
                                          paste(st$obsm, collapse = ", "), call. = FALSE)
    p <- file.path(st$path, "obsm", spatialCoords)
    if (length(.store_dim(p)) != 2L || .store_dim(p)[2] < 2L) stop("obsm '", spatialCoords, "' needs two columns.", call. = FALSE)
    return(list(obsm = TRUE, path = p, label = spatialCoords, both = .store_chunk(p)[2] >= 2L))
  }
  list(obsm = FALSE, x = .store_column(st, spatialCoords[1]), y = .store_column(st, spatialCoords[2]),
       label = spatialCoords)
}

## Whether v (rows of the table) equals, as as.character() sees it, ref at each row's image (codes im); NA equals NA.
.same_value <- function(v, im, ref) {
  d <- cell_rows_differing(v, im, ref, attr(im, "base"))
  !length(d) || identical(as.character(v[d]), as.character(ref[im[d] + 1L - attr(im, "base")]))
}

## The context of the cell method (as .cell_context(by_image = TRUE) gives it for the same cells as a data.frame)
## from one pass over the store's columns.
.cell_context_store <- function(st, condition, subject, covariates, image_id, cell_type, coords, ref = NULL) {
  N <- st$n
  if (!is.finite(N) || N < 1) stop("the table has no cells.", call. = FALSE)
  if (N > .Machine$integer.max) stop("tables of more than 2^31 - 1 cells are not supported.", call. = FALSE)
  N <- as.integer(N)
  ci <- .store_column(st, image_id); ct <- .store_column(st, cell_type)
  level_cols <- unique(c(condition, subject))   # constant within images (checked) and read at each image's first row
  cl <- lapply(stats::setNames(level_cols, level_cols), .store_column, st = st)
  for (cc in c(list(ci, ct), cl)) if (cc$n != N) stop("column '", cc$name, "' has ", cc$n, " rows, not ", N, ".", call. = FALSE)

  # string (non-categorical) image IDs and cell types: levels in order of first appearance, as for a data.frame
  levi <- new.env(); levt <- new.env(); levi$values <- levt$values <- NULL
  if (ci$kind != "categorical" || ct$kind != "categorical") {
    for (a in seq.int(1L, N, by = as.integer(.store_block()))) {
      b <- as.integer(min(N, a + .store_block() - 1))
      if (ci$kind != "categorical") .store_codes(ci, a, b, levi)
      if (ct$kind != "categorical") .store_codes(ct, a, b, levt)
    }
  }
  li <- if (ci$kind == "categorical") ci$levels else as.character(levi$values)
  lt <- if (ct$kind == "categorical") ct$levels else as.character(levt$values)
  nli <- length(li); nlt <- length(lt)

  # one pass: cells and first row per (image level, type level), runs of images, image-level values at the first row
  counts <- matrix(0L, nli, nlt + 1L); first_rows <- matrix(NA_integer_, nli, nlt + 1L)
  runs <- 0; prev <- NA_integer_
  ref_val <- lapply(cl, function(cc) NULL); constant <- vapply(cl, function(cc) TRUE, TRUE)
  seen <- logical(nli); sub_all <- NULL
  for (a in seq.int(1L, N, by = as.integer(.store_block()))) {
    b <- as.integer(min(N, a + .store_block() - 1))
    im <- .store_codes(ci, a, b, levi); ty <- .store_codes(ct, a, b, levt)
    sc <- cell_block_scan(im, attr(im, "base"), nli, ty, attr(ty, "base"), nlt, a - 1L, prev)
    rm(ty)
    counts <- counts + sc$counts
    fill <- is.na(first_rows) & !is.na(sc$first); first_rows[fill] <- sc$first[fill]
    runs <- runs + sc$runs; prev <- sc$last
    new <- which(!seen & !is.na(sc$pos))
    for (nm in names(cl)) {
      v <- .store_values(cl[[nm]], a, b)
      if (is.null(ref_val[[nm]])) ref_val[[nm]] <- v[rep(NA_integer_, nli)]
      ref_val[[nm]][new] <- v[sc$pos[new]]
      if (constant[[nm]] && !.same_value(v, im, ref_val[[nm]])) constant[[nm]] <- FALSE
      if (!is.null(subject) && nm == subject) sub_all <- unique(c(sub_all, unique(v)))
    }
    seen[new] <- TRUE
  }
  ref_val <- Map(.store_as_values, cl, ref_val)
  if (!is.null(subject)) sub_all <- unique(as.character(.store_as_values(cl[[subject]], sub_all)))
  lay <- list(counts = counts, first = first_rows)

  image_labels <- sort(li[rowSums(lay$counts) > 0], method = "radix")   # C-locale order, as the Python twin
  n_images <- length(image_labels)
  first_level <- apply(lay$first, 1, function(z) if (all(is.na(z))) NA_integer_ else min(z, na.rm = TRUE))
  il <- match(image_labels, li)
  first <- first_level[il]
  # image-level columns at each image's first row (the condition and subject from the pass)
  pheno_cols <- unique(c(condition, subject, covariates))
  pheno <- lapply(stats::setNames(pheno_cols, pheno_cols), function(nm)
    if (nm %in% names(ref_val)) ref_val[[nm]][il] else {
      cc <- .store_column(st, nm); o <- order(first); v <- .store_read(cc, first[o]); v[order(o)] })
  pheno <- structure(pheno, class = "data.frame", row.names = first)

  ctx <- list(pheno = pheno, image_labels = image_labels, n_images = n_images, first = first)
  # the subjects over all the cells: unique(as.character()) of the pass
  ctx <- .cell_units(ctx, if (!is.null(subject)) pheno[[subject]], length(sub_all))
  if (!is.null(condition))
    ctx <- .cell_condition(ctx, pheno[[condition]], condition, ref, check = function()
      if (!constant[[condition]]) stop("'", condition, "' must be constant within each image.", call. = FALSE))

  # types present (the last column of the layout is a missing type), by first appearance with the cells sorted by image
  tcol <- which(colSums(lay$counts) > 0)
  key <- apply(lay$first[, tcol, drop = FALSE], 2, function(z) { i <- which(!is.na(z))[1]; c(i, z[i]) })
  tcol <- tcol[do.call(order, as.data.frame(t(key)))]
  ctx$type_labels <- c(lt, NA_character_)[tcol]
  type_map <- rep(-1L, nlt + 1L); type_map[tcol] <- seq_along(tcol) - 1L
  tb <- table(factor(integer(0), levels = seq_len(n_images) - 1L), factor(integer(0), levels = seq_along(tcol) - 1L))
  tb[] <- lay$counts[il, tcol, drop = FALSE]
  ctx$counts <- unclass(tb)
  # cells per image and type level (all levels; strings as factors with levels in order of first appearance)
  nc <- ctx$counts[match(li, image_labels), match(lt, ctx$type_labels), drop = FALSE]; nc[is.na(nc)] <- 0L
  ctx$n_cells <- structure(array(as.integer(nc), dim = dim(nc), dimnames = stats::setNames(list(li, lt), c("", ""))), class = "table")

  # how to read the cells of a batch of images: each image a run of rows (one block of rows per batch), or a scan
  ctx$store <- list(st = st, ci = ci, ct = ct, coords = coords, type_map = type_map, levi = levi, levt = levt,
                    contiguous = runs == n_images, image_level = il, n_levels = nli, n_rows = N)
  ctx
}

## The neighbour totals of every radius (or k) from the store, a batch of whole images at a time; the same list as
## .cell_stream() gives for the cells in memory.
.store_stream <- function(ctx, radii = NULL, k = NULL, allocation = TRUE, n_threads = 1L) {
  s <- ctx$store
  n <- ctx$n_images; Tn <- length(ctx$type_labels)
  K <- if (is.null(k)) length(radii) else 1L
  A <- lapply(seq_len(K), function(j) array(0, c(Tn, Tn, n)))
  B <- lapply(seq_len(K), function(j) if (allocation) matrix(0, Tn, n) else array(0, c(Tn, Tn, n)))
  size <- rowSums(ctx$counts)
  budget <- min(.store_batch() * max(1L, as.integer(n_threads)), max(.store_batch(), .store_batch_max()))
  # batches of consecutive images (in row order when each image is a run of rows)
  ord <- if (s$contiguous) order(ctx$first) else seq_len(n)
  batch <- integer(n); acc <- 0; id <- 1L
  for (j in seq_along(ord)) {
    if (acc > 0 && acc + size[ord[j]] > budget) { id <- id + 1L; acc <- 0 }
    batch[j] <- id; acc <- acc + size[ord[j]]
  }
  tmap <- function(ty) {   # codes (base 0 or 1, missing as -1 or NA) to the context's 0-based types
    base <- attr(ty, "base"); ty <- as.integer(ty) + (1L - base); ty[is.na(ty) | ty < 1L] <- length(s$type_map)
    s$type_map[ty] }
  for (bi in unique(batch)) {
    imgs <- ord[batch == bi]
    if (s$contiguous) {
      a <- ctx$first[imgs[1]]; b <- a + sum(size[imgs]) - 1L
      # read a block of rows at a time into the batch (the reader's working memory grows with the rows read)
      blk <- as.integer(.store_block()); n_b <- b - a + 1L
      m <- matrix(0, n_b, 2L); ty <- integer(n_b)
      for (u in seq.int(a, b, by = blk)) {
        w <- as.integer(min(b, u + blk - 1L)); at <- (u - a + 1L):(w - a + 1L)
        z <- .store_coords(s$st, s$coords, u, w)
        m[at, 1L] <- z$x; m[at, 2L] <- z$y
        ty[at] <- tmap(.store_codes(s$ct, u, w, s$levt))
      }
      xy <- list(x = m, y = NULL); rm(m, z)
      invisible(gc(verbose = FALSE))   # the reader's working memory, before the neighbour sums
    } else {
      # scan the rows for the batch's images, then sort them stably by image
      img_of_level <- rep(NA_integer_, s$n_levels); img_of_level[s$image_level[imgs]] <- seq_along(imgs)
      parts <- list()
      for (a in seq.int(1L, s$n_rows, by = as.integer(.store_block()))) {
        b <- as.integer(min(s$n_rows, a + .store_block() - 1))
        ic <- .store_codes(s$ci, a, b, s$levi); im <- img_of_level[as.integer(ic) + (1L - attr(ic, "base"))]
        keep <- which(!is.na(im))
        if (!length(keep)) next
        xy <- .store_coords(s$st, s$coords, a, b)
        ty <- tmap(.store_codes(s$ct, a, b, s$levt))
        parts[[length(parts) + 1L]] <- list(im = im[keep], x = xy$x[keep], y = xy$y[keep], ty = ty[keep])
      }
      im <- unlist(lapply(parts, `[[`, "im")); o <- order(im, method = "radix")
      xy <- list(x = unlist(lapply(parts, `[[`, "x"))[o], y = unlist(lapply(parts, `[[`, "y"))[o])
      ty <- unlist(lapply(parts, `[[`, "ty"))[o]
      rm(parts, im, o)
    }
    res <- cell_batch_totals(as.numeric(xy$x), if (is.null(xy$y)) numeric(0) else as.numeric(xy$y), ty, as.integer(c(0L, cumsum(size[imgs]))), Tn,
                             as.numeric(if (is.null(radii)) numeric(0) else radii), if (is.null(k)) 0L else as.integer(k),
                             allocation, as.integer(n_threads))
    rm(xy, ty)
    for (j in seq_len(K)) {
      A[[j]][, , imgs] <- res[[j]]$a
      if (allocation) B[[j]][, imgs] <- res[[j]]$b else B[[j]][, , imgs] <- res[[j]]$b
    }
  }
  lapply(seq_len(K), function(j)
    if (allocation) list(any = A[[j]], self_expected = B[[j]]) else list(totals = A[[j]], sq = B[[j]]))
}

## A SpatialData object (spatialdataR): its table's store when it was read from one (list(path, table = NULL)), else the
## table as a data.frame of its colData, its region column and the coordinates (reducedDim "spatial" as x and y, with
## `coords` TRUE). `info`: the name of a per-image table of the object, returned as a data.frame of its colData.
.spatialdata_input <- function(x, table = NULL, info = NULL) {
  tabs <- methods::slot(x, "tables")
  if (!length(tabs)) stop("the SpatialData object has no table; give the path of its store instead.", call. = FALSE)
  if (is.null(table)) table <- if (length(tabs) == 1L) names(tabs) else if ("table" %in% names(tabs)) "table" else
    stop("the SpatialData object has several tables (", paste(names(tabs), collapse = ", "), "); choose one with `table`.",
         call. = FALSE)
  if (!table %in% names(tabs)) stop("table '", table, "' not found; the object has: ", paste(names(tabs), collapse = ", "),
                                    call. = FALSE)
  sce <- tabs[[table]]
  .need("SingleCellExperiment", "to read a SpatialData object's table")
  meta <- SingleCellExperiment::int_metadata(sce)
  src <- meta$source_path
  if (is.character(src) && length(src) == 1L && dir.exists(src) && requireNamespace("Rarr", quietly = TRUE))
    return(list(path = src, table = NULL))
  d <- as.data.frame(SummarizedExperiment::colData(sce))
  attrs <- meta$spatialdata_attrs
  icd <- SingleCellExperiment::int_colData(sce)
  for (k in c(attrs$region_key, attrs$instance_key)) if (!k %in% names(d) && k %in% names(icd)) d[[k]] <- icd[[k]]
  xy <- if ("spatial" %in% SingleCellExperiment::reducedDimNames(sce)) SingleCellExperiment::reducedDim(sce, "spatial")
  if (!is.null(xy) && !any(c("x", "y") %in% names(d))) { d$x <- as.vector(xy[, 1]); d$y <- as.vector(xy[, 2]) }
  if (is.character(info) && length(info) == 1L) {   # imageInfo: another table of the object, by name
    if (!info %in% names(tabs)) stop("table '", info, "' not found; the object has: ", paste(names(tabs), collapse = ", "),
                                     call. = FALSE)
    info <- as.data.frame(SummarizedExperiment::colData(tabs[[info]]))
  }
  list(cells = d, region_key = attrs$region_key, info = info)
}

#' Write cells as a SpatialData Zarr store
#'
#' Writes a `data.frame` of cells as the cell table of a SpatialData Zarr store (an AnnData table under `tables/`, in
#' Zarr v3), which [spicy()] reads a block of rows at a time and spatialdataR's `readSpatialData()` reads as a
#' `SpatialData` object. The coordinates become the `obsm` array `"spatial"`, the other columns the table's `obs`
#' (factors as categorical columns), and the image column the table's region key (always a categorical column, its
#' categories the levels of a factor, else the images in order of first appearance). A table with one row per image can
#' be added as a second table. Each column is written in chunks of 2^20 rows compressed with zstd, one chunk at a
#' time, so writing needs little memory beyond that of `cells`, however many cells it has.
#'
#' @param cells A `data.frame` with one row per cell. Its columns are factors, character, integer, double or
#'   logical vectors.
#' @param path The directory of the store (replaced if it exists).
#' @param imageID,spatialCoords The column of the image and the two columns of the coordinates.
#' @param tableKey The name of the cell table.
#' @param imageInfo Optional: a `data.frame` with one row per image (row names: the image IDs), written as the table
#'   `infoKey`.
#' @param infoKey The name of the per-image table.
#' @return `path`, invisibly.
#' @examples
#' cells <- data.frame(imageID = factor(rep(c("a", "b"), each = 50)), cellType = factor(rep(c("T", "B"), 50)),
#'                     x = runif(100, 0, 100), y = runif(100, 0, 100))
#' if (requireNamespace("Rarr", quietly = TRUE)) {
#'   path <- writeCellStore(cells, tempfile(fileext = ".zarr"))
#'   list.files(file.path(path, "tables"))
#' }
#' @export
writeCellStore <- function(cells, path, imageID = "imageID", spatialCoords = c("x", "y"), tableKey = "table",
                           imageInfo = NULL, infoKey = "images") {
  .need("Rarr", "to write a SpatialData store")
  cells <- as.data.frame(cells)
  miss <- setdiff(c(imageID, spatialCoords), names(cells))
  if (length(miss)) stop("columns not found in `cells`: ", paste(miss, collapse = ", "), call. = FALSE)
  obs <- cells[setdiff(names(cells), spatialCoords)]
  .zarr_check_columns(obs, "cells")
  if (!is.null(imageInfo)) { imageInfo <- as.data.frame(imageInfo); .zarr_check_columns(imageInfo, "imageInfo") }
  unlink(path, recursive = TRUE); dir.create(file.path(path, "tables"), recursive = TRUE)
  path <- normalizePath(path)
  Rarr::write_zarr_group(path, "", zarr_version = 3L)
  Rarr::write_zarr_group(path, "tables", zarr_version = 3L)
  Rarr::write_zarr_attributes(path, new.zattrs = list(spatialdata_attrs = list(version = "0.2")))
  n <- nrow(cells)
  img <- cells[[imageID]]
  regions <- if (is.factor(img)) integer(0) else character(0)   # in order of first appearance, a chunk at a time
  for (r in .zarr_chunks(n)) {
    u <- unique(if (is.factor(img)) .subset(img, r[1]:r[2]) else as.character(img[r[1]:r[2]]))
    regions <- c(regions, u[!u %in% regions])
    .zarr_free()
  }
  if (is.factor(img)) regions <- levels(img)[regions]
  # the region key as a categorical (as spatialdata requires): a factor keeps its levels, other IDs take theirs in
  # order of first appearance
  else obs[[imageID]] <- structure(match(as.character(img), regions), levels = regions, class = "factor")
  attrs <- list(region = regions, region_key = imageID, instance_key = "instance_id")
  cx <- cells[[spatialCoords[1]]]; cy <- cells[[spatialCoords[2]]]
  .zarr_write_anndata(file.path(path, "tables", tableKey), obs, n,
                      index = function(a, b) .vlen_ids((a:b) - 1L),
                      instance = TRUE,
                      spatial = function(a, b) as.vector(rbind(as.numeric(cx[a:b]), as.numeric(cy[a:b]))),
                      spatialdata_attrs = attrs,
                      attrs = c(list(`spatialdata-encoding-type` = "ngff:regions_table", version = "0.2"), attrs))
  if (!is.null(imageInfo)) {
    ids <- rownames(imageInfo)
    .zarr_write_anndata(file.path(path, "tables", infoKey), imageInfo, nrow(imageInfo),
                        index = function(a, b) ids[a:b],
                        attrs = list(`spatialdata-encoding-type` = "ngff:regions_table", version = "0.2"))
  }
  Rarr::zarr_consolidate_metadata(path)
  invisible(path)
}

## ---- writing a store ---------------------------------------------------------------------------------------------
## writeCellStore() writes each AnnData table in the encoding anndata gives it in a Zarr v3 store. Every array of the
## length of the table is cut into chunks of .store_write_chunk() rows; each chunk is taken from the input, encoded
## (little-endian numbers, or vlen-utf8 strings), compressed with zstd and written before the next is made, so no
## whole column is ever copied.

.store_write_chunk <- function() as.integer(getOption("spicyR.store_write_chunk", 2^20))

## Frees the last chunk's temporary vectors (a collection of the youngest generation): without it they would pile
## up until R's next full collection, whose threshold grows with the size of `cells`.
.zarr_free <- function() invisible(gc(verbose = FALSE, full = FALSE))

## The first and last rows of each chunk of n rows.
.zarr_chunks <- function(n, ch = .store_write_chunk()) {
  if (n < 1) return(list())
  a <- seq.int(1L, n, by = ch)
  .mapply(c, list(a, pmin(a + ch - 1L, n)), NULL)
}

## Columns of a type that can be written.
.zarr_check_columns <- function(d, what) {
  ok <- vapply(d, function(v) is.factor(v) || is.character(v) || is.logical(v) || is.numeric(v), TRUE)
  if (!all(ok)) stop("columns of `", what, "` that are not factors, character, numbers or logical: ",
                     paste(names(d)[!ok], collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

## The smallest integer type holding lo..hi (as anndata writes integers).
.zarr_int_type <- function(lo, hi) {
  if (lo < 0) return(if (lo >= -128 && hi <= 127) "int8" else if (lo >= -32768 && hi <= 32767) "int16" else "int32")
  if (hi <= 255) "uint8" else if (hi <= 65535) "uint16" else "uint32"
}

## A Zarr v3 group, with attributes.
.zarr_group <- function(p, attrs = NULL) {
  Rarr::write_zarr_group(dirname(p), basename(p), zarr_version = 3L)
  if (length(attrs)) Rarr::write_zarr_attributes(p, new.zattrs = attrs, zarr_version = 3L)
  invisible(p)
}

## A Zarr v3 array of n rows (and `ncol` columns) of type `type`, written a chunk of rows at a time: get(a, b) gives
## rows a..b as a vector (in row-major order when there are columns), or strings already in the vlen-utf8 encoding.
## `encoding`: its anndata encoding-type and -version, or NULL for none.
.zarr_write_rows <- function(p, n, type, get, ncol = NULL, encoding = c("array", "0.2.0")) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  ch <- max(1L, min(.store_write_chunk(), n))
  num <- function(x) paste(sprintf("%.0f", as.double(x)), collapse = ",")
  string <- type == "string"
  meta <- paste0(
    '{"zarr_format":3,"node_type":"array","shape":[', num(c(n, ncol)), '],"data_type":"', type, '",',
    '"chunk_grid":{"name":"regular","configuration":{"chunk_shape":[', num(c(ch, ncol)), ']}},',
    '"chunk_key_encoding":{"name":"default","configuration":{"separator":"/"}},',
    '"fill_value":', if (string) '""' else if (type == "bool") "false" else "0", ',',
    '"codecs":[', if (string) '{"name":"vlen-utf8"}' else '{"name":"bytes","configuration":{"endian":"little"}}',
    ',{"name":"zstd","configuration":{"level":0,"checksum":false}}]',
    if (!is.null(encoding)) paste0(',"attributes":{"encoding-type":"', encoding[1], '","encoding-version":"',
                                   encoding[2], '"}'), '}')
  writeLines(meta, file.path(p, "zarr.json"))
  size <- switch(type, bool = , int8 = , uint8 = 1L, int16 = , uint16 = 2L, int32 = , uint32 = 4L, float64 = 8L,
                 string = NA_integer_)
  per <- if (is.null(ncol)) 1L else as.integer(ncol)
  for (k in seq_along(rows <- .zarr_chunks(n, ch))) {
    v <- get(rows[[k]][1], rows[[k]][2])
    m <- length(v)
    if (is.raw(v)) {   # strings already encoded (.vlen_ids()): pad with empty strings to a full chunk
      m <- rows[[k]][2] - rows[[k]][1] + 1L
      bytes <- c(v, raw(4L * (ch - m)))
      bytes[1:4] <- writeBin(ch, raw(), size = 4L, endian = "little")
    } else {
      if (m < ch * per) v[(m + 1L):(ch * per)] <- vector(typeof(v), 1L)   # a full chunk: the fill value after the rows
      bytes <- if (string) .vlen_utf8(v) else
        writeBin(if (is.logical(v)) as.integer(v) else v, raw(), size = size, endian = "little")
    }
    f <- file.path(p, "c", if (is.null(ncol)) k - 1L else file.path(k - 1L, 0L))
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    writeBin(memCompress(bytes, "zstd"), f)
    rm(v, bytes); .zarr_free()
  }
  invisible(p)
}

## A string scalar: a zero-dimensional string array (anndata's encoding of a string in uns).
.zarr_write_string <- function(p, s) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  writeLines(paste0('{"zarr_format":3,"node_type":"array","shape":[],"data_type":"string",',
                    '"chunk_grid":{"name":"regular","configuration":{"chunk_shape":[]}},',
                    '"chunk_key_encoding":{"name":"default","configuration":{"separator":"/"}},"fill_value":"",',
                    '"codecs":[{"name":"vlen-utf8"},{"name":"zstd","configuration":{"level":0,"checksum":false}}],',
                    '"attributes":{"encoding-type":"string","encoding-version":"0.2.0"}}'), file.path(p, "zarr.json"))
  writeBin(memCompress(.vlen_utf8(s), "zstd"), file.path(p, "c"))
  invisible(p)
}

## Strings in the vlen-utf8 encoding: their number, then each string's length in bytes and its bytes. writeBin()
## gives the bytes of each string followed by a NUL; each byte moves 3 bytes further on per earlier string (whose
## 4-byte length replaces its 1-byte NUL), and the NUL before each string's bytes becomes the first of its length.
.vlen_utf8 <- function(s) {
  s <- enc2utf8(s)
  n <- length(s)
  len <- nchar(s, type = "bytes")
  w <- writeBin(s, raw())
  out <- raw(4L + 4L * n + length(w) - n)
  if (n) {
    k <- seq_len(length(w) - 1L)   # the NUL after the last string has no place
    out[k + 3L * rep.int(seq_len(n), len + 1L)[k] + 5L] <- w[k]
    at <- 4L + 4L * seq_len(n) - 3L + c(0L, cumsum(len)[-n])   # where each string's length goes
    out[rep(at, each = 4L) + 0:3] <- writeBin(len, raw(), size = 4L, endian = "little")
  }
  out[1:4] <- writeBin(n, raw(), size = 4L, endian = "little")
  out
}

## The vlen-utf8 encoding of the decimal strings of the integers `i` (>= 0), as .vlen_utf8(as.character(i)) gives
## it, made from the digits without making the strings (the row index of a table: 0, 1, ...).
.vlen_ids <- function(i) {
  n <- length(i)
  nd <- rep.int(1L, n)
  for (k in 1:9) nd <- nd + (i >= 10L^k)   # the number of digits
  at <- 4L + 4L * seq_len(n) - 3L + c(0L, cumsum(nd)[-n])   # where each string's length goes
  out <- raw(4L + 4L * n + sum(nd))
  out[1:4] <- writeBin(n, raw(), size = 4L, endian = "little")
  if (!n) return(out)
  out[rep(at, each = 4L) + 0:3] <- writeBin(nd, raw(), size = 4L, endian = "little")
  last <- at + 3L + nd   # the last digit's byte
  for (k in 0:9) {
    if (all(w <- nd > k)) w <- TRUE else if (!any(w)) break
    out[last[w] - k] <- as.raw(48L + (i[w] %/% as.integer(10^k)) %% 10L)
  }
  out
}

## One obs column, as anndata encodes it: a factor as a categorical (codes and categories), a vector with missing
## values as a nullable array (values and mask), any other vector as an array.
.zarr_write_column <- function(p, x, n) {
  if (is.factor(x)) {
    .zarr_group(p, list(`encoding-type` = "categorical", `encoding-version` = "0.2.0", ordered = is.ordered(x)))
    lev <- levels(x)
    .zarr_write_rows(file.path(p, "categories"), length(lev), "string", function(a, b) lev[a:b],
                     encoding = c("string-array", "0.2.0"))
    type <- .zarr_int_type(if (anyNA(x)) -1 else 0, length(lev) - 1)
    .zarr_write_rows(file.path(p, "codes"), n, type, function(a, b) {
      v <- .subset(x, a:b) - 1L; v[is.na(v)] <- -1L; v
    })
    return(invisible(p))
  }
  if (is.character(x)) {
    type <- "string"; fill <- ""; enc <- c("string-array", "0.2.0"); nullable <- "nullable-string-array"
  } else if (is.logical(x)) {
    type <- "bool"; fill <- FALSE; enc <- c("array", "0.2.0"); nullable <- "nullable-boolean"
  } else if (is.integer(x)) {
    na <- anyNA(x)
    rng <- if (n && !all(is.na(x))) c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)) else c(0L, 0L)
    type <- .zarr_int_type(if (na) min(rng[1], -1L) else rng[1], rng[2])
    fill <- -1L; enc <- c("array", "0.2.0"); nullable <- "nullable-integer"
  } else {
    .zarr_write_rows(p, n, "float64", function(a, b) as.double(x[a:b]))
    return(invisible(p))
  }
  if (!anyNA(x)) return(invisible(.zarr_write_rows(p, n, type, function(a, b) x[a:b], encoding = enc)))
  .zarr_group(p, c(list(`encoding-type` = nullable, `encoding-version` = "0.1.0"),
                   if (is.character(x)) list(`na-value` = "NA")))
  .zarr_write_rows(file.path(p, "values"), n, type, function(a, b) { v <- x[a:b]; v[is.na(v)] <- fill; v },
                   encoding = enc)
  .zarr_write_rows(file.path(p, "mask"), n, "bool", function(a, b) is.na(x[a:b]))
  invisible(p)
}

## An AnnData table of n rows at `tab`: obs holds the columns of `obs` (and, with `instance`, 0-based instance IDs as
## "instance_id"), its index given a chunk at a time by index(a, b); obsm "spatial" by spatial(a, b) (row-major x
## and y); uns "spatialdata_attrs" from `spatialdata_attrs`; `attrs` added to the table's attributes.
.zarr_write_anndata <- function(tab, obs, n, index, instance = FALSE, spatial = NULL, spatialdata_attrs = NULL,
                                attrs = list()) {
  .zarr_group(tab, c(list(`encoding-type` = "anndata", `encoding-version` = "0.1.0"), attrs))
  cols <- c(names(obs), if (instance) "instance_id")
  op <- file.path(tab, "obs")
  .zarr_group(op, list(`encoding-type` = "dataframe", `encoding-version` = "0.2.0", `_index` = "_index",
                       `column-order` = as.list(cols)))
  .zarr_write_rows(file.path(op, "_index"), n, "string", index, encoding = c("string-array", "0.2.0"))
  for (cn in names(obs)) .zarr_write_column(file.path(op, cn), obs[[cn]], n)
  # instance IDs in 16 bits at least (spatialdata refuses an 8-bit instance key)
  if (instance) .zarr_write_rows(file.path(op, "instance_id"), n, .zarr_int_type(0, max(n - 1, 256)),
                                 function(a, b) (a:b) - 1L)
  vp <- file.path(tab, "var")
  .zarr_group(vp, list(`encoding-type` = "dataframe", `encoding-version` = "0.2.0", `_index` = "_index",
                       `column-order` = list()))
  .zarr_write_rows(file.path(vp, "_index"), 0L, "string", NULL, encoding = c("string-array", "0.2.0"))
  dict <- list(`encoding-type` = "dict", `encoding-version` = "0.1.0")
  for (g in c("layers", "obsm", "obsp", "uns", "varm", "varp")) .zarr_group(file.path(tab, g), dict)
  if (!is.null(spatial)) .zarr_write_rows(file.path(tab, "obsm", "spatial"), n, "float64", spatial, ncol = 2L)
  if (!is.null(spatialdata_attrs)) {
    sp <- file.path(tab, "uns", "spatialdata_attrs")
    .zarr_group(sp, dict)
    reg <- spatialdata_attrs$region
    .zarr_write_rows(file.path(sp, "region"), length(reg), "string", function(a, b) reg[a:b],
                     encoding = c("string-array", "0.2.0"))
    for (k in c("region_key", "instance_key")) .zarr_write_string(file.path(sp, k), spatialdata_attrs[[k]])
  }
  invisible(tab)
}
