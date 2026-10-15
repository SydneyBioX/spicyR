## ---- the cell method from a SpatialData (Zarr) store ------------------------------------------------------------
## A SpatialData store's cell table is an AnnData group (tables/<name>) whose columns are Zarr arrays, read here with
## Rarr. Nothing of the length of the table is held in R: one pass over the image, cell type, condition and subject
## columns, a block of rows at a time, gives the cells per image and type, the first row of each image (and of each
## image and type) and the image-level values; the coordinates and cell types are then read a batch of whole images
## at a time for the neighbour sums (cell_batch_totals). The context has the fields .cell_context() gives a
## data.frame of the same cells, so every later step is shared.

## Whether `x` names a Zarr store (a directory).
.is_store <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && dir.exists(x)

## Rows read per block in the pass over the columns, and cells per thread in a batch of whole images.
.store_block <- function() getOption("spicyR.store_block", 2^23)
.store_batch <- function() getOption("spicyR.store_batch", 2^24)

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

## Integer codes (1-based, NA where missing) of rows a..b of a categorical column; for any other column, the codes
## against `lev`, a growing record of the values in order of first appearance (an environment).
.store_codes <- function(col, a, b, lev = NULL) {
  if (col$kind == "categorical") {
    v <- as.integer(as.vector(Rarr::read_zarr_array(col$path, index = list(a:b)))) + 1L
    v[v < 1L] <- NA_integer_
    return(v)
  }
  v <- .store_read(col, a:b)
  u <- unique(v); u <- u[!is.na(u)]
  lev$values <- c(lev$values, u[!u %in% lev$values])
  match(v, lev$values)
}

## Coordinates of rows a..b: an obsm array (first two columns) or two obs columns.
.store_coords <- function(st, coords, a, b) {
  if (coords$obsm) {
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
    return(list(obsm = TRUE, path = p, label = spatialCoords))
  }
  list(obsm = FALSE, x = .store_column(st, spatialCoords[1]), y = .store_column(st, spatialCoords[2]),
       label = spatialCoords)
}

## Whether v (rows of the table) equals, as as.character() sees it, ref at each row's image (codes im); NA equals NA.
.same_value <- function(v, im, ref) {
  if (is.factor(v)) { v <- unclass(v); attributes(v) <- NULL; ref <- unclass(ref); attributes(ref) <- NULL }   # codes
  d <- cell_rows_differing(v, im, ref)
  !length(d) || identical(as.character(v[d]), as.character(ref[im[d]]))
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
    im <- .store_codes(ci, a, b, levi)
    if (anyNA(im)) stop("missing image ID in row ", a - 1L + which(is.na(im))[1], ".", call. = FALSE)
    ty <- .store_codes(ct, a, b, levt)
    lay <- cell_layout(im, nli, ty, nlt)
    counts <- counts + lay$counts
    fr <- lay$first + (a - 1L); fill <- is.na(first_rows) & !is.na(fr); first_rows[fill] <- fr[fill]
    n <- length(im)
    runs <- runs + sum(im[-1L] != im[-n]) + (is.na(prev) || prev != im[1L]); prev <- im[n]
    new <- which(!seen & tabulate(im, nli) > 0L); pos <- match(new, im)
    for (nm in names(cl)) {
      v <- .store_read(cl[[nm]], a:b)
      if (is.null(ref_val[[nm]])) ref_val[[nm]] <- v[rep(NA_integer_, nli)]
      ref_val[[nm]][new] <- v[pos]
      if (constant[[nm]] && !.same_value(v, im, ref_val[[nm]])) constant[[nm]] <- FALSE
      if (!is.null(subject) && nm == subject && !constant[[nm]])
        sub_all <- unique(c(sub_all, unique(as.character(v))))
    }
    seen[new] <- TRUE
  }
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

  if (!is.null(subject)) {
    sv <- pheno[[subject]]
    n_sub <- if (constant[[subject]]) length(unique(as.character(sv))) else length(sub_all)
  }
  if (is.null(subject) || n_sub == n_images) {
    unit_labels <- image_labels; image_unit <- seq_len(n_images) - 1L
  } else {
    sub <- as.character(sv); unit_labels <- unique(sub); image_unit <- match(sub, unit_labels) - 1L
  }
  ctx <- list(pheno = pheno, image_labels = image_labels, n_images = n_images, unit_labels = unit_labels,
              image_unit = image_unit, first = first)
  cv <- pheno[[condition]]
  lev <- .cell_levels(cv)
  if (length(lev) < 2L) stop("`condition` needs at least two levels; found ", length(lev), ".", call. = FALSE)
  if (!is.null(ref)) {
    if (!as.character(ref) %in% lev) stop("`ref` is not a level of `condition`.", call. = FALSE)
    lev <- c(as.character(ref), setdiff(lev, as.character(ref)))
  }
  if (!constant[[condition]]) stop("'", condition, "' must be constant within each image.", call. = FALSE)
  ctx$levels <- lev
  ctx$image_group <- match(as.character(cv), lev) - 1L
  if (any(tapply(ctx$image_group, image_unit, function(z) length(unique(z))) > 1L))
    stop("each subject must belong to a single '", condition, "' level.", call. = FALSE)

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
  budget <- .store_batch() * max(1L, as.integer(n_threads))
  # batches of consecutive images (in row order when each image is a run of rows)
  ord <- if (s$contiguous) order(ctx$first) else seq_len(n)
  batch <- integer(n); acc <- 0; id <- 1L
  for (j in seq_along(ord)) {
    if (acc > 0 && acc + size[ord[j]] > budget) { id <- id + 1L; acc <- 0 }
    batch[j] <- id; acc <- acc + size[ord[j]]
  }
  tmap <- function(ty) { ty[is.na(ty)] <- length(s$type_map); s$type_map[ty] }
  for (bi in unique(batch)) {
    imgs <- ord[batch == bi]
    if (s$contiguous) {
      a <- ctx$first[imgs[1]]; b <- a + sum(size[imgs]) - 1L
      xy <- .store_coords(s$st, s$coords, a, b)
      ty <- tmap(.store_codes(s$ct, a, b, s$levt))
    } else {
      # scan the rows for the batch's images, then sort them stably by image
      img_of_level <- rep(NA_integer_, s$n_levels); img_of_level[s$image_level[imgs]] <- seq_along(imgs)
      parts <- list()
      for (a in seq.int(1L, s$n_rows, by = as.integer(.store_block()))) {
        b <- as.integer(min(s$n_rows, a + .store_block() - 1))
        im <- img_of_level[.store_codes(s$ci, a, b, s$levi)]
        keep <- which(!is.na(im))
        if (!length(keep)) next
        xy <- .store_coords(s$st, s$coords, a, b)
        ty <- .store_codes(s$ct, a, b, s$levt)
        parts[[length(parts) + 1L]] <- list(im = im[keep], x = xy$x[keep], y = xy$y[keep], ty = ty[keep])
      }
      im <- unlist(lapply(parts, `[[`, "im")); o <- order(im, method = "radix")
      xy <- list(x = unlist(lapply(parts, `[[`, "x"))[o], y = unlist(lapply(parts, `[[`, "y"))[o])
      ty <- tmap(unlist(lapply(parts, `[[`, "ty"))[o])
      rm(parts, im, o)
    }
    res <- cell_batch_totals(as.numeric(xy$x), as.numeric(xy$y), ty, as.integer(c(0L, cumsum(size[imgs]))), Tn,
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
