## ---- writing a SpatialData store ---------------------------------------------------------------------------------
## A SpatialExperiment, SingleCellExperiment or data.frame of cells written as a SpatialData Zarr store: Zarr v3, the
## SpatialData 0.2 container, tables in anndata's encoding and images as OME-NGFF image elements. Every array is cut
## into chunks; each chunk is taken from the input, encoded, compressed with zstd and written before the next is made,
## so writing needs little memory beyond that of the input, however many cells it has.
##
## The code follows the idiom of spatialdataR's draft writers (writeTable(x, name, path, ...), writeImage(...)), so
## that it can move there: spatialdataR does not yet write tables in Zarr v3. Numbers and logicals are written with
## Rarr; strings are written here, in the vlen-utf8 encoding that anndata and Python's zarr use and Rarr does not.

.store_write_chunk <- \() as.integer(getOption("spicyR.store_write_chunk", 2^20))

## Frees the last chunk's temporary vectors (a collection of the youngest generation): without it they would pile
## up until R's next full collection, whose threshold grows with the size of the input.
.zarr_free <- \() invisible(gc(verbose=FALSE, full=FALSE))

## The first and last rows of each chunk of n rows.
.zarr_chunks <- function(n, ch=.store_write_chunk()) {
    if (n < 1) return(list())
    a <- seq.int(1L, n, by=ch)
    .mapply(c, list(a, pmin(a + ch - 1L, n)), NULL)
}

## Columns of a type that can be written.
.zarr_check_columns <- function(d, what) {
    ok <- vapply(d, \(v) is.factor(v) || is.character(v) || is.logical(v) || is.numeric(v), TRUE)
    if (!all(ok)) stop("columns of `", what, "` that are not factors, character, numbers or logical: ",
        paste(names(d)[!ok], collapse=", "), call.=FALSE)
    invisible(TRUE)
}

## The smallest integer type holding lo..hi (as anndata writes integers).
.zarr_int_type <- function(lo, hi) {
    if (lo < 0) return(if (lo >= -128 && hi <= 127) "int8" else if (lo >= -32768 && hi <= 32767) "int16" else "int32")
    if (hi <= 255) "uint8" else if (hi <= 65535) "uint16" else "uint32"
}

## Rarr's name of a Zarr data type.
.rarr_type <- \(type) switch(type,
    bool="|b1", int8="|i1", uint8="|u1", int16="<i2", uint16="<u2", int32="<i4", uint32="<u4", float64="<f8",
    stop("no Rarr data type for ", type))

## A Zarr v3 group, with attributes.
.zarr_group <- function(p, attrs=NULL) {
    Rarr::write_zarr_group(dirname(p), basename(p), zarr_version=3L)
    if (length(attrs)) Rarr::write_zarr_attributes(p, new.zattrs=attrs, zarr_version=3L)
    invisible(p)
}

## A Zarr v3 array of n rows (and `ncol` columns) of type `type`, written a chunk of rows at a time: get(a, b) gives
## rows a..b as a vector or a matrix, or strings already in the vlen-utf8 encoding. A matrix is cut into chunks of
## whole rows, about as many values per chunk as a column has. `encoding`: its anndata encoding-type and -version, or
## NULL for none.
.zarr_write_rows <- function(p, n, type, get, ncol=NULL, encoding=c("array", "0.2.0")) {
    if (type == "string") return(.zarr_write_strings(p, n, get, encoding))
    ch <- .store_write_chunk()
    if (!is.null(ncol)) ch <- min(ch, (2L * ch) %/% max(1L, ncol))
    ch <- max(1L, min(ch, n))
    dir.create(dirname(p), recursive=TRUE, showWarnings=FALSE)
    withCallingHandlers(   # an empty array has a chunk larger than itself, which Rarr warns of
        Rarr::create_empty_zarr_array(p, dim=as.integer(c(n, ncol)), chunk_dim=as.integer(c(ch, ncol)),
            data_type=.rarr_type(type), order="C", zarr_version=3L),
        warning=\(w) if (n == 0) invokeRestart("muffleWarning"))
    for (r in .zarr_chunks(n, ch)) {
        v <- get(r[1], r[2])
        if (type == "bool") storage.mode(v) <- "logical"
        idx <- if (is.null(ncol)) list(r[1]:r[2]) else list(r[1]:r[2], seq_len(ncol))
        Rarr::update_zarr_array(p, v, index=idx)
        rm(v); .zarr_free()
    }
    if (!is.null(encoding)) Rarr::write_zarr_attributes(p, new.zattrs=list(
        `encoding-type`=encoding[1], `encoding-version`=encoding[2]), zarr_version=3L)
    invisible(p)
}

## A string array of n rows in the vlen-utf8 encoding, a chunk at a time: get(a, b) gives strings a..b, or their
## encoding already made (.vlen_ids()).
.zarr_write_strings <- function(p, n, get, encoding) {
    dir.create(p, recursive=TRUE, showWarnings=FALSE)
    ch <- max(1L, min(.store_write_chunk(), n))
    writeLines(paste0(
        '{"zarr_format":3,"node_type":"array","shape":[', sprintf("%.0f", n), '],"data_type":"string",',
        '"chunk_grid":{"name":"regular","configuration":{"chunk_shape":[', sprintf("%.0f", ch), ']}},',
        '"chunk_key_encoding":{"name":"default","configuration":{"separator":"/"}},"fill_value":"",',
        '"codecs":[{"name":"vlen-utf8"},{"name":"zstd","configuration":{"level":0,"checksum":false}}]',
        if (!is.null(encoding)) paste0(',"attributes":{"encoding-type":"', encoding[1],
            '","encoding-version":"', encoding[2], '"}'), '}'), file.path(p, "zarr.json"))
    for (k in seq_along(rows <- .zarr_chunks(n, ch))) {
        v <- get(rows[[k]][1], rows[[k]][2])
        m <- rows[[k]][2] - rows[[k]][1] + 1L
        if (is.raw(v)) {   # already encoded: pad with empty strings to a full chunk
            bytes <- c(v, raw(4L * (ch - m)))
            bytes[1:4] <- writeBin(ch, raw(), size=4L, endian="little")
        } else {
            if (m < ch) v[(m + 1L):ch] <- ""   # a full chunk: the fill value after the rows
            bytes <- .vlen_utf8(v)
        }
        dir.create(file.path(p, "c"), showWarnings=FALSE)
        writeBin(memCompress(bytes, "zstd"), file.path(p, "c", k - 1L))
        rm(v, bytes); .zarr_free()
    }
    invisible(p)
}

## A string scalar: a zero-dimensional string array (anndata's encoding of a string in uns).
.zarr_write_string <- function(p, s) {
    dir.create(p, recursive=TRUE, showWarnings=FALSE)
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
    len <- nchar(s, type="bytes")
    w <- writeBin(s, raw())
    out <- raw(4L + 4L * n + length(w) - n)
    if (n) {
        k <- seq_len(length(w) - 1L)   # the NUL after the last string has no place
        out[k + 3L * rep.int(seq_len(n), len + 1L)[k] + 5L] <- w[k]
        at <- 4L + 4L * seq_len(n) - 3L + c(0L, cumsum(len)[-n])   # where each string's length goes
        out[rep(at, each=4L) + 0:3] <- writeBin(len, raw(), size=4L, endian="little")
    }
    out[1:4] <- writeBin(n, raw(), size=4L, endian="little")
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
    out[1:4] <- writeBin(n, raw(), size=4L, endian="little")
    if (!n) return(out)
    out[rep(at, each=4L) + 0:3] <- writeBin(nd, raw(), size=4L, endian="little")
    last <- at + 3L + nd   # the last digit's byte
    for (k in 0:9) {
        if (all(w <- nd > k)) w <- TRUE else if (!any(w)) break
        out[last[w] - k] <- as.raw(48L + (i[w] %/% as.integer(10^k)) %% 10L)
    }
    out
}

## One column of a data frame, as anndata encodes it: a factor as a categorical (codes and categories), a vector with
## missing values as a nullable array (values and mask), any other vector as an array.
.zarr_write_column <- function(p, x, n) {
    if (is.factor(x)) {
        .zarr_group(p, list(`encoding-type`="categorical", `encoding-version`="0.2.0", ordered=is.ordered(x)))
        lev <- levels(x)
        .zarr_write_rows(file.path(p, "categories"), length(lev), "string", \(a, b) lev[a:b],
            encoding=c("string-array", "0.2.0"))
        type <- .zarr_int_type(if (anyNA(x)) -1 else 0, length(lev) - 1)
        .zarr_write_rows(file.path(p, "codes"), n, type, \(a, b) {
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
        rng <- if (n && !all(is.na(x))) c(min(x, na.rm=TRUE), max(x, na.rm=TRUE)) else c(0L, 0L)
        type <- .zarr_int_type(if (na) min(rng[1], -1L) else rng[1], rng[2])
        fill <- -1L; enc <- c("array", "0.2.0"); nullable <- "nullable-integer"
    } else {
        .zarr_write_rows(p, n, "float64", \(a, b) as.double(x[a:b]))
        return(invisible(p))
    }
    if (!anyNA(x)) return(invisible(.zarr_write_rows(p, n, type, \(a, b) x[a:b], encoding=enc)))
    .zarr_group(p, c(list(`encoding-type`=nullable, `encoding-version`="0.1.0"),
        if (is.character(x)) list(`na-value`="NA")))
    .zarr_write_rows(file.path(p, "values"), n, type, \(a, b) { v <- x[a:b]; v[is.na(v)] <- fill; v },
        encoding=enc)
    .zarr_write_rows(file.path(p, "mask"), n, "bool", \(a, b) is.na(x[a:b]))
    invisible(p)
}

## A data frame of n rows (obs or var), its index given a chunk at a time by index(a, b); `extra`: more columns,
## each a function(a, b) with attributes "type" and "name", written after those of `d`.
.zarr_write_dataframe <- function(p, d, n, index, extra=list()) {
    cols <- c(names(d), vapply(extra, attr, "", "name"))
    .zarr_group(p, list(`encoding-type`="dataframe", `encoding-version`="0.2.0", `_index`="_index",
        `column-order`=as.list(cols)))
    .zarr_write_rows(file.path(p, "_index"), n, "string", index, encoding=c("string-array", "0.2.0"))
    for (cn in names(d)) .zarr_write_column(file.path(p, cn), d[[cn]], n)
    for (e in extra) .zarr_write_rows(file.path(p, attr(e, "name")), n, attr(e, "type"), e)
    invisible(p)
}

## An assay (features x cells) as anndata's matrix of cells x features. A sparse matrix is written as a csr_matrix:
## the columns of a dgCMatrix are the rows of the CSR form of its transpose, so its slots are written as they are.
## Any other matrix (dense, or a DelayedMatrix) is written as a dense array, a block of cells at a time.
.zarr_write_matrix <- function(p, x) {
    n <- ncol(x); m <- nrow(x)
    if (.is_class(x, "sparseMatrix")) {
        if (!.is_class(x, "dgCMatrix")) {
            .need("Matrix", "to write a sparse assay")
            x <- methods::as(methods::as(methods::as(x, "CsparseMatrix"), "generalMatrix"), "dMatrix")
        }
        .zarr_group(p, list(`encoding-type`="csr_matrix", `encoding-version`="0.1.0", shape=c(n, m)))
        nnz <- length(x@x)
        .zarr_write_rows(file.path(p, "data"), nnz, "float64", \(a, b) x@x[a:b], encoding=NULL)
        .zarr_write_rows(file.path(p, "indices"), nnz, "int32", \(a, b) x@i[a:b], encoding=NULL)
        .zarr_write_rows(file.path(p, "indptr"), n + 1L, "int32", \(a, b) x@p[a:b], encoding=NULL)
        return(invisible(p))
    }
    type <- switch(typeof(as.matrix(x[0, 0, drop=FALSE])), integer="int32", logical="bool", "float64")
    .zarr_write_rows(p, n, type, \(a, b) t(as.matrix(x[, a:b, drop=FALSE])), ncol=m)
}

## The categorical region key and the regions it names: a factor keeps its levels, other IDs are given theirs in
## order of first appearance (found a chunk at a time).
.region_key <- function(img) {
    if (is.factor(img)) return(list(key=img, regions=levels(img)))
    regions <- character(0)
    for (r in .zarr_chunks(length(img))) {
        u <- unique(as.character(img[r[1]:r[2]]))
        regions <- c(regions, u[!u %in% regions])
        .zarr_free()
    }
    list(key=structure(match(as.character(img), regions), levels=regions, class="factor"), regions=regions)
}

## A table: an AnnData group `name` under `path`/tables. `x` is a SummarizedExperiment (a SingleCellExperiment or
## SpatialExperiment, say), or a
## data.frame with one row per observation. obs holds the colData (or the columns) other than the coordinates, and
## 0-based instance IDs as `instance_key`; X the first assay and layers the others; var the rowData; obsm "spatial"
## the coordinates (a SpatialExperiment's spatialCoords(), or the two `spatialCoords` columns) and the reducedDims.
## With `region_key`, the table annotates the regions named in that column (as a categorical), as SpatialData
## requires; without it, its index is the row names of `x` (a table of images, say).
.writeTable <- function(x, name, path, region_key=NULL, instance_key="instance_id", spatialCoords=NULL) {
    tab <- file.path(path, "tables", name)
    sce <- .is_class(x, "SummarizedExperiment")
    obs <- if (sce) as.data.frame(SummarizedExperiment::colData(x), optional=TRUE) else as.data.frame(x)
    n <- nrow(obs)
    xy <- NULL
    if (!is.null(spatialCoords) && all(spatialCoords %in% names(obs))) {
        xy <- spatialCoords; cx <- obs[[xy[1]]]; cy <- obs[[xy[2]]]
        obs <- obs[setdiff(names(obs), xy)]
        spatial <- \(a, b) cbind(as.numeric(cx[a:b]), as.numeric(cy[a:b]))
    } else if (.is_class(x, "SpatialExperiment")) {
        sc <- SpatialExperiment::spatialCoords(x)
        spatial <- \(a, b) { v <- sc[a:b, 1:2, drop=FALSE]; storage.mode(v) <- "double"; unname(v) }
        xy <- colnames(sc)[1:2]
    } else if (!is.null(spatialCoords)) {
        stop("columns not found: ", paste(setdiff(spatialCoords, names(obs)), collapse=", "), call.=FALSE)
    }
    .zarr_check_columns(obs, if (sce) "colData" else "cells")
    rds <- if (.is_class(x, "SingleCellExperiment")) as.list(SingleCellExperiment::reducedDims(x)) else list()
    if (!is.null(xy) && "spatial" %in% names(rds))
        stop("a reducedDim is called \"spatial\", the name of the coordinates in obsm", call.=FALSE)
    attrs <- list(`encoding-type`="anndata", `encoding-version`="0.1.0")
    sdattrs <- NULL
    if (!is.null(region_key)) {
        rk <- .region_key(obs[[region_key]])
        obs[[region_key]] <- rk$key
        sdattrs <- list(region=rk$regions, region_key=region_key, instance_key=instance_key)
        attrs <- c(attrs, list(`spatialdata-encoding-type`="ngff:regions_table", version="0.2"), sdattrs)
    } else {
        attrs <- c(attrs, list(`spatialdata-encoding-type`="ngff:regions_table", version="0.2"))
    }
    .zarr_group(tab, attrs)

    # obs: the row index is the column names of an object, the row names of a data.frame of images, else 0, 1, ...
    ids <- if (sce) colnames(x) else if (is.null(region_key)) rownames(obs)
    index <- if (is.null(ids)) \(a, b) .vlen_ids((a:b) - 1L) else \(a, b) ids[a:b]
    extra <- list()
    if (!is.null(region_key)) extra <- list(structure(\(a, b) (a:b) - 1L, name=instance_key,
        type=.zarr_int_type(0, max(n - 1, 256))))   # 16 bits at least: spatialdata refuses an 8-bit instance key
    .zarr_write_dataframe(file.path(tab, "obs"), obs, n, index, extra)

    # X, layers and var
    dict <- list(`encoding-type`="dict", `encoding-version`="0.1.0")
    for (g in c("layers", "obsm", "obsp", "uns", "varm", "varp")) .zarr_group(file.path(tab, g), dict)
    if (sce && length(SummarizedExperiment::assayNames(x)) < length(SummarizedExperiment::assays(x)))
        stop("every assay needs a name", call.=FALSE)
    an <- if (sce) SummarizedExperiment::assayNames(x) else character(0)
    if (length(an)) {
        .zarr_write_matrix(file.path(tab, "X"), SummarizedExperiment::assay(x, an[1], withDimnames=FALSE))
        for (a in an[-1])
            .zarr_write_matrix(file.path(tab, "layers", a), SummarizedExperiment::assay(x, a, withDimnames=FALSE))
    }
    var <- if (sce) as.data.frame(SummarizedExperiment::rowData(x), optional=TRUE) else data.frame()
    .zarr_check_columns(var, "rowData")
    vids <- if (sce) rownames(x) else NULL
    if (is.null(vids)) vids <- as.character(seq_len(nrow(var)) - 1L)
    .zarr_write_dataframe(file.path(tab, "var"), var, nrow(var), \(a, b) vids[a:b])

    # obsm: the coordinates and the reducedDims
    if (!is.null(xy)) .zarr_write_rows(file.path(tab, "obsm", "spatial"), n, "float64", spatial, ncol=2L)
    for (r in names(rds)) {
        v <- as.matrix(rds[[r]])
        .zarr_write_rows(file.path(tab, "obsm", r), n, if (is.integer(v)) "int32" else "float64",
            \(a, b) unname(v[a:b, , drop=FALSE]), ncol=ncol(v))
    }

    # uns: the SpatialData attributes
    if (!is.null(sdattrs)) {
        sp <- file.path(tab, "uns", "spatialdata_attrs")
        .zarr_group(sp, dict)
        reg <- sdattrs$region
        .zarr_write_rows(file.path(sp, "region"), length(reg), "string", \(a, b) reg[a:b],
            encoding=c("string-array", "0.2.0"))
        for (k in c("region_key", "instance_key")) .zarr_write_string(file.path(sp, k), sdattrs[[k]])
    }
    invisible(tab)
}

## An image: an OME-NGFF image element `name` under `path`/images, from a SpatialImage (imgData() of a
## SpatialExperiment), its pixels as uint8 RGB channels (c, y, x). `scale` takes image pixels to the coordinate
## system `cs`: a SpatialExperiment's coordinates are in full-resolution pixels and its image is scaled by its
## scaleFactor, so `scale` is 1 / scaleFactor.
.writeImage <- function(x, name, path, cs="global", scale=1) {
    ras <- as.matrix(SpatialExperiment::imgRaster(x))
    h <- nrow(ras); w <- ncol(ras)
    arr <- aperm(array(t(grDevices::col2rgb(as.vector(ras))), c(h, w, 3L)), c(3L, 1L, 2L))
    storage.mode(arr) <- "integer"
    dimnames(arr) <- list(c=NULL, y=NULL, x=NULL)
    p <- file.path(path, "images", name)
    axes <- list(list(name="c", type="channel"), list(name="y", type="space"), list(name="x", type="space"))
    sys <- \(nm) list(name=nm, axes=list(list(name="c", type="channel"), list(name="y", type="space", unit="unit"),
        list(name="x", type="space", unit="unit")))
    .zarr_group(p, list(
        ome=list(
            omero=list(channels=list(list(label="r"), list(label="g"), list(label="b"))),
            version="0.5-dev-spatialdata",
            multiscales=list(list(
                datasets=list(list(path="s0", coordinateTransformations=list(list(type="scale", scale=c(1, 1, 1))))),
                name=paste0("/images/", name),
                axes=axes,
                coordinateTransformations=list(list(type="scale", scale=c(1, scale, scale),
                    input=sys("cyx"), output=sys(cs)))))),
        spatialdata_attrs=list(version="0.3")))
    Rarr::write_zarr_array(arr, file.path(p, "s0"), chunk_dim=c(3L, min(h, 1024L), min(w, 1024L)),
        data_type="|u1", order="C", zarr_version=3L)
    invisible(p)
}

## The images of a SpatialExperiment (its imgData()), one image element per image, named <sample_id>_<image_id>,
## each in the coordinate system of its sample.
.writeImages <- function(x, path) {
    img <- SpatialExperiment::imgData(x)
    if (!nrow(img)) return(invisible(character(0)))
    .zarr_group(file.path(path, "images"))
    nms <- gsub("[^A-Za-z0-9_.-]", "_", paste0(img$sample_id, "_", img$image_id))
    if (anyDuplicated(nms)) stop("images with the same <sample_id>_<image_id>: ",
        paste(unique(nms[duplicated(nms)]), collapse=", "), call.=FALSE)
    for (i in seq_len(nrow(img)))
        .writeImage(img$data[[i]], nms[i], path, cs=as.character(img$sample_id[i]), scale=1 / img$scaleFactor[i])
    invisible(nms)
}

## A new, empty SpatialData store at `path` (replaced if it exists): its root group and tables group.
.create_store <- function(path) {
    unlink(path, recursive=TRUE)
    dir.create(file.path(path, "tables"), recursive=TRUE)
    path <- normalizePath(path)
    Rarr::write_zarr_group(path, "", zarr_version=3L)
    Rarr::write_zarr_group(path, "tables", zarr_version=3L)
    Rarr::write_zarr_attributes(path, new.zattrs=list(spatialdata_attrs=list(version="0.2")))
    path
}
