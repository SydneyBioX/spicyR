## A SpatialData-style Zarr store of a data.frame of cells, written with Rarr: tables/<table> is an AnnData group whose
## obs holds the columns (factors as categoricals, character as string arrays, numbers as arrays) and whose obsm
## holds the coordinates ("spatial", unless `coords_in_obs`). `chunk` rows per chunk.
write_cell_store <- function(d, path, table = "table", region_key = "imageID", xy = c("x", "y"), coords_in_obs = FALSE,
                             chunk = 1000L, zarr_version = 3L) {
  w_attrs <- function(p, a) Rarr::write_zarr_attributes(p, new.zattrs = a, zarr_version = zarr_version)
  w_group <- function(p)   # p under `path`
    Rarr::write_zarr_group(path, if (p == path) "" else substring(p, nchar(path) + 2L), zarr_version = zarr_version)
  w_array <- function(x, p, ch = chunk) {
    x <- as.array(x); dims <- dim(x)
    Rarr::write_zarr_array(x, p, chunk_dim = pmin(c(ch, dims[-1]), dims), zarr_version = zarr_version,
                           nchar = if (is.character(x)) max(1L, nchar(x, "bytes"), na.rm = TRUE))
  }
  unlink(path, recursive = TRUE); dir.create(path); w_group(path)
  w_attrs(path, list(spatialdata_attrs = list(version = "0.2")))
  tab <- file.path(path, "tables", table)
  w_group(tab)
  w_attrs(tab, list(`encoding-type` = "anndata", `encoding-version` = "0.1.0", region = as.list(unique(as.character(d[[region_key]]))),
                    region_key = region_key, instance_key = "instance_id"))
  obs <- file.path(tab, "obs"); w_group(obs)
  cols <- setdiff(names(d), if (!coords_in_obs) xy)
  w_attrs(obs, list(`column-order` = as.list(cols), `_index` = "_index", `encoding-type` = "dataframe",
                    `encoding-version` = "0.2.0"))
  w_array(as.character(seq_len(nrow(d)) - 1L), file.path(obs, "_index"))
  for (cn in cols) {
    v <- d[[cn]]; p <- file.path(obs, cn)
    if (is.factor(v)) {
      w_group(p)
      w_attrs(p, list(ordered = is.ordered(v), `encoding-type` = "categorical", `encoding-version` = "0.2.0"))
      codes <- as.integer(v) - 1L; codes[is.na(codes)] <- -1L
      w_array(codes, file.path(p, "codes")); w_array(levels(v), file.path(p, "categories"))
    } else {
      w_array(v, p)
    }
  }
  w_group(file.path(tab, "obsm"))
  if (!coords_in_obs) w_array(cbind(as.numeric(d[[xy[1]]]), as.numeric(d[[xy[2]]])), file.path(tab, "obsm", "spatial"))
  invisible(path)
}
