# The ratio as the original R implementation computed it, one image at a time.
propReference <- function(df, k, from, to, includeZeroCells) {
  types <- as.character(df$cellType)
  labels <- paste(rep(from, times = length(to)), rep(to, each = length(from)), sep = "__")
  out <- t(vapply(split(df, factor(df$imageID, unique(df$imageID))), function(z) {
    ty <- as.character(z$cellType)
    v <- rep(NA_real_, length(labels))
    if (nrow(z) <= k) return(v)
    pp <- spatstat.geom::ppp(z$x, z$y, range(z$x), range(z$y), check = FALSE)
    nnType <- matrix(ty[spatstat.geom::nnwhich(pp, k = seq_len(k))], ncol = k)
    for (i in seq_along(labels)) {
      A <- rep(from, times = length(to))[i]; B <- rep(to, each = length(from))[i]
      nA <- sum(ty == A); nB <- sum(ty == B)
      if (nA == 0 || nB == 0 || nB == nrow(z)) {
        if (includeZeroCells && nB > 0) v[i] <- 0
        next
      }
      v[i] <- mean(rowSums(nnType[ty == A, , drop = FALSE] == B)) / k / (nB / nrow(z))
    }
    v
  }, numeric(length(labels))))
  colnames(out) <- labels
  out
}

test_that("getPairwiseProp() matches the per-image R computation", {
  cd <- as.data.frame(SummarizedExperiment::colData(diabetesData))
  cd <- cd[cd$imageID %in% unique(cd$imageID)[1:6], ]
  df <- data.frame(imageID = cd$imageID, cellType = cd$cellType, x = cd$x, y = cd$y)
  allTypes <- unique(as.character(df$cellType))

  for (args in list(
    list(k = 15, from = allTypes, to = allTypes, includeZeroCells = FALSE),
    list(k = 5, from = c("Tc", "beta"), to = c("Th", "alpha", "beta"), includeZeroCells = TRUE)
  )) {
    res <- getPairwiseProp(df, k = args$k, from = args$from, to = args$to,
                           includeZeroCells = args$includeZeroCells, cores = 2)
    ref <- propReference(df, args$k, args$from, args$to, args$includeZeroCells)
    expect_equal(unname(res), unname(ref), tolerance = 1e-10)
    expect_identical(colnames(res), colnames(ref))
  }
})

test_that("an image with at most k cells is all NA", {
  df <- data.frame(imageID = rep(c("a", "b"), c(5, 40)),
                   cellType = rep(c("u", "v"), length.out = 45),
                   x = c(1:5, runif(40)), y = c(5:1, runif(40)))
  res <- getPairwiseProp(df, k = 5)
  expect_true(all(is.na(res["a", ])))
  expect_false(anyNA(res["b", ]))
})
