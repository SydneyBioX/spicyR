## Kontextual designs (`parent`). The Poisson design's per-image estimate must
## equal Statial's Kontextual statistic, here recomputed by brute force from its
## definition (Statial's calcKontextual.R, .Kontext); the binomial design must
## reduce to the ordinary one when the context is every cell type.

make_context_cells <- function(seed = 5) {
  set.seed(seed)
  do.call(rbind, lapply(1:8, function(i) {
    n <- rpois(4, c(150, 120, 90, 60)); tt <- rep(c("Tum", "CD8", "CD4", "Mac"), n)
    left <- runif(length(tt)) < ifelse(tt == "Tum", 0.3, 0.7)
    data.frame(x = ifelse(left, runif(length(tt), 0, 500), runif(length(tt), 500, 1000)),
               y = runif(length(tt), 0, 1000), cellType = tt, imageID = sprintf("i%d", i),
               response = if (i <= 4) "alpha" else "beta", stringsAsFactors = FALSE)
  }))
}

immune <- c("CD8", "CD4", "Mac")

## Statial's statistic for one image, without edge correction: every count
## includes the cell itself (closepairs(distinct = FALSE)), lambda_c(x) is the
## context count within r over pi r^2, and K / (pi r^2) is returned.
brute_kontextual <- function(img, from, to, parent, r) {
  d <- as.matrix(stats::dist(img[, c("x", "y")]))
  close <- d <= r
  lam <- rowSums(close[, img$cellType %in% parent, drop = FALSE]) / (pi * r^2)
  a <- which(img$cellType == from); b <- which(img$cellType == to)
  share <- length(b) / sum(img$cellType %in% parent)
  num <- sum(close[a, b] * outer(lam[a], lam[b] * share, "/"))
  num / sum(lam[a]) / (pi * r^2)
}

test_that("the Poisson design's per-image estimate is Statial's Kontextual statistic", {
  cells <- make_context_cells()
  ord <- order(cells$imageID); cells <- cells[ord, ]
  labs <- sort(unique(cells$imageID)); img <- match(cells$imageID, labs) - 1L
  tl <- unique(cells$cellType); tc <- match(cells$cellType, tl) - 1L
  off <- as.integer(c(0, cumsum(tabulate(img + 1L, nbins = length(labs)))))
  d <- dataset_create(cells$x, cells$y, tc, off, length(tl))
  dataset_build_radius_index(d, 60)
  dataset_build_context(d, as.integer(match(immune, tl) - 1L), "convex", FALSE)
  for (pr in list(c("Tum", "CD8"), c("CD4", "Mac"))) {
    md <- dataset_kontextual_model_data(d, match(pr[1], tl) - 1L, match(pr[2], tl) - 1L)
    mine <- as.numeric(tapply(md$n, md$image, sum) / tapply(md$density, md$image, sum))
    ref <- vapply(labs, function(l) brute_kontextual(cells[cells$imageID == l, ], pr[1], pr[2], immune, 60),
                  numeric(1))
    expect_equal(mine, unname(ref), tolerance = 1e-12)
  }
})

test_that("edge correction changes only the cells whose disc leaves the window", {
  cells <- make_context_cells()
  cells <- cells[order(cells$imageID), ]
  one <- cells[cells$imageID == "i1", ]
  tl <- unique(one$cellType)
  d <- dataset_create(one$x, one$y, match(one$cellType, tl) - 1L, c(0L, nrow(one)), length(tl))
  dataset_build_radius_index(d, 60)
  from <- match("Tum", tl) - 1L; to <- match("CD8", tl) - 1L
  dataset_build_context(d, as.integer(match(immune, tl) - 1L), "rectangle", FALSE)
  plain <- dataset_kontextual_model_data(d, from, to)
  dataset_build_context(d, as.integer(match(immune, tl) - 1L), "rectangle", TRUE)
  corrected <- dataset_kontextual_model_data(d, from, to)
  expect_identical(plain$row, corrected$row)
  x <- one$x[plain$row + 1L]; y <- one$y[plain$row + 1L]
  inside <- pmin(x - min(one$x), max(one$x) - x, y - min(one$y), max(one$y) - y) > 60
  expect_equal(corrected$density[inside], plain$density[inside])
  expect_true(all(corrected$weight[!inside] >= plain$weight[!inside]))
})

test_that("Kontextual pairs are ordered, with `to` in the context", {
  out <- spicy_glm(make_context_cells(), condition = "response", r = 60, parent = immune)
  fitted <- rbind(out$results[, c("from", "to")], out$skipped[, c("from", "to")])
  expect_equal(nrow(fitted), 4 * 3)
  expect_true(all(fitted$to %in% immune))
  expect_identical(out$parent, immune)
  expect_error(spicy_glm(make_context_cells(), condition = "response", r = 60, parent = immune,
                         from = "CD8", to = "Tum"), "inside `parent`")
  expect_error(spicy_glm(make_context_cells(), condition = "response", r = 60, parent = immune,
                         sigma = 200), "alternative")
})

test_that("a binomial context of every cell type is the ordinary binomial design", {
  cells <- make_context_cells()
  all_types <- unique(cells$cellType)
  ctx <- spicy_glm(cells, condition = "response", family = "binomial", k = 10, parent = all_types,
                   from = "Tum", to = "CD8")$results
  plain <- spicy_glm(cells, condition = "response", family = "binomial", k = 10,
                     from = "Tum", to = "CD8")$results
  expect_equal(ctx$log_odds_ratio, plain$log_odds_ratio, tolerance = 1e-12)
  expect_equal(ctx$p_value, plain$p_value, tolerance = 1e-10)
})

test_that("the binomial design's MLE with per-cell trials is glm's", {
  cells <- make_context_cells()
  ord <- order(cells$imageID); cells <- cells[ord, ]
  labs <- sort(unique(cells$imageID)); img <- match(cells$imageID, labs) - 1L
  tl <- unique(cells$cellType); tc <- match(cells$cellType, tl) - 1L
  off <- as.integer(c(0, cumsum(tabulate(img + 1L, nbins = length(labs)))))
  d <- dataset_create(cells$x, cells$y, tc, off, length(tl))
  dataset_build_knn(d, 10L)
  dataset_build_context(d, as.integer(match(immune, tl) - 1L), "convex", FALSE)
  md <- dataset_kontextual_binomial_model_data(d, match("Tum", tl) - 1L, match("CD8", tl) - 1L)
  group <- as.integer(md$image >= 4L)
  expect_true(all(md$n <= md$trials) && all(md$trials > 0))
  fit <- fit_pair_binomial_trials_cpp(md$image, md$image, group, md$n, md$trials, md$p0, "mle", "fast")
  ref <- stats::glm(cbind(md$n, md$trials - md$n) ~ 0 + factor(group), offset = stats::qlogis(md$p0),
                    family = stats::binomial())
  expect_equal(fit$beta, unname(stats::coef(ref)), tolerance = 1e-8)
})
