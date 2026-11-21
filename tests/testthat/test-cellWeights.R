## cellWeights: 0/1 weights give the results of the cell types, in each method (stage has three levels, so
## the tables have one row per pair and level)
data("diabetesData")
cells <- spicyR:::.format_data(diabetesData, "imageID", "cellType", c("x", "y"), FALSE)
lev <- levels(droplevels(cells$cellType))
onehot <- matrix(0, nrow(cells), length(lev), dimnames = list(NULL, lev))
onehot[cbind(seq_len(nrow(cells)), match(as.character(cells$cellType), lev))] <- 1

test_that("getPairwise() with 0/1 weights equals the cell types", {
  cf <- droplevels(cells)
  a <- getPairwise(cf, r = c(20, 50), window = "square")
  b <- getPairwise(cf, r = c(20, 50), window = "square", cellWeights = onehot)
  expect_identical(unname(a), unname(b[rownames(a), colnames(a)]))
})

test_that("spicy() cell method with 0/1 weights equals the cell types", {
  for (eff in c("allocation", "count")) {
    a <- spicy(cells, condition = "stage", subject = "case", r = 50, from = c("Tc", "Th"), to = c("beta", "alpha", "Th"), effect = eff)$cellResults
    b <- spicy(cells, condition = "stage", subject = "case", r = 50, from = c("Tc", "Th"), to = c("beta", "alpha", "Th"), effect = eff,
               cellWeights = onehot)$cellResults
    a <- a[a$from != a$to, ]; b <- b[match(paste(a$from, a$to, a$level), paste(b$from, b$to, b$level)), ]
    expect_equal(b$excess_difference, a$excess_difference, tolerance = 1e-10)
    expect_equal(b$p_value, a$p_value, tolerance = 1e-10)
  }
})

test_that("kontextualEngine() with 0/1 weights equals the cell types", {
  pdf <- data.frame(from = c("alpha", "beta"), to = "Tc", parent_name = "tcells")
  pdf$parent <- list(c("Tc", "Th"), c("Tc", "Th"))
  a <- kontextualEngine(cells, pdf, condition = "stage", subject = "case", r = 50)$cellResults
  b <- kontextualEngine(cells, pdf, condition = "stage", subject = "case", r = 50, cellWeights = onehot)$cellResults
  b <- b[match(paste(a$from, a$to, a$level), paste(b$from, b$to, b$level)), ]
  expect_equal(b$excess_difference, a$excess_difference, tolerance = 1e-10)
  expect_equal(b$p_value, a$p_value, tolerance = 1e-10)
})

test_that("soft weights are checked", {
  expect_error(spicy(cells, condition = "stage", subject = "case", r = 50, cellWeights = onehot[-1, ]), "rows")
  expect_error(spicy(cells, condition = "stage", subject = "case", r = 50, cellWeights = 2 * onehot), "between 0 and 1")
  expect_error(spicy(cells, condition = "stage", subject = "case", r = 50, from = "nope", to = "Tc", cellWeights = onehot),
               "not a column")
})

test_that("a vector of weights equals the matrix with each cell's weight on its own type", {
  set.seed(1); q <- runif(nrow(cells), 0.5, 1)
  m <- onehot * q
  a <- getPairwise(cells, r = 50, window = "square", cellWeights = q)
  b <- getPairwise(cells, r = 50, window = "square", cellWeights = m)
  expect_equal(a, b[, colnames(a)])
  a <- spicy(cells, condition = "stage", subject = "case", r = 50, from = "Tc", to = c("beta", "Th"), cellWeights = q)$cellResults
  b <- spicy(cells, condition = "stage", subject = "case", r = 50, from = "Tc", to = c("beta", "Th"), cellWeights = m)$cellResults
  expect_equal(a$p_value, b$p_value)
  expect_error(spicy(cells, condition = "stage", subject = "case", r = 50, cellWeights = q[-1]), "values but there are")
})

## Labels with 0/1 weights that share cells (a label with itself, or a group and one of its subtypes), allocation:
## random labelling of all the image's cells.

## All distinct arrangements of a multiset of category codes (one row each).
arrangements <- function(cats) {
  if (length(cats) <= 1L) return(matrix(cats, nrow = 1L))
  u <- unique(cats)
  do.call(rbind, lapply(u, function(v) cbind(v, arrangements(cats[-match(v, cats)]))))
}
## O = the number of TO cells with a FROM cell other than themselves within r, for 0/1 vectors Fs and Ts.
overlap_count <- function(nb, Fs, Ts) sum(vapply(seq_along(nb), function(i) Ts[i] && any(Fs[nb[[i]]]), TRUE))
neighbours <- function(x, y, r) lapply(seq_along(x), function(i) setdiff(which((x - x[i])^2 + (y - y[i])^2 <= r^2), i))

test_that("the random-labelling mean and variance of labels that share cells are exact", {
  ns <- asNamespace("spicyR")
  set.seed(1)
  check <- function(N, r, nF, nT, nFT, self = nF == nT && nFT == nT) {
    x <- runif(N); y <- runif(N)
    d <- ns$dataset_create(x, y, rep(0L, N), c(0L, N), 1L)
    ns$dataset_build_radius_index(d, r)
    nb <- neighbours(x, y, r)
    # 1 = FROM and TO, 2 = FROM only, 3 = TO only, 4 = neither
    A <- arrangements(c(rep(1, nFT), rep(2, nF - nFT), rep(3, nT - nFT), rep(4, N - nF - nT + nFT)))
    O <- apply(A, 1, function(lab) overlap_count(nb, lab <= 2, lab %in% c(1, 3)))
    m <- ns$dataset_overlap_moments(d, matrix(nF), matrix(nT), matrix(nFT), self, 0L, 50L)
    c(m$E[1] - mean(O), m$V[1] - (mean(O^2) - mean(O)^2))
  }
  # self-pairs (nF = nT = nFT), a group whose cells in an image are all of one subtype (the same counts, not a
  # self-pair), a group and its subtype both ways, partial and no overlap, a single FROM cell, every cell labelled
  d <- rbind(check(9, 0.4, 3, 3, 3), check(10, 0.6, 4, 4, 4), check(9, 0.4, 3, 3, 3, self = FALSE),
             check(9, 0.5, 2, 5, 2), check(9, 0.5, 5, 2, 2), check(10, 0.35, 3, 6, 3), check(9, 0.5, 3, 4, 1),
             check(9, 0.5, 2, 3, 0), check(10, 0.5, 1, 5, 1), check(8, 0.5, 8, 3, 3), check(8, 0.5, 3, 8, 3))
  expect_lt(max(abs(d)), 1e-9)
})

test_that("a group and its subtypes: each image's excess counts every cell, around the exact mean", {
  set.seed(2)
  imgs <- lapply(seq_len(12), function(i) {
    z <- data.frame(x = runif(9, 0, 100), y = runif(9, 0, 100), cellType = sample(c("a", "b", "c", "c"), 9, TRUE))
    z$imageID <- sprintf("i%02d", i); z$patient <- sprintf("p%02d", (i + 1) %/% 2)
    z$condition <- if (i <= 6) "x" else "y"
    z })
  d <- do.call(rbind, imgs)
  M <- outer(d$cellType, c("a", "b", "c"), "==") * 1
  colnames(M) <- c("a", "b", "c")
  M <- cbind(M, G = M[, "a"] + M[, "b"])
  r <- 40
  res <- spicy(d, condition = "condition", subject = "patient", r = r, cellWeights = M, from = c("G", "a"),
               to = c("G", "a"))
  for (p in list(c("G", "a"), c("a", "G"), c("G", "G"), c("a", "a"))) {
    rows <- lapply(imgs, function(z) {
      Mi <- M[d$imageID == z$imageID[1], ]
      nb <- neighbours(z$x, z$y, r)
      A <- arrangements(z$cellType)   # random labelling of the image's cells
      lab <- function(L, l) if (l == "G") L %in% c("a", "b") else L == l
      O0 <- apply(A, 1, function(L) overlap_count(nb, lab(L, p[1]), lab(L, p[2])))
      c(O = overlap_count(nb, Mi[, p[1]] == 1, Mi[, p[2]] == 1), E = mean(O0), V = mean(O0^2) - mean(O0)^2,
        nT = sum(Mi[, p[2]]))
    })
    rows <- as.data.frame(do.call(rbind, rows))
    keep <- rows$nT > 0 & rows$V > 1e-9
    avoid <- sum(rows$O[keep] - rows$E[keep]) < 0
    n <- if (avoid) rows$E else rows$nT - rows$E
    keep <- keep & n >= 1
    y <- ifelse(keep, (rows$O - rows$E) / n, NA_real_)
    got <- res$pairwiseAssoc[[paste0(p[1], "__", p[2])]]
    expect_equal(got, y, tolerance = 1e-9, info = paste(p, collapse = " -> "))
  }
})

test_that("with 0/1 weights a label with itself gives the cell type's self-pair; other pairs are unchanged", {
  sc <- sim_cells(seed = 5, n_pat = 5)
  lev <- sort(unique(sc$cellType))
  M <- outer(sc$cellType, lev, "==") * 1
  colnames(M) <- lev
  Mg <- cbind(M, lymphocyte = M[, "T"] + M[, "B"])
  a <- spicy(sc, condition = "condition", subject = "patient", r = 30, from = c("tumour", "T"),
             to = c("tumour", "T", "B"))$cellResults
  for (W in list(M, Mg)) {
    b <- spicy(sc, condition = "condition", subject = "patient", r = 30, from = c("tumour", "T"),
               to = c("tumour", "T", "B"), cellWeights = W)$cellResults
    expect_identical(paste(b$from, b$to), paste(a$from, a$to))
    expect_true(any(a$from == a$to & !is.na(a$p_value)))
    for (col in c("excess_difference", "se", "df", "p_value", "tau2"))
      expect_equal(b[[col]], a[[col]], tolerance = 1e-10, info = col)
  }
  # survival: the same self-pairs
  sc$os <- survival::Surv(stats::setNames(seq_along(unique(sc$patient)) + 0.5, unique(sc$patient))[sc$patient],
                          rep(c(1, 0, 1), length.out = length(unique(sc$patient)))[match(sc$patient, unique(sc$patient))])
  a <- spicy(sc, condition = "os", subject = "patient", r = 30, from = "T", to = c("T", "B"))$cellResults
  b <- spicy(sc, condition = "os", subject = "patient", r = 30, from = "T", to = c("T", "B"), cellWeights = M)$cellResults
  expect_identical(paste(b$from, b$to), paste(a$from, a$to))
  expect_equal(b$p_value, a$p_value, tolerance = 1e-10)
  expect_equal(b$score_coefficient, a$score_coefficient, tolerance = 1e-10)
})

test_that("soft weights and the count effect: the moments under permutation of each cell's pair of weights are exact", {
  ns <- asNamespace("spicyR")
  perms <- function(n) { if (n == 1L) return(matrix(1L)); p <- perms(n - 1L)
    do.call(rbind, lapply(seq_len(n), function(i) cbind(i, p + (p >= i)))) }
  check <- function(N, r, a, b, seed) {
    set.seed(seed); x <- runif(N); y <- runif(N)
    d <- ns$dataset_create(x, y, rep(0L, N), c(0L, N), 1L)
    ns$dataset_build_radius_index(d, r)
    nb <- neighbours(x, y, r)
    P <- perms(N)
    out <- numeric(0)
    for (alloc in c(TRUE, FALSE)) {
      stat <- function(pi) { A <- a[pi]; B <- b[pi]
        h <- vapply(seq_len(N), function(l) if (alloc) 1 - prod(1 - A[nb[[l]]]) else sum(A[nb[[l]]]), 0)
        c(sum(B * h), sum(h)) }
      S <- apply(P, 1, stat)
      m <- ns$dataset_overlap_soft_rows(d, a, matrix(b), alloc)
      obs <- stat(seq_len(N))
      if (alloc) out <- c(out, m[1] - obs[1], m[2] - mean(S[1, ]), m[3] - mean((S[1, ] - mean(S[1, ]))^2))
      else {   # the count effect: O - p L, with E = p L
        p <- m[2] / obs[2]; D <- S[1, ] - p * S[2, ]
        out <- c(out, m[1] - obs[1], mean(D), m[3] - mean((D - mean(D))^2))
      }
    }
    out
  }
  set.seed(10)
  w <- runif(7)
  d <- rbind(check(7, 0.5, runif(7), runif(7), 1),                             # soft, partly shared
             check(7, 0.6, w, w, 2),                                           # a soft label with itself
             check(7, 0.5, c(1, 1, .7, .2, 0, 0, .4), c(.5, 1, 0, .3, 0, .8, .4), 3),   # some `from` weights 1
             check(7, 0.45, c(1, 1, 1, 1, 0, 0, 0), c(1, 1, 0, 0, 0, 0, 0), 4),        # 0/1 group -> subtype
             check(7, 0.9, runif(7), runif(7), 6))                             # every cell a neighbour
  expect_lt(max(abs(d)), 1e-12)
})

test_that("with 0/1 weights the general moments are those of the counts", {
  ns <- asNamespace("spicyR")
  set.seed(3); N <- 400; x <- runif(N, 0, 200); y <- runif(N, 0, 200)
  d <- ns$dataset_create(x, y, rep(0L, N), c(0L, N), 1L)
  ns$dataset_build_radius_index(d, 20)
  F <- rbinom(N, 1, 0.3); T <- F * rbinom(N, 1, 0.5)
  for (p in list(list(F, T, FALSE), list(T, F, FALSE), list(F, F, TRUE))) {
    s <- ns$dataset_overlap_soft_rows(d, p[[1]], matrix(p[[2]]), TRUE)
    m <- ns$dataset_overlap_moments(d, matrix(sum(p[[1]])), matrix(sum(p[[2]])), matrix(sum(p[[1]] * p[[2]])), p[[3]],
                                    0L, 50L)
    expect_equal(s[2:3], c(m$E, m$V), tolerance = 1e-9)
  }
  # several pairs in one call, with images of 0, 1 and 400 cells: as one pair at a time
  d2 <- ns$dataset_create(c(x, 5, x[1:3]), c(y, 5, y[1:3]), rep(0L, N + 4), c(0L, N, N, N + 1L, N + 4L), 1L)
  ns$dataset_build_radius_index(d2, 20)
  a <- c(F * runif(N), 1, 0.5, 1, 0); B <- cbind(c(T, 1, 1, 0, 1), c(runif(N), 0.3, 0.2, 0.2, 0.9))
  for (alloc in c(TRUE, FALSE)) {
    both <- ns$dataset_overlap_soft_rows(d2, a, B, alloc)
    one <- lapply(1:2, function(q) ns$dataset_overlap_soft_rows(d2, a, B[, q, drop = FALSE], alloc))
    expect_identical(both[, c(1, 3, 5, 7)], one[[1]])
    expect_identical(both[, c(2, 4, 6, 8)], one[[2]])
  }
})

## Soft weights of sim_cells(seed = 5, n_pat = 5): a dense matrix of probabilities, a confidence for each cell's own
## type, and the dense matrix with the B cells' rows set to 0/1 (as tests/testthat/reference_1.99.59_soft.rds).
soft_inputs <- function() {
  sc <- sim_cells(seed = 5, n_pat = 5)
  lev <- sort(unique(sc$cellType), method = "radix")   # the same order in every locale
  onehot <- outer(sc$cellType, lev, "==") * 1; colnames(onehot) <- lev
  set.seed(8)
  dense <- onehot * 3 + matrix(runif(length(onehot)), nrow(onehot)); dense <- dense / rowSums(dense)
  q <- runif(nrow(sc), 0.5, 1)
  mixed <- dense; mixed[sc$cellType == "B", ] <- onehot[sc$cellType == "B", ]
  list(cells = sc, dense = dense, q = q, mixed = mixed)
}

test_that("soft pairs whose `to` cells all are candidates give the results of spicyR 1.99.59", {
  # reference_1.99.59_soft.rds: spicy() of spicyR 1.99.59 with these weights (no self-pairs; every `to` cell has
  # `from` weight below 1)
  ref <- readRDS(test_path("reference_1.99.59_soft.rds"))
  s <- soft_inputs()
  for (k in names(ref)) {
    nm <- strsplit(k, " ")[[1]]
    now <- spicy(s$cells, condition = "condition", subject = "patient", r = 30, from = c("tumour", "T", "B"),
                 to = c("T", "B", "macro"), cellWeights = s[[nm[1]]], effect = nm[2])$cellResults
    # p_adj is over every tested pair, and the self-pairs are now among them
    cols <- setdiff(names(ref[[k]]), "p_adj")
    expect_equal(now[now$from != now$to, cols], ref[[k]][, cols], tolerance = 1e-12, label = k)
  }
})

test_that("soft self-pairs, a soft subtype of a 0/1 group, and the count effect are tested", {
  s <- soft_inputs()
  W <- cbind(s$dense, lymphocyte = (s$cells$cellType %in% c("T", "B")) * 1)
  for (eff in c("allocation", "count")) {
    b <- spicy(s$cells, condition = "condition", subject = "patient", r = 30, from = c("T", "lymphocyte"),
               to = c("T", "lymphocyte"), cellWeights = W, effect = eff)$cellResults
    expect_identical(paste(b$from, b$to), c("T T", "T lymphocyte", "lymphocyte T", "lymphocyte lymphocyte"))
    expect_true(all(is.finite(b$p_value)), label = eff)
  }
})

test_that("count effect: with 0/1 weights a label with itself gives the cell type's self-pair", {
  sc <- sim_cells(seed = 5, n_pat = 5)
  lev <- sort(unique(sc$cellType))
  M <- outer(sc$cellType, lev, "==") * 1
  colnames(M) <- lev
  a <- spicy(sc, condition = "condition", subject = "patient", r = 30, from = c("tumour", "T"), to = c("tumour", "T", "B"),
             effect = "count")$cellResults
  b <- spicy(sc, condition = "condition", subject = "patient", r = 30, from = c("tumour", "T"), to = c("tumour", "T", "B"),
             effect = "count", cellWeights = M)$cellResults
  expect_identical(paste(b$from, b$to), paste(a$from, a$to))
  for (col in c("excess_difference", "se", "df", "p_value", "tau2"))
    expect_equal(b[[col]], a[[col]], tolerance = 1e-10, info = col)
})
