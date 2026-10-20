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
