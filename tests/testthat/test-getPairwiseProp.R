## Tests for getPairwiseProp(): the fixed-k observed/expected proportion ratio
## and its use as a spicyR alternateResult.

make_cells <- function(seed = 1, n_subject = 8, img_per_subject = 2,
                       cells_per_img = 500) {
  set.seed(seed)
  blocks <- list()
  for (s in seq_len(n_subject)) {
    cond <- if (s <= n_subject / 2) "onset" else "long"
    for (j in seq_len(img_per_subject)) {
      nc <- cells_per_img
      ct <- sample(c("A", "B", "C"), nc, replace = TRUE, prob = c(0.3, 0.4, 0.3))
      x <- runif(nc, 0, 1000)
      y <- runif(nc, 0, 1000)
      if (cond == "onset") {
        # pull some A cells toward B cells to create a real condition effect
        bx <- mean(x[ct == "B"])
        by <- mean(y[ct == "B"])
        mv <- sample(which(ct == "A"), sum(ct == "A") %/% 2)
        x[mv] <- 0.5 * x[mv] + 0.5 * bx
        y[mv] <- 0.5 * y[mv] + 0.5 * by
      }
      blocks[[length(blocks) + 1]] <- data.frame(
        imageID = paste0("s", s, "_img", j),
        case = paste0("case", s),
        stage = cond,
        cellType = ct, x = x, y = y,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, blocks)
}

test_that("output is shaped and ordered for spicy(alternateResult=)", {
  cells <- make_cells()
  fmt <- spicyR:::.format_data(cells, "imageID", "cellType", c("x", "y"), FALSE)
  m <- getPairwiseProp(cells, k = 15)
  g <- suppressWarnings(suppressMessages(getPairwise(fmt, r = c(20, 50))))

  expect_true(is.matrix(m))
  # rows: one per image, in getImagePheno() order
  expect_identical(rownames(m), as.character(unique(fmt$imageID)))
  # columns: same set as getPairwise(), ordered as spicy() builds its labels
  expect_setequal(colnames(m), colnames(g))
  allTypes <- as.character(unique(fmt$cellType))
  expect_identical(
    colnames(m),
    paste(rep(allTypes, times = length(allTypes)),
      rep(allTypes, each = length(allTypes)),
      sep = "__"
    )
  )
})

test_that("all non-NA ratios are finite and non-negative", {
  cells <- make_cells()
  m <- getPairwiseProp(cells, k = 15)
  v <- m[!is.na(m)]
  expect_true(length(v) > 0)
  expect_true(all(is.finite(v) & v >= 0))
})

test_that("ratio is ~1 under a homogeneous random target thinning", {
  set.seed(42)
  reps <- 40
  vals <- vapply(seq_len(reps), function(i) {
    nc <- 800
    ct <- sample(c("A", "B"), nc, replace = TRUE, prob = c(0.5, 0.5))
    cells <- data.frame(
      imageID = "img1", cellType = ct,
      x = runif(nc, 0, 1000), y = runif(nc, 0, 1000)
    )
    getPairwiseProp(cells, k = 15, from = "A", to = "B")[1, 1]
  }, numeric(1))

  se <- sd(vals) / sqrt(reps)
  expect_lt(abs(mean(vals) - 1), 4 * se)
})

test_that("degenerate image-pairs return NA, or a floor when includeZeroCells", {
  set.seed(7)
  nc <- 400
  # image with no "C" cells at all
  ct <- sample(c("A", "B"), nc, replace = TRUE)
  cells <- data.frame(
    imageID = "img1", cellType = ct,
    x = runif(nc, 0, 1000), y = runif(nc, 0, 1000)
  )

  m_na <- getPairwiseProp(cells, k = 15, from = c("A", "B", "C"), to = c("A", "B", "C"))
  expect_true(is.na(m_na[1, "A__C"]))
  expect_true(is.na(m_na[1, "C__B"]))
  expect_false(is.na(m_na[1, "A__B"]))

  m_zero <- getPairwiseProp(cells,
    k = 15, from = c("A", "B", "C"), to = c("A", "B", "C"),
    includeZeroCells = TRUE
  )
  # target C absent -> p0 = 0 -> ratio undefined, stays NA even with includeZeroCells
  expect_true(is.na(m_zero[1, "A__C"]))
  # reference C absent but target A present -> floored at pi_hat = 0, i.e. v = 0
  expect_equal(unname(m_zero[1, "C__A"]), 0)
  # A and B both present -> real ratio
  expect_false(is.na(m_zero[1, "A__B"]))
})

test_that("images with <= k cells return all-NA rows", {
  set.seed(3)
  small <- data.frame(
    imageID = "tiny", cellType = rep(c("A", "B"), 5),
    x = runif(10), y = runif(10)
  )
  big <- data.frame(
    imageID = "big",
    cellType = sample(c("A", "B"), 400, replace = TRUE),
    x = runif(400), y = runif(400)
  )
  cells <- rbind(small, big)

  m <- getPairwiseProp(cells, k = 15)
  expect_true(all(is.na(m["tiny", ])))
  expect_false(all(is.na(m["big", ])))
})

test_that("feeds spicy() as alternateResult (weights = FALSE)", {
  cells <- make_cells()
  m <- getPairwiseProp(cells, k = 15)

  finite_condition_p <- function(res) {
    pv <- as.data.frame(res$p.value)
    cond <- pv[, grep("condition", colnames(pv)), drop = FALSE]
    any(vapply(cond, function(x) any(is.finite(x)), logical(1)))
  }

  res_mem <- suppressWarnings(suppressMessages(
    spicy(cells,
      condition = "stage", subject = "case",
      alternateResult = m, weights = FALSE
    )
  ))
  res_lm <- suppressWarnings(suppressMessages(
    spicy(cells, condition = "stage", alternateResult = m, weights = FALSE)
  ))

  expect_s4_class(res_mem, "SpicyResults")
  expect_true(finite_condition_p(res_mem))
  expect_true(finite_condition_p(res_lm))
})

test_that("weights = TRUE works end-to-end on diabetesData", {
  skip_if_not_installed("lme4")
  data("diabetesData")
  m <- getPairwiseProp(diabetesData, k = 15)

  res <- suppressWarnings(suppressMessages(
    spicy(diabetesData,
      condition = "stage", subject = "case",
      alternateResult = m, weights = TRUE
    )
  ))

  pv <- as.matrix(res$p.value)
  expect_true("conditionOnset" %in% colnames(pv))
  # the weight model does not degenerate: the vast majority of pairs fit
  # (a handful are constant-0 ratios that no linear model can test)
  expect_gt(mean(is.finite(pv)), 0.95)
})
