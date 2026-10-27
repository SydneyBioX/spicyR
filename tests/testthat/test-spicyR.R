test_that(
    "the image method is unchanged from spicyR 1.x",
    {
        original_result <- readRDS(
            system.file("testdata/original_result.rds", package = "spicyR")
        )
        expect_equal(
            suppressWarnings(spicy(diabetesData,
                condition = "stage", subject = "case",
                from = "Tc", to = "Th", method = "image"
            )),
            original_result,
            tolerance = 0.01
        )
    }
)

test_that("the image method recycles a single from over several to", {
  res <- suppressWarnings(suppressMessages(spicy(diabetesData, condition = "stage", subject = "case",
                                                 from = "Tc", to = c("Th", "beta"), method = "image")))
  expect_equal(sort(rownames(res$p.value)), c("Tc__Th", "Tc__beta"))
})

test_that("the image method does not test a pair with the same value in every image", {
    # A fit to a constant response gives a p-value from rounding error alone, so the pair's p-value is NA (and it
    # stays out of the multiple-testing correction).
    set.seed(1)
    pheno <- data.frame(stage = rep(c("a", "b"), each = 12), case = rep(seq_len(8), each = 3))
    w <- runif(24, 0.5, 2)
    mem <- function(y) spatialMEM(y, cells = NULL, subject = "case", condition = "stage", covariates = NULL,
                                  weightFunction = w, cellCounts = NULL, pheno = pheno)
    expect_identical(mem(rep(-25, 24)), NA)
    expect_identical(mem(c(NA, rep(0, 23))), NA)
    expect_true(is.matrix(mem(rnorm(24))))
    lmFit <- function(y) spatialLM(y, from = 1, to = 1, cells = NULL, condition = "stage", covariates = NULL,
                                   weightFunction = w, cellCounts = matrix(1, 24, 1), pheno = pheno)
    expect_identical(lmFit(rep(-25, 24)), NA)
    expect_s3_class(lmFit(rnorm(24)), "lm")
    p <- cleanMEM(list(a__b = mem(rep(-25, 24)), b__a = mem(rnorm(24))))$p.value
    expect_true(is.na(p["a__b", "conditionb"]) && !is.na(p["b__a", "conditionb"]))
})

test_that("the weight model follows scam()'s smoothing-parameter search", {
    # The GCV surface is flat near its minimum, so a search that strays from
    # scam()'s path (e.g. with a finite-difference gradient) stops elsewhere
    # within tolerance and moves the weights; the fit must match scam() closely.
    skip_if_not_installed("scam")
    res <- suppressWarnings(suppressMessages(spicy(diabetesData,
        condition = "stage", subject = "case", method = "image"
    )))
    pa <- do.call(cbind, res$pairwiseAssoc)
    rS <- as.vector(apply(pa, 2, function(x) {
        if (sum(!is.na(x)) > 1 && stats::sd(x, na.rm = TRUE) > 0) (x - mean(x, na.rm = TRUE))^2 else rep(NA, length(x))
    }))
    ok <- !is.na(rS)
    y <- log10(rS[ok] + 1)
    x1 <- log10(as.vector(res$nCells[, res$comparisons$from])[ok] + 1)
    x2 <- log10(as.vector(res$nCells[, res$comparisons$to])[ok] + 1)
    z <- mpdWeightFit(y, x1, x2, x1, x2)
    fit <- scam::scam(y ~ s(x1, bs = "mpd") + s(x2, bs = "mpd"))
    expect_equal(unname(z), unname(stats::fitted(fit)), tolerance = 1e-10)
})

test_that("the weight model fits a cell type paired with itself, as scam() does", {
    set.seed(3)
    x <- log10(rpois(60, 200) + 1)
    y <- 2 - 0.8 * x + rnorm(60, sd = 0.1)
    xnew <- seq(min(x), max(x), length.out = 25)
    z <- mpdWeightFit(y, x, x, xnew, xnew)
    expect_false(is.null(z))
    expect_true(all(diff(z) <= 1e-8))
    skip_if_not_installed("scam")
    fit <- scam::scam(y ~ s(x1, bs = "mpd") + s(x2, bs = "mpd"), data = data.frame(y = y, x1 = x, x2 = x))
    expected <- stats::predict(fit, data.frame(x1 = xnew, x2 = xnew))
    expect_equal(unname(z), as.vector(expected), tolerance = 1e-3)
})

test_that("topPairs() lists an image-method survival result, and chooses pairs", {
  data("diabetesData", package = "spicyR")
  d <- diabetesData
  set.seed(3)
  ev <- tapply(d$case, d$case, function(z) rexp(1))
  d$survival <- survival::Surv(ev[as.character(d$case)], rep(1, ncol(d)))
  res <- suppressMessages(spicy(d, condition = "survival", method = "image", from = "Tc", to = c("beta", "Th")))
  tp <- topPairs(res)
  expect_identical(names(tp), c("coefficient", "se", "p.value", "adj.pvalue", "from", "to"))
  expect_setequal(rownames(tp), c("Tc__beta", "Tc__Th"))
  expect_equal(tp$p.value, res$survivalResults$p.value[match(rownames(tp), res$survivalResults$test)])
  expect_identical(rownames(topPairs(res, to = "Th")), "Tc__Th")
})

test_that("a survival model that cannot be fitted for one pair gives NA for that pair, with a warning", {
  set.seed(1)
  pheno <- data.frame(subject = rep(seq_len(20), each = 2))
  pheno$surv <- survival::Surv(rep(rexp(20), each = 2), rep(rbinom(20, 1, 0.7), each = 2))
  M <- cbind(a__b = rnorm(40), a__c = NA_real_, b__c = rnorm(40))
  for (subject in list(NULL, "subject")) {
    if (!is.null(subject)) skip_if_not_installed("coxme")
    expect_warning(res <- spatialSurv(M, "surv", pheno, subject = subject),
                   "could not be fitted for 1 of 3 pairs; their results are NA: a__c")
    expect_setequal(res$test, colnames(M))
    expect_true(all(is.na(res[res$test == "a__c", c("coef", "se.coef", "p.value")])))
    expect_false(anyNA(res[res$test != "a__c", c("coef", "se.coef", "p.value")]))
  }
})

test_that("the image method stops when an image ID is reused across patients", {
  cells <- sim_cells(seed = 2, n_pat = 4, n_cells = 150)
  cells$roi <- sub(".*_", "ROI_", cells$imageID)   # ROI_1 and ROI_2 in every patient
  expect_error(suppressWarnings(suppressMessages(
    spicy(cells, condition = "condition", subject = "patient", imageID = "roi", from = "tumour", to = "T",
          method = "image"))), "'condition' must be constant within each image")
})

test_that("topPairs() of an image-method result with no pair tested says so", {
  cells <- sim_cells(seed = 2, n_pat = 4, n_cells = 150)
  cells$condition <- "A"
  res <- suppressWarnings(suppressMessages(spicy(cells, condition = "condition", from = "tumour", to = "T",
                                                 method = "image")))
  expect_error(topPairs(res), "no pair could be tested")
})

test_that("the image method gives NA, without stopping, for an image whose cells span no area", {
  cells <- sim_cells(seed = 2, n_pat = 4, n_cells = 150)
  one <- cells$imageID == "A01_1"
  cells <- cells[!one | !duplicated(one), ]                 # a single cell in A01_1
  two <- which(cells$imageID == "A02_1")
  cells <- cells[-two[-(1:2)], ]                           # two cells in A02_1
  line <- which(cells$imageID == "A03_1")
  cells <- cells[-line[-(1:3)], ]                          # three cells on a line in A03_1
  cells$y[cells$imageID == "A03_1"] <- cells$x[cells$imageID == "A03_1"]
  for (window in c("square", "convex", "concave")) for (sigma in list(NULL, 50)) {
    pa <- suppressWarnings(suppressMessages(getPairwise(cells, r = c(10, 30), window = window, sigma = sigma,
                                                        from = "tumour", to = c("T", "tumour"))))
    expect_true(all(is.na(pa["A01_1", ])))
    if (window == "convex") expect_true(all(is.na(pa[c("A02_1", "A03_1"), ])))
    expect_false(anyNA(pa["B01_1", ]))
  }
  res <- suppressWarnings(suppressMessages(spicy(cells, condition = "condition", subject = "patient", sigma = 50,
                                                 from = "tumour", to = "T", method = "image")))
  expect_false(is.na(topPairs(res)$p.value))
})

test_that("getPairwise() reads a data.frame's columns by the names given", {
  cells <- sim_cells(seed = 2, n_pat = 2, n_cells = 150)
  renamed <- cells
  names(renamed)[match(c("imageID", "cellType", "x", "y"), names(renamed))] <- c("sample", "type", "px", "py")
  for (window in c("convex", "concave")) {
    a <- suppressMessages(getPairwise(cells, r = c(10, 30), window = window))
    b <- suppressMessages(getPairwise(renamed, imageID = "sample", cellType = "type", spatialCoords = c("px", "py"),
                                      r = c(10, 30), window = window))
    expect_equal(b, a)
  }
})

test_that("getPairwise() matches cellWeights to the cells when the images are not in order", {
  skip_if_not_installed("SingleCellExperiment")
  cells <- sim_cells(seed = 2, n_pat = 2, n_cells = 150)
  cells <- cells[sample(nrow(cells)), ]
  sce <- SingleCellExperiment::SingleCellExperiment(colData = S4Vectors::DataFrame(cells))
  lev <- sort(unique(cells$cellType))
  onehot <- outer(cells$cellType, lev, "==") * 1
  colnames(onehot) <- lev
  a <- getPairwise(sce, r = c(10, 30), window = "square")
  b <- getPairwise(sce, r = c(10, 30), window = "square", cellWeights = onehot)
  expect_equal(b[rownames(a), colnames(a)], a)
})
