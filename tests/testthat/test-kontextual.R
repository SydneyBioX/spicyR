## kontextualEngine() (the engine of Statial::kontextualTest()): the C++ core against a plain-R reference written
## from Statial's Kontextual statistic, the
## exactness of its random-labelling moments, and the interface.

cells <- sim_cells(seed = 4)
parent <- c("T", "B", "macro")

## Statial::parentCombinations()'s format: every type as `from`, every type of the parent as `to`
make_parent_df <- function(types, parent, name) {
  d <- expand.grid(from = types, to = parent, stringsAsFactors = FALSE)
  d$parent <- rep(list(parent), nrow(d)); d$parent_name <- name; d
}
pdf <- make_parent_df(unique(cells$cellType), parent, "lymphoid")

## Plain R, one image, no edge correction: Statial's weights (parent counts include the cell itself) and the
## Kontextual excess's image row on Statial's scale.
ref_kontextual_row <- function(z, from, to, parent, r) {
  nb <- ref_neighbours(z$x, z$y, r); near <- nb$adj; diag(near) <- TRUE
  isP <- z$cellType %in% parent; isF <- z$cellType == from; isT <- z$cellType == to
  lam <- colSums(near[isP, , drop = FALSE]) / (pi * r^2)
  nT <- sum(isT); nP <- sum(isP)
  s <- colSums(nb$adj[isF, , drop = FALSE] * lam[isF]) / (lam * nT / nP)
  cand <- if (from %in% parent) isP & !isF else isP
  M <- sum(cand); p <- nT / M; sm <- s[cand]
  list(O = sum(s[isT]), E = p * sum(sm), n = sum(lam[isF]), v = p * (1 - p) * M / (M - 1) * sum((sm - mean(sm))^2),
       s = sm, nT = nT, L = sqrt(sum(s[isT]) / sum(lam[isF]) / pi) - r)
}

test_that("image rows match the plain-R reference (no edge correction, no label clustering)", {
  ctx <- spicyR:::.cell_context(cells, "condition", "patient", "imageID", "cellType", c("x", "y"))
  code <- function(z) match(z, ctx$type_labels) - 1L
  spicyR:::dataset_build_radius_index(ctx$data, 30)
  spicyR:::dataset_build_context(ctx$data, code(parent), "convex", FALSE)
  for (pr in list(c("tumour", "T"), c("macro", "B"))) {
    S <- spicyR:::dataset_kontextual_sums(ctx$data, code(pr[1]), code(pr[2]))
    rows <- spicyR:::stats_kontextual_image_rows(S, ctx$counts, code(pr[1]), code(pr[2]), matrix(numeric(0), 0, 0))
    for (k in c(1, 7)) {
      z <- ctx$df[ctx$image_codes == rows$img[k], ]
      ref <- ref_kontextual_row(z, pr[1], pr[2], parent, 30)
      expect_equal(c(rows$O[k], rows$E[k], rows$n[k], rows$v[k]), c(ref$O, ref$E, ref$n, ref$v), tolerance = 1e-10)
    }
  }
})

test_that("the Kontextual statistic's random-labelling moments are exact (enumeration)", {
  set.seed(7)
  z <- data.frame(x = runif(16, 0, 60), y = runif(16, 0, 60), cellType = sample(c("A", "B", "C", "D"), 16, TRUE))
  z$cellType[1:3] <- c("A", "B", "C")
  ref <- ref_kontextual_row(z, "A", "B", c("B", "C"), 25)
  en <- ref_rl_enumerate(ref$s, ref$nT)
  expect_equal(unname(en[["mean"]]), ref$E, tolerance = 1e-12)
  expect_equal(unname(en[["var"]]), ref$v, tolerance = 1e-12)
})

test_that("kontextualEngine (psi and abundance off) is the frailty GEE on the reference rows", {
  res <- kontextualEngine(cells, pdf, condition = "condition", subject = "patient",
                        r = 30, from = "tumour", to = "T", adjustAbundance = FALSE, labelClustering = FALSE,
                        edgeCorrect = FALSE)$cellResults
  ids <- sort(unique(cells$imageID), method = "radix")
  rows <- do.call(rbind, lapply(seq_along(ids), function(i) {
    r <- ref_kontextual_row(cells[cells$imageID == ids[i], ], "tumour", "T", parent, 30)
    data.frame(img = i - 1L, O = r$O, E = r$E, n = r$n, v = r$v) }))
  pat <- cells$patient[match(ids, cells$imageID)]; grp <- cells$condition[match(ids, cells$imageID)]
  rows$unit <- match(pat, unique(pat)) - 1L; rows$group <- as.integer(grp == "B")
  ref <- spicyR:::stats_excess_test(rows, rows$unit, rows$group, length(unique(pat)), TRUE, "cr2")
  expect_equal(res$p_value, ref$p, tolerance = 1e-10)
  expect_equal(res$excess_difference, ref$difference, tolerance = 1e-10)
  # T cells were moved next to tumour cells in group B: more tumour cells around T than around the other lymphoid cells
  expect_gt(res$excess_difference, 0); expect_lt(res$p_value, 0.05)
})

test_that("the interface: triples are labelled, plots work", {
  pl <- kontextualEngine(cells, pdf, condition = "condition", subject = "patient", r = 30)
  expect_true(all(pl$cellResults$parent == "lymphoid"))
  expect_true("tumour__T__lymphoid" %in% rownames(pl$cellResults))
  # no parent_name: the parent is named by its types
  un <- pdf; un$parent_name <- NULL
  expect_true(paste0("tumour__T__", paste(sort(parent), collapse = "+")) %in% rownames(kontextualEngine(cells, un, condition = "condition", subject = "patient",
                                                                    r = 30, from = "tumour")$cellResults))
  tp <- topPairs(pl, n = 3)
  expect_true("parent" %in% names(tp))
  expect_s3_class(spicyBoxPlot(pl, from = "tumour", to = "T", parent = "lymphoid"), "ggplot")
  expect_s3_class(signifPlot(pl), "ggplot")
  expect_error(kontextualEngine(cells, make_parent_df("B", c("B", "macro"), "p"), condition = "condition", r = 30, to = "T"),
               "no triple")
})

test_that("survival outcomes: the score test on the reference rows", {
  set.seed(9)
  pats <- unique(cells$patient); tm <- stats::setNames(rexp(length(pats), 0.1), pats); ev <- stats::setNames(rbinom(length(pats), 1, 0.7), pats)
  cs <- cells; cs$surv <- survival::Surv(tm[cs$patient], ev[cs$patient])
  res <- kontextualEngine(cs, pdf, condition = "surv", subject = "patient", r = 30,
                        from = "tumour", to = "T", adjustAbundance = FALSE, labelClustering = FALSE, edgeCorrect = FALSE)
  tab <- res$cellResults
  expect_equal(tab$parent, "lymphoid")
  ids <- sort(unique(cs$imageID), method = "radix")
  rows <- do.call(rbind, lapply(seq_along(ids), function(i) {
    r <- ref_kontextual_row(cs[cs$imageID == ids[i], ], "tumour", "T", parent, 30)
    data.frame(img = i - 1L, O = r$O, E = r$E, n = r$n, v = r$v) }))
  pat <- cs$patient[match(ids, cs$imageID)]; up <- unique(pat); rows$unit <- match(pat, up) - 1L
  null <- spicyR:::stats_cox_fit(unname(tm[up]), unname(ev[up]), matrix(0, length(up), 0))
  ref <- spicyR:::stats_survival_test(rows, rows$unit, length(up), null$martingale, unname(tm[up]), unname(ev[up]), numeric(0))
  expect_equal(tab$p_value, ref$score_p, tolerance = 1e-10)
})

test_that("signifPlot draws survival results of the cell method", {
  set.seed(9)
  pats <- unique(cells$patient); cs <- cells
  cs$surv <- survival::Surv(stats::setNames(rexp(length(pats), 0.1), pats)[cs$patient], rep(1, nrow(cs)))
  expect_s3_class(signifPlot(spicy(cs, condition = "surv", subject = "patient", r = 30, from = "tumour")), "ggplot")
  expect_s3_class(signifPlot(kontextualEngine(cs, pdf, condition = "surv",
                                            subject = "patient", r = 30)), "ggplot")
})

test_that("a triple's result does not depend on the other triples requested (psi over every from type)", {
  full <- kontextualEngine(cells, pdf[pdf$from != pdf$to, ], condition = "condition", subject = "patient", r = 30)$cellResults
  one <- kontextualEngine(cells, pdf, condition = "condition", subject = "patient", r = 30, from = "tumour", to = "T")$cellResults
  expect_equal(one["tumour__T__lymphoid", "p_value"], full["tumour__T__lymphoid", "p_value"], tolerance = 1e-12)
  expect_equal(one["tumour__T__lymphoid", "excess_difference"], full["tumour__T__lymphoid", "excess_difference"], tolerance = 1e-12)
})
