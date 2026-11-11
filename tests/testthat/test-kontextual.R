## kontextualEngine() (the engine of Statial::kontextualTest()): the C++ core against a plain-R reference written
## from Statial's Kontextual statistic, the exactness of its random-labelling moments (self triples included), the
## image random effect and abundance adjustment as spicy()'s, and the interface.

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
  ctx <- spicyR:::.cell_context(cells, "condition", "patient", "imageID", "cellType", c("x", "y"), dataset = TRUE)
  code <- function(z) match(z, ctx$type_labels) - 1L
  spicyR:::dataset_build_radius_index(ctx$data, 30)
  spicyR:::dataset_build_context(ctx$data, code(parent), "convex", FALSE)
  for (pr in list(c("tumour", "T"), c("macro", "B"))) {
    S <- spicyR:::dataset_kontextual_sums(ctx$data, code(pr[1]), code(pr[2]))
    rows <- spicyR:::stats_kontextual_image_rows(S, ctx$counts, code(pr[1]), code(pr[2]), matrix(numeric(0), 0, 0),
                                                 matrix(numeric(0), 0, 0))
    for (k in c(1, 7)) {
      z <- cells[as.character(cells$imageID) == ctx$image_labels[rows$img[k] + 1L], ]   # the image's cells, in input order
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

## The reference rows of every image (input order of the image IDs as the package sorts them), with unit and group.
ref_kontextual_rows <- function(cells, from, to, parent, r) {
  ids <- sort(unique(cells$imageID), method = "radix")
  rows <- do.call(rbind, lapply(seq_along(ids), function(i) {
    z <- cells[cells$imageID == ids[i], ]
    x <- ref_kontextual_row(z, from, to, parent, r)
    data.frame(img = i - 1L, O = x$O, E = x$E, n = x$n, v = x$v, n_to = sum(z$cellType == to),
               n_parent = sum(z$cellType %in% parent)) }))
  pat <- cells$patient[match(ids, cells$imageID)]; grp <- cells$condition[match(ids, cells$imageID)]
  rows$unit <- match(pat, unique(pat)) - 1L; rows$group <- as.integer(grp == "B")
  rows
}

test_that("kontextualEngine (psi and abundance off) is the frailty GEE on the reference rows, with the image random effect", {
  res <- kontextualEngine(cells, pdf, condition = "condition", subject = "patient",
                        r = 30, from = "tumour", to = "T", labelClustering = FALSE, edgeCorrect = FALSE)$cellResults
  rows <- ref_kontextual_rows(cells, "tumour", "T", parent, 30)
  # each image's effect (O - E) / D with variance v / D^2 + omega2 (Paule-Mandel within patients), as spicy()
  h <- ref_image_heterogeneity(rows, rows$unit)
  expect_gt(attr(h, "tau2_image"), 0)
  expect_equal(res$tau2_image, attr(h, "tau2_image"), tolerance = 1e-8)
  ref <- spicyR:::stats_excess_test(h, h$unit, h$group, length(unique(rows$unit)), TRUE, "cr2", FALSE)
  expect_equal(res$p_value, ref$p, tolerance = 1e-8)
  expect_equal(res$excess_difference, ref$difference, tolerance = 1e-8)
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
  rows <- ref_image_heterogeneity(rows, rows$unit)   # the image random effect, as spicy()'s survival test
  null <- spicyR:::stats_cox_fit(unname(tm[up]), unname(ev[up]), matrix(0, length(up), 0))
  ref <- spicyR:::stats_survival_test(rows, rows$unit, length(up), null$martingale, unname(tm[up]), unname(ev[up]),
                                      matrix(0, nrow(rows), 0), matrix(0, length(up), 0))
  expect_equal(tab$p_value, ref$score_p, tolerance = 1e-8)
})

test_that("signifPlot draws survival results of the cell method", {
  set.seed(9)
  pats <- unique(cells$patient); cs <- cells
  cs$surv <- survival::Surv(stats::setNames(rexp(length(pats), 0.1), pats)[cs$patient], rep(1, nrow(cs)))
  expect_s3_class(signifPlot(spicy(cs, condition = "surv", subject = "patient", r = 30, from = "tumour")), "ggplot")
  expect_s3_class(signifPlot(kontextualEngine(cs, pdf, condition = "surv",
                                            subject = "patient", r = 30)), "ggplot")
})

test_that("with labelClustering, a triple's result does not depend on the other triples requested (psi over every from type)", {
  full <- kontextualEngine(cells, pdf[pdf$from != pdf$to, ], condition = "condition", subject = "patient", r = 30,
                           labelClustering = TRUE)$cellResults
  one <- kontextualEngine(cells, pdf, condition = "condition", subject = "patient", r = 30, from = "tumour", to = "T",
                          labelClustering = TRUE)$cellResults
  expect_equal(one["tumour__T__lymphoid", "p_value"], full["tumour__T__lymphoid", "p_value"], tolerance = 1e-12)
  expect_equal(one["tumour__T__lymphoid", "excess_difference"], full["tumour__T__lymphoid", "excess_difference"], tolerance = 1e-12)
})

test_that("abundance adjustment: off by default; when on, the log share of the `to` type within the parent", {
  expect_false(formals(kontextualEngine)$adjustAbundance)
  res <- kontextualEngine(cells, pdf, condition = "condition", subject = "patient", r = 30, from = "tumour", to = "T",
                          edgeCorrect = FALSE, adjustAbundance = TRUE)$cellResults
  expect_equal(res$adjusted_for, "abundance")
  rows <- ref_kontextual_rows(cells, "tumour", "T", parent, 30)
  h <- ref_image_heterogeneity(rows, rows$unit)
  share <- log(pmax(rows$n_to, 0.5) / rows$n_parent)
  # the adjusted design of spicy() (ext-AB): one indicator per condition and the centred share, which varies within
  # patients, so it gets its patient mean, omega2 by B2 and tau2 on the patient summaries
  Z <- cbind(rows$group == 0, rows$group == 1, share - mean(share)) * 1
  a <- ref_adjusted_design(rows, rows$unit, Z, 2)
  cv <- c(-1, 1, rep(0, ncol(a$Z) - 2))
  ref <- ref_design(a$rows, rows$unit, a$Z, cv, a$tau2)
  w <- ref_design(a$rows, rows$unit, a$Z, replace(0 * cv, 3, 1), a$tau2)
  expect_equal(res$p_value, unname(ref[4]), tolerance = 1e-8)
  expect_equal(res$abundance_effect, unname(w[1]), tolerance = 1e-8)
  expect_equal(c(res$tau2, res$tau2_image), c(a$tau2, a$omega2), tolerance = 1e-8)
})

test_that("one radius: several radii are not supported", {
  expect_error(kontextualEngine(cells, pdf, condition = "condition", subject = "patient", r = c(20, 30)),
               "several radii are not supported")
})

## Every labelling of n `A` cells among the parent's cells of one small image, each labelling as its own image (same
## positions, so the same window, weights and exact variance): the package's v of a self triple against the variance of
## f (O - p L) over the labellings.
self_enumeration <- function(x, y, par, r, n, window, edge) {
  M <- sum(par); S <- utils::combn(M, n); if (ncol(S) == 1L) S <- cbind(S, S)
  z <- do.call(rbind, lapply(seq_len(ncol(S)), function(k) {
    ty <- ifelse(par, "P", "Z"); ty[which(par)[S[, k]]] <- "A"
    data.frame(imageID = sprintf("i%04d", k), cellType = ty, x = x, y = y, patient = sprintf("i%04d", k),
               condition = c("a", "b")[k %% 2 + 1]) }))
  ctx <- spicyR:::.cell_context(z, "condition", "patient", "imageID", "cellType", c("x", "y"), dataset = TRUE)
  code <- function(t) match(t, ctx$type_labels) - 1L
  spicyR:::dataset_build_radius_index(ctx$data, r)
  spicyR:::dataset_build_context(ctx$data, code(intersect(c("A", "P"), ctx$type_labels)), window, edge)
  s <- spicyR:::dataset_kontextual_sums(ctx$data, code("A"), code("A"))
  rows <- spicyR:::stats_kontextual_image_rows(s, ctx$counts, code("A"), code("A"), matrix(numeric(0), 0, 0),
                                               spicyR:::dataset_kontextual_self_sums(ctx$data))
  R <- M / n * (s[1, ] - (n - 1) / (M - 1) * s[2, ])
  list(rows = rows, var = mean((R - mean(R))^2), mean = mean(R), images = ncol(S))
}

test_that("self triples: the exact random-labelling variance of O - p L (enumeration), and the zero-variance rule", {
  set.seed(3)
  x <- c(runif(5, 0, 30), runif(5, 25, 80), runif(4, 0, 80)); y <- c(runif(5, 0, 30), runif(5, 0, 60), runif(4, 0, 60))
  par <- rep(c(TRUE, FALSE), c(10, 4))
  for (cfg in list(list("convex", TRUE, 2), list("convex", TRUE, 5), list("rectangle", TRUE, 8), list("convex", FALSE, 4))) {
    en <- self_enumeration(x, y, par, 30, cfg[[3]], cfg[[1]], cfg[[2]])
    expect_equal(nrow(en$rows), en$images)
    expect_equal(en$rows$v, rep(en$var, en$images), tolerance = 1e-12)
    expect_lt(abs(en$mean), 1e-10 * sqrt(en$var))   # O - p L has mean 0: E is exact
  }
  # every parent cell an A cell (n = M): every labelling gives the same O - p L; exact variance 0, so no row
  expect_equal(nrow(self_enumeration(x, y, par, 30, 10, "convex", TRUE)$rows), 0L)
  # n = 1: O = 0 for every labelling; no row
  expect_equal(nrow(self_enumeration(x, y, par, 30, 1, "convex", TRUE)$rows), 0L)
})

test_that("self triples through the engine: the rows with the exact variance, the test as for other triples", {
  sp <- make_parent_df("T", parent, "lymphoid")
  res <- kontextualEngine(cells, sp[sp$to == "T", ], condition = "condition", subject = "patient", r = 30,
                          edgeCorrect = FALSE)$cellResults
  expect_true("T__T__lymphoid" %in% rownames(res))
  ctx <- spicyR:::.cell_context(cells, "condition", "patient", "imageID", "cellType", c("x", "y"), dataset = TRUE)
  code <- function(z) match(z, ctx$type_labels) - 1L
  spicyR:::dataset_build_radius_index(ctx$data, 30)
  spicyR:::dataset_build_context(ctx$data, code(parent), "convex", FALSE)
  S <- spicyR:::dataset_kontextual_sums(ctx$data, code("T"), code("T"))
  expect_error(spicyR:::stats_kontextual_image_rows(S, ctx$counts, code("T"), code("T"), matrix(numeric(0), 0, 0),
                                                    matrix(numeric(0), 0, 0)), "self triple")
  rows <- spicyR:::stats_kontextual_image_rows(S, ctx$counts, code("T"), code("T"), matrix(numeric(0), 0, 0),
                                               spicyR:::dataset_kontextual_self_sums(ctx$data))
  rows$unit <- ctx$image_unit[rows$img + 1L]; rows$group <- ctx$image_group[rows$img + 1L]
  h <- ref_image_heterogeneity(rows, rows$unit)
  ref <- spicyR:::stats_excess_test(h, h$unit, h$group, length(ctx$unit_labels), TRUE, "cr2", FALSE)
  expect_equal(res["T__T__lymphoid", "p_value"], ref$p, tolerance = 1e-8)
})

test_that("more than two conditions: a follow-up is kontextualEngine() on the patients of its two levels", {
  # Kontextual rows go through spicy()'s image random effect (.cell_image_effect), so each follow-up re-estimates
  # tau2_image on its two levels' patients, and the abundance share stays Kontextual's
  k3 <- sim_cells(n_pat = 16, seed = 4); k3$condition <- as.character(k3$condition)   # (enough power to call a pair)
  k3$condition[k3$patient %in% sprintf("A%02d", 1:8)] <- "C"
  p3 <- make_parent_df(unique(k3$cellType), parent, "lymphoid")
  for (adjust in c(FALSE, TRUE)) {
    res <- kontextualEngine(k3, p3, condition = "condition", subject = "patient", r = 30, adjustAbundance = adjust)
    fu <- res$followUpResults
    expect_gt(nrow(fu), 0)
    for (lv in list(c("A", "B"), c("A", "C"), c("B", "C"))) {
      sub <- k3[k3$condition %in% lv, ]; sub$condition <- factor(sub$condition, levels = lv)
      two <- kontextualEngine(sub, p3, condition = "condition", subject = "patient", r = 30,
                              adjustAbundance = adjust)$cellResults
      f <- fu[fu$level_ref == lv[1] & fu$level_comp == lv[2], ]
      if (!nrow(f)) next
      t2 <- two[match(paste(f$from, f$to, f$parent), paste(two$from, two$to, two$parent)), ]
      for (col in c("excess_difference", "se", "df", "p_value", "tau2", "tau2_image"))
        expect_equal(f[[col]], t2[[col]], tolerance = 1e-8, label = paste(adjust, paste(lv, collapse = ""), col))
    }
  }
})
