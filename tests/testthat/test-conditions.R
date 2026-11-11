## More than two conditions (DECISIONS 4 Nov 2026, problem C): the overall test of each pair (approximate Hotelling's
## T^2, CR2), BH over the pairs, and the gated follow-ups (every two levels by the two-condition test on their
## patients, Bonferroni within the pair at the Benjamini-Bogomolov level alpha R / m).

## Three conditions: A (4 patients, no effect), C (4, no effect), B (8, T cells moved next to tumour cells).
c3 <- sim_cells(n_pat = 8, seed = 4); c3$condition <- as.character(c3$condition)
c3$condition[c3$patient %in% sprintf("A%02d", 1:4)] <- "C"
run3 <- function(...) suppressMessages(spicy(c3, "condition", subject = "patient", r = 30, ...))

## The image rows of a pair with the image random effect, as spicy() takes them
rows3 <- function(cells, f, t, effect = "allocation", ref = NULL) {
  ctx <- spicyR:::.cell_context(cells, "condition", "patient", "imageID", "cellType", c("x", "y"), ref = ref)
  g <- spicyR:::.cell_graph(ctx, list(c(f, t)), r = 30, label_clustering = FALSE, effect = effect)
  g$image_heterogeneity <- TRUE
  list(ctx = ctx, rows = spicyR:::.cell_rows(ctx, g, f, t))
}

test_that("the overall test is the approximate Hotelling's T^2 of equal means, whatever the reference", {
  res <- run3()
  tab <- res$cellResults; first <- tab[tab$level == tab$level[1], ]
  tested <- first[!is.na(first$overall_p_value), ]
  expect_gt(nrow(tested), 10)
  # the closed form of the prototype (cell-means design, per-level CR2, AHT), pair by pair
  for (i in seq_len(nrow(tested))) {
    x <- rows3(c3, tested$from[i], tested$to[i])$rows
    a <- ref_aht_cellmeans(x, x$group, 3L)
    expect_equal(unname(c(tested$overall_F[i], tested$overall_df2[i], tested$overall_p_value[i])),
                 unname(a[c("F", "df2", "p")]), tolerance = 1e-7)
  }
  # reference-free: another reference level gives the same overall test (and the same pairs tested)
  other <- run3(ref = "C")$cellResults
  o <- other[other$level == other$level[1], ]
  o <- o[match(paste(first$from, first$to), paste(o$from, o$to)), ]
  expect_equal(o$overall_p_value, first$overall_p_value, tolerance = 1e-10)
  expect_equal(o$overall_F, first$overall_F, tolerance = 1e-10)
  # two levels: the closed form is the two-condition t test (F = t^2, eta its Satterthwaite df)
  two <- suppressMessages(spicy(sim_cells(), "condition", subject = "patient", r = 30, from = "tumour", to = "T"))$cellResults
  x <- rows3(sim_cells(), "tumour", "T")$rows
  a <- ref_aht_cellmeans(x, x$group, 2L)
  expect_equal(unname(a["p"]), two$p_value, tolerance = 1e-8)
  expect_equal(unname(a["F"]), (two$excess_difference / two$se)^2, tolerance = 1e-8)
  expect_equal(unname(a["df2"]), two$df, tolerance = 1e-8)
})

test_that("follow-ups: the pairs called by the overall test, every two levels, at the Benjamini-Bogomolov level", {
  for (alpha in c(0.05, 0.5)) {
    res <- run3(alpha = alpha)
    tab <- res$cellResults; first <- tab[tab$level == tab$level[1], ]
    called <- which(first$overall_p_adj <= alpha)
    R <- length(called); m <- sum(!is.na(first$overall_p_value))
    expect_gt(R, 0)
    fu <- res$followUpResults
    expect_identical(res$alpha, alpha)
    expect_setequal(unique(paste(fu$from, fu$to)), paste(first$from, first$to)[called])
    expect_equal(nrow(fu), 3 * R)
    expect_setequal(paste(fu$level_ref, fu$level_comp), c("A B", "A C", "B C"))
    # the gate: alpha R / (m K), K the comparisons tested in the pair
    K <- tapply(!is.na(fu$p_value), paste(fu$from, fu$to), sum)[paste(fu$from, fu$to)]
    expect_equal(fu$n_comparisons, as.vector(K))
    expect_equal(fu$p_threshold, alpha * R / (m * fu$n_comparisons))
    expect_identical(fu$significant, fu$p_value <= fu$p_threshold)
    expect_equal(fu$overall_p_adj, first$overall_p_adj[match(paste(fu$from, fu$to), paste(first$from, first$to))])
  }
  # tumour -> T: B differs from A and C, A and C do not differ
  fu <- run3()$followUpResults
  tt <- fu[fu$from == "tumour" & fu$to == "T", ]
  expect_true(all(tt$significant[tt$level_comp == "B" | tt$level_ref == "B"]))
  expect_false(tt$significant[tt$level_ref == "A" & tt$level_comp == "C"])
  # nothing called: the columns, no rows
  none <- suppressMessages(spicy(c3, "condition", subject = "patient", r = 30, from = "B", to = "macro"))
  expect_equal(nrow(none$followUpResults), 0L)
  expect_true(all(c("p_value", "p_threshold", "significant") %in% names(none$followUpResults)))
})

test_that("a follow-up is the two-condition test of spicy() on the patients of its two levels", {
  for (adjust in c(FALSE, TRUE)) for (effect in c("count", "allocation")) {
    res <- run3(alpha = 0.5, effect = effect, adjustAbundance = adjust)
    fu <- res$followUpResults
    expect_gt(nrow(fu), 0)
    for (lv in list(c("A", "B"), c("A", "C"), c("B", "C"))) {
      sub <- c3[c3$condition %in% lv, ]; sub$condition <- factor(sub$condition, levels = lv)
      two <- suppressMessages(spicy(sub, "condition", subject = "patient", r = 30, effect = effect,
                                    adjustAbundance = adjust))$cellResults
      f <- fu[fu$level_ref == lv[1] & fu$level_comp == lv[2], ]
      t2 <- two[match(paste(f$from, f$to), paste(two$from, two$to)), ]
      # allocation: the side of the two levels alone may differ from the side over every image (which the
      # follow-up keeps, as the overall test); compare where they agree
      same <- if (effect == "allocation") f$side == t2$side else rep(TRUE, nrow(f))
      expect_gt(sum(same), 0)
      for (col in c("excess_ref", "excess_comp", "excess_difference", "se", "df", "p_value", "tau2"))
        expect_equal(f[[col]][same], t2[[col]][same], tolerance = 1e-10, label = paste(effect, adjust, col))
      if (adjust) expect_equal(f$unadjusted_p_value[same], t2$unadjusted_p_value[same], tolerance = 1e-10)
    }
  }
})

test_that("follow-ups match the prototype's two-group test (cg_lib.R twogroup_pair)", {
  res <- run3(alpha = 0.5)
  fu <- res$followUpResults
  levels3 <- c("A", "B", "C")
  for (i in seq_len(nrow(fu))) {
    x <- rows3(c3, fu$from[i], fu$to[i])
    raw <- attr(x$rows, "raw"); ctx <- x$ctx
    a <- match(fu$level_ref[i], levels3) - 1L; b <- match(fu$level_comp[i], levels3) - 1L
    # the prototype: the patients of the two levels, recoded in order; omega2 again on them; the unadjusted test
    lev_unit <- ctx$image_group[match(seq_along(ctx$unit_labels) - 1L, ctx$image_unit)]
    sel <- which(lev_unit %in% c(a, b)) - 1L
    s <- raw[raw$group %in% c(a, b), , drop = FALSE]
    unit <- match(s$unit, sel) - 1L
    h <- spicyR:::stats_image_heterogeneity_rows(s, unit); h$unit <- unit
    grp <- as.integer(lev_unit[sel + 1L][unit + 1L] == b); h$group <- grp
    r <- spicyR:::stats_excess_test(h, h$unit, grp, length(sel), TRUE, "cr2", FALSE)
    expect_equal(fu$p_value[i], r$p, tolerance = 1e-12)
    expect_equal(fu$excess_difference[i], r$difference, tolerance = 1e-12)
  }
})

test_that("topFollowUps() lists the follow-ups; print() counts them; alpha is checked", {
  res <- run3(alpha = 0.5)
  top <- topFollowUps(res)
  expect_equal(nrow(top), nrow(res$followUpResults))
  expect_false(is.unsorted(top$overall_p_value))
  expect_true(all(topFollowUps(res, significant = TRUE)$significant))
  expect_true(all(topFollowUps(res, from = "tumour")$from == "tumour"))
  expect_equal(nrow(topFollowUps(res, n = 2)), 2L)
  expect_output(print(res), "Follow-up comparisons")
  expect_error(run3(alpha = 2), "alpha")
  two <- suppressMessages(spicy(sim_cells(), "condition", subject = "patient", r = 30, from = "tumour", to = "T"))
  expect_null(two$followUpResults)
  expect_error(topFollowUps(two), "more than two conditions")
  # several radii with more than two conditions stay refused
  expect_error(suppressMessages(spicy(c3, "condition", subject = "patient", r = c(20, 30))), "several radii")
})

test_that("two and three conditions: the results of spicyR 1.99.56 are unchanged", {
  ref <- readRDS(test_path("reference_1.99.56_conditions.rds"))
  cells <- sim_cells()
  run <- function(...) { r <- suppressMessages(spicy(...)); list(cellResults = r$cellResults, radiusResults = r$radiusResults) }
  expect_identical(run(cells, "condition", subject = "patient", r = 30), ref$two_allocation)
  expect_identical(run(cells, "condition", subject = "patient", r = 30, effect = "count"), ref$two_count)
  # with the abundance share (it varies between a patient's images, so ext-AB's adjusted design changes the adjusted
  # test): the unadjusted test and every other column are unchanged; the adjusted design is tested in
  # test-adjustedDesign.R
  same_except <- function(a, b, changed) {
    # (several radii with max-T: the new column maxT_converged)
    expect_identical(setdiff(names(a$cellResults), "maxT_converged"), names(b$cellResults))
    keep <- setdiff(names(b$cellResults), changed)
    expect_identical(a$cellResults[keep], b$cellResults[keep])
    expect_identical(a$radiusResults, b$radiusResults)
  }
  ab <- c("excess_ref", "excess_comp", "excess_difference", "se", "df", "p_value", "p_adj", "tau2", "tau2_image",
          "abundance_effect", "abundance_p_value", "overall_F", "overall_df2", "overall_p_value", "overall_p_adj")
  same_except(run(cells, "condition", subject = "patient", r = 30, adjustAbundance = TRUE), ref$two_adjusted, ab)
  # several radii: ext-D's max-T p-value (relative error about 1e-3; 1.99.56's had an absolute error of about 1e-4 and
  # was up to a few percent off at p near 0.003); every radius's test is unchanged
  rr <- run(cells, "condition", subject = "patient", r = c(20, 30, 40))
  same_except(rr, ref$two_radii, c("p_value", "p_adj"))
  expect_identical(rr$cellResults$maxT_converged, rep(TRUE, nrow(rr$cellResults)))
  expect_equal(rr$cellResults$p_value, ref$two_radii$cellResults$p_value, tolerance = 0.05)
  expect_identical(run(c3, "condition", subject = "patient", r = 30), ref$three_allocation)
  same_except(run(c3, "condition", subject = "patient", r = 30, effect = "count", adjustAbundance = TRUE),
              ref$three_count_adjusted, ab)
})
