## Results saved by an earlier spicyR (reference_*.rds) predate the skewness correction of the two-condition CR2 t,
## which changes only p-values. Every other column must be unchanged; each p_value must equal the reference either
## as it is (tests without the correction: adjusted designs, more conditions, continuous) or as the uncorrected
## p = 2 pt(-|difference / se|, df) of the same row; the other p columns (p_adj, combined and unadjusted p) keep their
## missing values. combined = TRUE (p_value combined over several radii, from corrected t's): p_value only to within
## a relative 0.1.
expect_unchanged_but_skewness <- function(now, ref, label = "", tolerance = 1e-12, combined = FALSE) {
  pcols <- grep("p_value|p_adj", names(ref), value = TRUE)
  other <- setdiff(names(ref), pcols)
  expect_identical(intersect(names(ref), names(now)), names(ref), label = label)
  expect_equal(now[other], ref[other], tolerance = tolerance, label = label)
  for (cl in pcols) expect_identical(is.na(now[[cl]]), is.na(ref[[cl]]), label = paste(label, cl))
  if (combined) {
    expect_equal(now$p_value, ref$p_value, tolerance = 0.1, label = paste(label, "p_value"))
  } else if ("p_value" %in% pcols && all(c("excess_difference", "se", "df") %in% names(now))) {
    p_unc <- 2 * stats::pt(-abs(now$excess_difference / now$se), now$df)
    same <- abs(now$p_value - ref$p_value) <= 1e-10 | abs(p_unc - ref$p_value) <= 1e-10
    expect_true(all(same | is.na(ref$p_value)), label = paste(label, "p_value"))
  }
}

## The skewness of the two-condition test, shrunk by its own uncertainty, from the units' standardised residuals z and
## leverages h (plain R, from the definitions): g = (sum z^3 / sum (1 - h)(1 - 2h)) / (sum z^2 / sum (1 - h))^1.5; its
## delta-method variance V = mean(IF^2) / n, IF = t^3 - 3t - g - 1.5 g (t^2 - 1), t = z / sigma, n = sum (1 - h) + 1;
## shrunk g^3 / (g^2 + V); 0 when a condition (g, if given) has fewer than 3 units.
ref_shrunk_skewness <- function(z, h, g = NULL) {
  if (!is.null(g) && (length(unique(g)) < 2 || any(table(g) < 3))) return(0)
  d3 <- sum((1 - h) * (1 - 2 * h)); if (!(d3 > 0)) return(0)
  s2 <- sum(z^2) / sum(1 - h); g <- (sum(z^3) / d3) / s2^1.5
  t <- z / sqrt(s2); IF <- t^3 - 3 * t - g - 1.5 * g * (t^2 - 1)
  V <- mean(IF^2) / (sum(1 - h) + 1)
  g^3 / (g^2 + V)
}
