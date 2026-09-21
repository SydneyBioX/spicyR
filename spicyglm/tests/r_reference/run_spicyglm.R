# Shared R runner for the reference fixtures and the benchmarks.
#
# Fits every cell-type pair with the R implementation (spicyR, gee branch)
# through the same modelDataGen() / buildGLM() path spicyGLM() uses, and writes
# results_<name>.csv, skipped_<name>.csv and, with diagnostics,
# diag_<table>_<name>.csv. Returns the elapsed fitting time in seconds
# (writing CSVs excluded).

load_spicyR <- function(spicyR) {
  suppressPackageStartupMessages({
    library(methods)
    library(spatstat.geom)
    library(dplyr)
  })
  for (f in c("spicyGLM.R", "spicyGLM_meanModel.R", "cr2_fast_multi.R", "dpr1eig.R")) {
    # the final setMethod() in spicyGLM.R needs a class this runner never uses
    tryCatch(source(file.path(spicyR, "R", f)), error = function(e) NULL)
  }
}

# For validation only: refit iteratively fitted pairs (binomial, and Poisson
# through glm / brglm2, which cr2Method = "naive" uses) to a tight convergence
# tolerance so the comparison is not limited by glm's / brglmFit's default
# stopping rule. Same models and estimators; timings are always taken without this.
tighten_iterative_fits <- function() {
  assign("fit_mle_glm", function(dfPair) {
    glm(n ~ 0 + condition, offset = log(density), family = poisson(), data = dfPair,
        epsilon = 1e-14, maxit = 1000)
  }, envir = globalenv())
  assign("fit_firth_brglm2", function(dfPair) {
    glm(n ~ 0 + condition, offset = log(density), family = poisson(), data = dfPair,
        method = brglm2::brglmFit, type = "AS_mean", epsilon = 1e-14, maxit = 1000)
  }, envir = globalenv())
  assign("fit_mle_binom_glm", function(dfPair) {
    glm(cbind(n, k - n) ~ 0 + condition, offset = qlogis(p0), family = binomial(), data = dfPair,
        epsilon = 1e-14, maxit = 1000)
  }, envir = globalenv())
  assign("fit_firth_binom_brglm2", function(dfPair) {
    glm(cbind(n, k - n) ~ 0 + condition, offset = qlogis(p0), family = binomial(), data = dfPair,
        method = brglm2::brglmFit, type = "AS_mean", epsilon = 1e-14, maxit = 1000)
  }, envir = globalenv())
}

run_spicyglm <- function(cells, out_dir, name, subject, window = "convex", estimator = "firth", r = NULL,
                         diagnostics = FALSE, family = "poisson", k = NULL, cr2Method = "fast") {
  start <- Sys.time()
  one_to_one <- is.null(subject) || length(unique(cells[[subject]])) == length(unique(cells$imageID))
  presence <- computeCellTypePresence(cells, condition = "condition", imageID = "imageID", cellType = "cellType")
  types <- unique(cells$cellType)
  pairs <- rbind(do.call(rbind, lapply(combn(types, 2, simplify = FALSE),
                                       function(p) data.frame(from = p[1], to = p[2]))),
                 data.frame(from = types, to = types))
  fits <- list(); skips <- list(); diags <- list()
  for (i in seq_len(nrow(pairs))) {
    from <- pairs$from[i]; to <- pairs$to[i]
    df <- suppressWarnings(modelDataGen(cells, condition = "condition", subject = subject, from = from, to = to,
                                        r = r, imageID = "imageID", cellType = "cellType",
                                        spatialCoords = c("x", "y"), window = window, oneToOne = one_to_one,
                                        cellTypePresence = presence, family = family, k = k))
    if (nrow(df) == 0) {
      skips[[length(skips) + 1]] <- data.frame(from = from, to = to, reason = attr(df, "skipReason"))
      next
    }
    df$from <- from; df$to <- to
    res <- suppressWarnings(buildGLM(df, oneToOne = one_to_one, subject = subject, estimator = estimator,
                                     computeDiagnostics = diagnostics, cellTypePresence = presence,
                                     family = family, k = k, cr2Method = cr2Method))
    if ("reason" %in% names(res)) {
      skips[[length(skips) + 1]] <- data.frame(from = from, to = to, reason = res$reason)
    } else {
      res$from <- from; res$to <- to
      fits[[length(fits) + 1]] <- res
      if (diagnostics) diags[[length(diags) + 1]] <- attr(res, "diagnostics")
    }
  }
  fit <- bind_rows(fits)
  fit$p.adj <- p.adjust(fit$p.value, method = "fdr")
  if (diagnostics) {
    all <- list(pair    = bind_rows(lapply(diags, `[[`, "pair")),
                patient = bind_rows(lapply(diags, `[[`, "patient")),
                image   = bind_rows(lapply(diags, `[[`, "image")))
    patient <- computeRelativeDiagnostics(all, level = "patient")
    image <- computeRelativeDiagnostics(all, level = "image")
    tables <- list(pair = all$pair, patient = patient, image = image,
                   cross_patient = crossPairDiagnostics(patient), cross_image = crossPairDiagnostics(image))
  }
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  write.csv(fit, file.path(out_dir, paste0("results_", name, ".csv")), row.names = FALSE)
  write.csv(bind_rows(skips), file.path(out_dir, paste0("skipped_", name, ".csv")), row.names = FALSE)
  if (diagnostics) {
    for (tbl in names(tables)) {
      write.csv(tables[[tbl]], file.path(out_dir, paste0("diag_", tbl, "_", name, ".csv")), row.names = FALSE)
    }
  }
  elapsed
}
