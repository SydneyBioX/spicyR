# ---- Libraries ------------------------------------------------------------

# Core tidyverse & plotting
library(tidyverse)
library(patchwork)
library(cowplot)
library(gridExtra)
library(plotROC)
library(fields)

# Parallel & future
library(parallel)
library(BiocParallel)
library(future)
library(future.apply)
library(progressr)

# Spatial statistics
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)

# Spatial omics frameworks
library(SpatialExperiment)
library(spicyR)

# Modeling
library(glmmTMB)
library(lme4)
library(geepack)
library(Matrix)
library(emmeans)
library(clubSandwich)
library(sandwich)
library(fixest)

# Your custom functions (must include model_data_gen etc.)
source("/Users/shreyarao/Desktop/functions.R")

# ---- Constants ------------------------------------------------------------

set.seed(51773)
window <- owin(c(0, 1000), c(0, 1000))

nPatients <- 40      # half control, half case
nIm       <- 3       # images per patient
nSim      <- 500     # replicates per lambda

counts   <- seq(20, 400, by = 10)   # per-image expected counts
Rs_fixed <- 50                      # fixed-r for all methods
lambdas  <- 40                      # true scale control knob

REF    <- "B"  # anchor/reference throughout
TARGET <- "A"  # counted / target

# ---- Small-sample test choice ---------------------------------------------

CS_TEST <- "EDT"  # lock test choice (EDT vs HTZ etc.)

wald_p <- function(fit, L, V, test = CS_TEST) {
  as.numeric(
    clubSandwich::Wald_test(
      fit, constraints = L, vcov = V,
      test = test,
      tidy = TRUE
    )$p_val[1]
  )
}

# ---- Diagnostics helper ---------------------------------------------------

diagnose_dat_r <- function(dat, r, lambda = NA, rep_id = NA, msg = NULL, debug_dir = NULL) {
  message(sprintf(
    "[λ=%s rep=%s r=%s] %s | rows=%d",
    as.character(lambda),
    as.character(rep_id),
    as.character(r),
    ifelse(is.null(msg), "No message", msg),
    if (is.null(dat)) 0L else nrow(dat)
  ))
  
  if (!is.null(debug_dir)) {
    dir.create(debug_dir, showWarnings = FALSE, recursive = TRUE)
    fn <- sprintf("%s/model_diag_true_lambda%03s_rep%04s_r%03s.rds",
                  debug_dir, as.character(lambda), as.character(rep_id), as.character(r))
    saveRDS(list(
      data = dat,
      summary = list(lambda = lambda, replicate = rep_id, r = r, msg = msg)
    ), fn)
  }
  invisible(NULL)
}

# ---- Simulation helpers ---------------------------------------------------

simulate_once <- function(i, lambda, scenario = c("signal", "null_equal_coloc", "null_indep")) {
  scenario <- match.arg(scenario)
  set.seed(i)
  
  # bandwidth/scale per subject used for constructing B from A
  if (scenario == "signal") {
    g1 <- rpois(nPatients/2, lambda)
    g2 <- rpois(nPatients/2, lambda + lambda/3)
  } else if (scenario == "null_equal_coloc") {
    g1 <- rpois(nPatients/2, lambda)
    g2 <- rpois(nPatients/2, lambda)         # equal across groups
  } else { # null_indep
    g1 <- g2 <- integer(nPatients/2)         # unused, but define for completeness
  }
  adjustSigma <- if (scenario == "null_indep") rep(NA_real_, nPatients) else c(g1, g2) + 1
  
  x <- y <- cellType <- imageID <- condition <- subject <- NULL
  for (p in seq_len(nPatients)) {
    for (j in seq_len(nIm)) {
      repeat { sA <- rpois(1, sample(counts, 1)); if (sA > 0) break }
      repeat { sB <- rpois(1, sample(counts, 1)); if (sB > 0) break }
      
      A <- rpoispp(sA / area.owin(window), win = window)
      
      if (scenario %in% c("signal","null_equal_coloc")) {
        aDens <- spatstat.explore::density.ppp(A, sigma = adjustSigma[p], kernel = "disc")
        T <- spatstat.geom::integral.im(aDens)
        aDens$v <- pmax(aDens$v, 0) * (sB / T)
        B <- rpoispp(aDens)                   # B attracted to A with same/different strength
      } else { # null_indep
        B <- rpoispp(sB / area.owin(window), win = window)  # CSR independent of A
      }
      
      lbl <- paste(p, j, sep = "_")
      x   <- c(x, A$x, B$x)
      y   <- c(y, A$y, B$y)
      cellType  <- c(cellType, rep("A", A$n), rep("B", B$n))
      imageID   <- c(imageID, rep(lbl, A$n + B$n))
      group     <- if (p <= nPatients/2) "Group1" else "Group2"
      condition <- c(condition, rep(group, A$n + B$n))
      subject   <- c(subject, rep(p, A$n + B$n))
    }
  }
  
  cellExp <- data.frame(
    x, y,
    cellType  = factor(cellType),
    imageID   = factor(imageID),
    condition = factor(condition),
    subject   = factor(subject)
  )
  
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = matrix(0, nrow = 1, ncol = nrow(cellExp))),
    spatialCoords = cbind(x = cellExp$x, y = cellExp$y),
    colData = cellExp[, c("cellType", "imageID", "subject", "condition")]
  )
  
  list(cellExp = cellExp, spe = spe)
}

simulate_pair <- function(i, lambda, null_type = c("null_equal_coloc", "null_indep")) {
  null_type <- match.arg(null_type)
  list(
    sig  = simulate_once(i, lambda, scenario = "signal"),
    null = simulate_once(i, lambda, scenario = null_type)
  )
}

# ---- VCOV helpers ---------------------------------------------------------

vcov_geeglm_type <- function(fit,
                             type = c("naive", "CR0", "CR2", "CR3"),
                             cluster) {
  type <- match.arg(type)
  
  if (type == "naive") {
    return(fit$geese$vbeta.naiv)
  }
  
  if (missing(cluster)) stop("Provide 'cluster' for CR0/CR2/CR3.")
  
  clubSandwich::vcovCR(fit, cluster = cluster, type = type)
}

vcov_lmglm_cluster <- function(fit,
                               type = c("naive", "CR0", "CR2", "CR3"),
                               cluster) {
  type <- match.arg(type)
  
  if (type == "naive") return(stats::vcov(fit))
  if (missing(cluster)) stop("Provide 'cluster' for clustered vcov.")
  
  clubSandwich::vcovCR(fit, cluster = cluster, type = type)
}

vcov_fixest_cluster <- function(fit,
                                type = c("naive", "CR0", "CR2", "CR3"),
                                cluster) {
  type <- match.arg(type)
  
  if (type == "naive") return(stats::vcov(fit))
  if (missing(cluster)) stop("Provide 'cluster' for clustered vcov.")
  
  clubSandwich::vcovCR(fit, cluster = cluster, type = type)
}

# ---- Pairwise p-value helpers --------------------------------------------

# ---- GEE ------------------------------------------------------------------

fit_gee_pairwise_p <- function(
    dat_at_r,
    corstr    = c("exchangeable", "independence"),
    cluster   = c("subject", "imageID"),
    vcov_type = c("naive", "CR0", "CR2", "CR3"),
    ref_level = "Group1"
) {
  corstr    <- match.arg(corstr)
  cluster   <- match.arg(cluster)
  vcov_type <- match.arg(vcov_type)
  
  dat <- droplevels(dat_at_r)
  dat$condition <- stats::relevel(dat$condition, ref = ref_level)
  
  fit <- geepack::geeglm(
    n_target ~ 0 + condition + offset(log(density)),
    id     = dat[[cluster]],
    family = poisson("log"),
    corstr = corstr,
    data   = dat
  )
  
  V <- vcov_geeglm_type(fit, type = vcov_type, cluster = dat[[cluster]])
  
  if (vcov_type == "naive") {
    p <- summary(
      emmeans::emmeans(fit, pairwise ~ condition, vcov. = V)
    )$contrasts$p.value[1]
    return(as.numeric(p))
  }
  
  L <- matrix(c(-1, 1), nrow = 1)
  wald_p(fit, L = L, V = V)
}

# ---- LM -------------------------------------------------------------------

fit_lm_pairwise_p <- function(
    dat,
    cluster   = c("subject", "imageID"),
    vcov_type = c("naive", "CR0", "CR2", "CR3"),
    ref_level = "Group1"
) {
  cluster   <- match.arg(cluster)
  vcov_type <- match.arg(vcov_type)
  
  dat <- droplevels(dat)
  dat$condition <- stats::relevel(dat$condition, ref = ref_level)
  
  fit <- lm(
    n_target ~ 0 + condition + offset(log(density)),
    data = dat
  )
  
  V <- vcov_lmglm_cluster(fit, type = vcov_type, cluster = dat[[cluster]])
  
  if (vcov_type == "naive") {
    p <- summary(
      emmeans::emmeans(fit, pairwise ~ condition, vcov. = V)
    )$contrasts$p.value[1]
    return(as.numeric(p))
  }
  
  L <- matrix(c(-1, 1), nrow = 1)
  wald_p(fit, L = L, V = V)
}

# ---- GLM ------------------------------------------------------------------

fit_glm_pairwise_p <- function(
    dat,
    cluster   = c("subject", "imageID"),
    vcov_type = c("naive", "CR0", "CR2", "CR3"),
    ref_level = "Group1"
) {
  cluster   <- match.arg(cluster)
  vcov_type <- match.arg(vcov_type)
  
  dat <- droplevels(dat)
  dat$condition <- stats::relevel(dat$condition, ref = ref_level)
  
  fit <- glm(
    n_target ~ 0 + condition + offset(log(density)),
    family = poisson("log"),
    data   = dat
  )
  
  V <- vcov_lmglm_cluster(fit, type = vcov_type, cluster = dat[[cluster]])
  
  if (vcov_type == "naive") {
    p <- summary(
      emmeans::emmeans(fit, pairwise ~ condition, vcov. = V)
    )$contrasts$p.value[1]
    return(as.numeric(p))
  }
  
  L <- matrix(c(-1, 1), nrow = 1)
  wald_p(fit, L = L, V = V)
}

# ---- FEOLS ----------------------------------------------------------------

fit_feols_pairwise_p <- function(
    dat,
    cluster   = c("subject", "imageID"),
    vcov_type = c("naive", "CR0", "CR2", "CR3"),
    ref_level = "Group1"
) {
  cluster   <- match.arg(cluster)
  vcov_type <- match.arg(vcov_type)
  
  dat <- droplevels(dat)
  dat$condition <- stats::relevel(dat$condition, ref = ref_level)
  
  fit <- fixest::feols(
    n_target ~ 0 + condition,
    offset = log(dat$density),
    data   = dat
  )
  
  V <- vcov_fixest_cluster(fit, type = vcov_type, cluster = dat[[cluster]])
  
  beta <- stats::coef(fit)
  if (length(beta) != 2) {
    stop("Expected exactly 2 condition coefficients, got: ", paste(names(beta), collapse = ", "))
  }
  
  L <- matrix(c(-1, 1), nrow = 1)
  
  if (vcov_type == "naive") {
    est <- as.numeric(L %*% beta)
    se  <- sqrt(as.numeric(L %*% V %*% t(L)))
    tval <- est / se
    df   <- fixest::degrees_freedom(fit, "t")
    return(as.numeric(2 * stats::pt(-abs(tval), df = df)))
  }
  
  wald_p(fit, L = L, V = V)
}

# ---- FEGLM ----------------------------------------------------------------

fit_feglm_pairwise_p <- function(
    dat,
    cluster   = c("subject", "imageID"),
    vcov_type = c("naive", "CR0", "CR2", "CR3"),
    ref_level = "Group1"
) {
  cluster   <- match.arg(cluster)
  vcov_type <- match.arg(vcov_type)
  
  dat <- droplevels(dat)
  dat$condition <- stats::relevel(dat$condition, ref = ref_level)
  
  fit <- fixest::feglm(
    n_target ~ 0 + condition,
    offset = log(dat$density),
    family = "poisson",
    data   = dat
  )
  
  V <- vcov_fixest_cluster(fit, type = vcov_type, cluster = dat[[cluster]])
  
  beta <- stats::coef(fit)
  if (length(beta) != 2) {
    stop("Expected exactly 2 condition coefficients, got: ", paste(names(beta), collapse = ", "))
  }
  
  L <- matrix(c(-1, 1), nrow = 1)
  
  if (vcov_type == "naive") {
    est <- as.numeric(L %*% beta)
    se  <- sqrt(as.numeric(L %*% V %*% t(L)))
    zval <- est / se
    return(as.numeric(2 * stats::pnorm(-abs(zval))))
  }
  
  wald_p(fit, L = L, V = V)
}

# ---- Evaluate all models per radius ---------------------------------------

eval_fixed_r_all_models <- function(
    spe_obj,
    Rs,
    cluster   = c("subject", "imageID"),
    win       = "convex",
    lambda    = NA,
    rep_id    = NA,
    debug_dir = NULL
) {
  cluster <- match.arg(cluster)
  vcov_types <- c("naive", "CR0", "CR2", "CR3")
  
  purrr::map_dfr(Rs, function(r) {
    t_radius0 <- Sys.time()
    
    dat_r <- model_data_gen(spe_obj, ref = REF, target = TARGET, rad = r, window = win)
    
    # diagnostics / skipping
    if (!is.data.frame(dat_r) || nrow(dat_r) == 0L) {
      diagnose_dat_r(dat_r, r, lambda, rep_id, "Skipped: empty dataset.", debug_dir)
      return(tibble(engine=NA_character_, vcov_type=NA_character_, method=NA_character_, r=r,
                    p=NA_real_, score=NA_real_, elapsed_method_sec=NA_real_,
                    elapsed_radius_sec=as.numeric(difftime(Sys.time(), t_radius0, units="secs"))))
    }
    if (dplyr::n_distinct(dat_r$condition) < 2L) {
      diagnose_dat_r(dat_r, r, lambda, rep_id, "Skipped: only one condition present.", debug_dir)
      return(tibble(engine=NA_character_, vcov_type=NA_character_, method=NA_character_, r=r,
                    p=NA_real_, score=NA_real_, elapsed_method_sec=NA_real_,
                    elapsed_radius_sec=as.numeric(difftime(Sys.time(), t_radius0, units="secs"))))
    }
    if (dplyr::n_distinct(dat_r[[cluster]]) < 2L) {
      diagnose_dat_r(dat_r, r, lambda, rep_id,
                     sprintf("Skipped: only one cluster (%s).", cluster), debug_dir)
      return(tibble(engine=NA_character_, vcov_type=NA_character_, method=NA_character_, r=r,
                    p=NA_real_, score=NA_real_, elapsed_method_sec=NA_real_,
                    elapsed_radius_sec=as.numeric(difftime(Sys.time(), t_radius0, units="secs"))))
    }
    if (any(!is.finite(dat_r$density)) || any(dat_r$density <= 0)) {
      diagnose_dat_r(dat_r, r, lambda, rep_id, "Skipped: invalid densities.", debug_dir)
      return(tibble(engine=NA_character_, vcov_type=NA_character_, method=NA_character_, r=r,
                    p=NA_real_, score=NA_real_, elapsed_method_sec=NA_real_,
                    elapsed_radius_sec=as.numeric(difftime(Sys.time(), t_radius0, units="secs"))))
    }
    
    res_list <- list()
    
    timed_fit <- function(engine, vt, fit_fun) {
      t0 <- Sys.time()
      out <- tryCatch({
        pval <- fit_fun()
        tibble(
          engine = engine, vcov_type = vt, r = r,
          p = pval,
          score = ifelse(is.na(pval), NA_real_, -log10(pval)),
          elapsed_method_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
        )
      }, error = function(e) {
        diagnose_dat_r(dat_r, r, lambda, rep_id,
                       paste(engine, vt, "error:", conditionMessage(e)), debug_dir)
        tibble(engine=engine, vcov_type=vt, r=r,
               p=NA_real_, score=NA_real_,
               elapsed_method_sec=as.numeric(difftime(Sys.time(), t0, units="secs")))
      })
      out
    }
    
    # LM
    for (vt in vcov_types) {
      res_list[[paste0("LM_", vt)]] <- timed_fit("LM", vt, function() {
        fit_lm_pairwise_p(dat_r, cluster = cluster, vcov_type = vt)
      })
    }
    
    # GLM
    for (vt in vcov_types) {
      res_list[[paste0("GLM_", vt)]] <- timed_fit("GLM", vt, function() {
        fit_glm_pairwise_p(dat_r, cluster = cluster, vcov_type = vt)
      })
    }
    
    # FEOLS
    for (vt in vcov_types) {
      res_list[[paste0("FEOLS_", vt)]] <- timed_fit("FEOLS", vt, function() {
        fit_feols_pairwise_p(dat_r, cluster = cluster, vcov_type = vt)
      })
    }
    
    # FEGLM
    for (vt in vcov_types) {
      res_list[[paste0("FEGLM_", vt)]] <- timed_fit("FEGLM", vt, function() {
        fit_feglm_pairwise_p(dat_r, cluster = cluster, vcov_type = vt)
      })
    }
    
    # GEE
    for (vt in vcov_types) {
      res_list[[paste0("GEE_", vt)]] <- timed_fit("GEE", vt, function() {
        fit_gee_pairwise_p(dat_r, corstr = "independence", cluster = cluster, vcov_type = vt)
      })
    }
    
    bind_rows(res_list) %>%
      mutate(
        method = paste(engine, vcov_type, paste0("r", r), sep = "_"),
        elapsed_radius_sec = as.numeric(difftime(Sys.time(), t_radius0, units = "secs"))
      )
  })
}

# ---- One replicate runner -------------------------------------------------

run_one_rep_combined <- function(i, lambda) {
  pr <- simulate_pair(i, lambda, null_type = "null_equal_coloc")
  
  all_sig <- eval_fixed_r_all_models(
    pr$sig$spe, Rs_fixed,
    cluster   = "subject",
    lambda    = lambda,
    rep_id    = i,
    debug_dir = "/Users/shreyarao/Desktop/"
  ) %>% mutate(scenario = "Signal")
  
  all_null <- eval_fixed_r_all_models(
    pr$null$spe, Rs_fixed,
    cluster   = "subject",
    lambda    = lambda,
    rep_id    = i,
    debug_dir = "/Users/shreyarao/Desktop/"
  ) %>% mutate(scenario = "NoSignal")
  
  bind_rows(all_sig, all_null) %>%
    mutate(lambda = lambda, replicate = i, cs_test = CS_TEST)
}

# ---- Outer loop over lambdas & replicates ---------------------------------

results_all_list <- list()

handlers("txtprogressbar")
plan(multisession, workers = 8)

for (lam in lambdas) {
  message(sprintf("\nλ = %d  --- running %d replicates (LM/GLM/FEOLS/FEGLM/GEE - Multiple Image) [test=%s]",
                  lam, nSim, CS_TEST))
  with_progress({
    p <- progressor(along = seq_len(nSim))
    res_lam <- future_lapply(seq_len(nSim), function(b) {
      p(sprintf("λ=%d | rep=%d", lam, b))
      run_one_rep_combined(i = b, lambda = lam)
    }, future.seed = TRUE)
    
    # Combine and store for this lambda
    results_all_list[[as.character(lam)]] <- bind_rows(res_lam)
    results_lam_df <- results_all_list[[as.character(lam)]]
    
    # ---- save results for this λ immediately (EDT-tagged) ----
    saveRDS(
      results_lam_df,
      file = sprintf("/Users/shreyarao/Desktop/results_multimodel_multi_EDT_lambda%03d.rds", lam)
    )
    
    saveRDS(
      results_lam_df %>%
        distinct(lambda, replicate, engine, vcov_type, method, r, scenario, elapsed_method_sec, cs_test),
      file = sprintf("/Users/shreyarao/Desktop/timings_multimodel_multi_EDT_lambda%03d.rds", lam)
    )
    
    # Free memory
    rm(res_lam, results_lam_df)
    gc()
  })
}

# Optionally combine all into one big RDS at the end
results_all <- bind_rows(results_all_list) %>%
  mutate(label = ifelse(scenario == "Signal", 1L, 0L))

saveRDS(results_all, "/Users/shreyarao/Desktop/results_multimodel_multi_EDT_all.rds")

timings <- results_all %>%
  distinct(lambda, replicate, engine, vcov_type, method, r, scenario, elapsed_method_sec, cs_test)

saveRDS(timings, "/Users/shreyarao/Desktop/timings_multimodel_multi_EDT_all.rds")
