## spicyFit: minimal S3 interface so closed-form fits are drop-in compatible with fast method code
new_spicyFit <- function(beta, mu_hat, y, estimator, backend) {
  structure(
    list(beta = beta, mu_hat = mu_hat, y = y,
         estimator = estimator, backend = backend),
    class = "spicyFit"
  )
}

#' @exportS3Method coef spicyFit
coef.spicyFit <- function(object, ...) object$beta

#' @exportS3Method fitted spicyFit
fitted.spicyFit <- function(object, ...) object$mu_hat

#' @exportS3Method residuals spicyFit
residuals.spicyFit <- function(object, type = "response", ...) {
  if (type != "response") stop("spicyFit only supports type = 'response'")
  object$y - object$mu_hat
}

check_rank1_design <- function(dfPair, estimator = c("mle", "firth")) {
  estimator <- match.arg(estimator)
  mm <- stats::model.matrix(~ 0 + condition, data = dfPair)
  ok <- ncol(mm) == 2 &&
    all(rowSums(mm) == 1) &&
    all(mm %in% c(0, 1))
  if (!ok) {
    fallback <- if (estimator == "firth") "backend = 'brglm2'" else "the glm() path"
    stop(
      "Closed-form fitting requires the no-covariate, rank-1-per-group design ",
      "(n ~ 0 + condition, offset = log(density)); this pair has additional ",
      "covariates or a different structure. Use ", fallback, " instead.",
      call. = FALSE
    )
  }
  invisible(mm)
}

## closed-form MLE: beta_g = log(Y_g / D_g)
fit_mle_closed_form <- function(dfPair) {
  mm <- check_rank1_design(dfPair, estimator = "mle")
  g <- apply(mm, 1, function(r) which(r == 1))
  
  Y_g <- tapply(dfPair$n, g, sum)
  D_g <- tapply(dfPair$density, g, sum)
  
  beta <- log(Y_g / D_g)
  names(beta) <- colnames(mm) 
  
  mu_hat <- exp(beta[g]) * dfPair$density
  
  new_spicyFit(beta = beta, mu_hat = unname(mu_hat), y = dfPair$n,
               estimator = "mle", backend = "closed_form")
}

## closed-form Firth: beta_g^Firth = log((Y_g + 0.5) / D_g)
fit_firth_closed_form <- function(dfPair) {
  mm <- check_rank1_design(dfPair, estimator = "firth")
  g <- apply(mm, 1, function(r) which(r == 1))
  
  Y_g <- tapply(dfPair$n, g, sum)
  D_g <- tapply(dfPair$density, g, sum)
  
  beta <- log((Y_g + 0.5) / D_g)
  names(beta) <- colnames(mm)
  
  mu_hat <- exp(beta[g]) * dfPair$density
  
  new_spicyFit(beta = beta, mu_hat = unname(mu_hat), y = dfPair$n,
               estimator = "firth", backend = "closed_form")
}

## Firth via brglm2::brglmFit -- general, no design restriction
fit_firth_brglm2 <- function(dfPair) {
  glm(n ~ 0 + condition, offset = log(density), family = poisson(),
      data = dfPair, method = "brglmFit", type = "AS_mean")
}

## ordinary MLE via glm() -- general, no design restriction
fit_mle_glm <- function(dfPair) {
  glm(n ~ 0 + condition, offset = log(density), family = poisson(),
      data = dfPair)
}

## dispatcher: routes to the right backend
fit_pair <- function(dfPair, estimator = c("mle", "firth"),
                     backend = c("closed_form", "glm", "brglm2")) {
  estimator <- match.arg(estimator)
  backend <- match.arg(backend)
  
  fit <- switch(
    paste(estimator, backend),
    "mle closed_form"   = fit_mle_closed_form(dfPair),
    "mle glm"           = fit_mle_glm(dfPair),
    "firth closed_form" = fit_firth_closed_form(dfPair),
    "firth brglm2"      = fit_firth_brglm2(dfPair),
    stop("Invalid estimator/backend combination: ", estimator, "/", backend,
         call. = FALSE)
  )
  
  list(fit = fit, estimator = estimator, backend = backend)
}