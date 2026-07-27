patient_total <- function(patient) {
  sum(patient$n_ij * patient$mu_hat)
}

compute_group_totals <- function(patients) {
  totals <- sapply(patients, patient_total)
  groups <- sapply(patients, function(p) p$group)
  tapply(totals, groups, sum)
}

# making reduced matrix elements
reduced_dpr1 <- function(patient) {
  Lambda_bar <- patient$mu_hat^2
  f_bar      <- sqrt(patient$n_ij) * patient$mu_hat^(3/2)
  list(Lambda_bar = Lambda_bar, f_bar = f_bar)
}

# r_i = residuals
Pi_adjoint <- function(r_i, image_index, n_ij) {
  n_i <- length(n_ij)
  image_index <- factor(image_index, levels = seq_len(n_i))
  sums <- tapply(r_i, image_index, sum)
  as.numeric(sums) / sqrt(n_ij)
}

# applying P_i
Pi_apply <- function(v, image_index, n_ij) {
  n_i <- length(n_ij)
  image_index <- factor(image_index, levels = seq_len(n_i))
  (v / sqrt(n_ij))[as.integer(image_index)]
}

# A_i bar

assemble_Abar <- function(mu_hat, eigenpairs) {
  n_i <- length(mu_hat)
  D_bar <- sqrt(mu_hat)
  
  V_bar <- sapply(eigenpairs, function(p) p$v)
  lambda_eig <- sapply(eigenpairs, function(p) p$lambda)
  
  inner <- V_bar %*% diag(1 / sqrt(lambda_eig), n_i) %*% t(V_bar)
  D_bar_diag <- diag(D_bar, n_i)
  
  D_bar_diag %*% inner %*% D_bar_diag - diag(n_i)
}

# making sure we flip signs and negate to get required results 
eig_Gbar_dpr1 <- function(Lambda_bar, f_bar, S_g) {
  pairs <- dpr1eig_full(D = -Lambda_bar, z = f_bar, rho = 1 / S_g)
  lapply(pairs, function(p) list(lambda = -p$lambda, v = p$v))
}

eig_Gbar_direct <- function(Lambda_bar, f_bar, S_g) {
  n_i <- length(Lambda_bar)
  Gbar <- diag(Lambda_bar, n_i) - (1/S_g) * (f_bar %*% t(f_bar))
  
  eig <- eigen(Gbar, symmetric = TRUE)
  
  lapply(seq_len(n_i), function(k) {
    list(lambda = eig$values[k], v = eig$vectors[, k])
  })
}

apply_Ai <- function(r_i, patient, image_index, S_g, method = c("direct", "dpr1")) {
  method <- match.arg(method)
  
  reduced <- reduced_dpr1(patient)
  eigenpairs <- if (method == "direct") {
    eig_Gbar_direct(reduced$Lambda_bar, reduced$f_bar, S_g)
  } else {
    eig_Gbar_dpr1(reduced$Lambda_bar, reduced$f_bar, S_g)
  }
  Abar <- assemble_Abar(patient$mu_hat, eigenpairs)
  
  d <- Pi_adjoint(r_i, image_index, patient$n_ij)
  correction <- Abar %*% d
  
  r_i + Pi_apply(as.numeric(correction), image_index, patient$n_ij)
}

vcovCR2_fast_multi <- function(patients, method = c("direct", "dpr1")) {
  method <- match.arg(method)
  
  S_g <- compute_group_totals(patients)
  
  M_tilde <- matrix(0, nrow = 2, ncol = 2)
  raw_sum_list <- numeric(length(patients))
  group_list <- integer(length(patients))
  
  for (idx in seq_along(patients)) {
    patient <- patients[[idx]]
    adjusted <- apply_Ai(patient$r_i, patient, patient$image_index, S_g[patient$group], method = method)
    score <- numeric(2)
    score[patient$group] <- sum(adjusted)
    M_tilde <- M_tilde + score %*% t(score)
    
    raw_sum_list[idx] <- sum(adjusted)
    group_list[idx] <- patient$group
  }
  
  B <- diag(1 / S_g)
  V <- B %*% M_tilde %*% B
  
  attr(V, "S_g") <- S_g
  attr(V, "group") <- group_list
  attr(V, "raw_sum") <- raw_sum_list
  
  V
}

build_patient <- function(cluster_id, cluster_vec, df_result, fit, condition_levels) {
  subj_rows <- cluster_vec == cluster_id
  subj_df <- df_result[subj_rows, ]
  subj_resid <- stats::residuals(fit, type = "response")[subj_rows]
  subj_fitted <- stats::fitted(fit)[subj_rows]
  
  image_order <- unique(subj_df$imageID)
  
  n_ij <- unname(sapply(image_order, function(im) sum(subj_df$imageID == im)))
  mu_hat <- unname(sapply(image_order, function(im) subj_fitted[subj_df$imageID == im][1]))
  image_index <- unname(match(subj_df$imageID, image_order))
  group <- as.integer(factor(subj_df$condition[1], levels = condition_levels))
  
  list(
    n_ij = n_ij,
    mu_hat = mu_hat,
    group = group,
    r_i = unname(subj_resid),
    image_index = image_index
  )
}

compute_e_i <- function(raw_sum, group, S_g) {
  w <- c(-1 / S_g[1], 1 / S_g[2])
  w[group] * raw_sum
}

compute_P_diag <- function(patient, S_g, method = "direct") {
  mu_i <- patient$mu_hat[patient$image_index]
  N_i <- length(mu_i)
  w <- c(-1 / S_g[1], 1 / S_g[2])
  
  A_ones <- apply_Ai(rep(1, N_i), patient, patient$image_index, S_g[patient$group], method = method)
  
  w[patient$group]^2 * sum(mu_i * A_ones^2)
}

compute_h_j <- function(patient, S_g, method = "direct") {
  mu_i <- patient$mu_hat[patient$image_index]
  
  A_mu <- apply_Ai(mu_i, patient, patient$image_index, S_g[patient$group], method = method)
  tau <- sum(A_mu)
  
  tau / S_g[patient$group]^1.5
}

compute_nu_hat <- function(patients, S_g, group, method = "direct") {
  m <- length(patients)
  
  P_diag <- sapply(patients, compute_P_diag, S_g = S_g, method = method)
  h <- sapply(patients, compute_h_j, S_g = S_g, method = method)
  
  P <- matrix(0, m, m)
  for (j in 1:m) for (k in 1:m) {
    if (group[j] == group[k]) P[j, k] <- -h[j] * h[k]
  }
  diag(P) <- P_diag - h^2
  
  Omega <- sum(diag(P))
  Omega^2 / sum(P^2)
}

waldTest_CR2_fast <- function(L_beta, V, patients, method = "direct") {
  raw_sum <- attr(V, "raw_sum")
  group <- attr(V, "group")
  S_g <- attr(V, "S_g")
  
  e_i <- compute_e_i(raw_sum, group, S_g)
  v_hat <- sum(e_i^2)
  
  nu_hat <- compute_nu_hat(patients, S_g, group, method = method)
  
  t_stat <- L_beta / sqrt(v_hat)
  p_value <- 2 * stats::pt(abs(t_stat), df = nu_hat, lower.tail = FALSE)
  
  list(t = t_stat, df = nu_hat, p.value = p_value, v_hat = v_hat)
}