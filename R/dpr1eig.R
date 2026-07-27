dpr1eig_shift_index <- function(D, z, rho, k) {
  if (k == 1) {
    return(list(i = 1, branch = "exterior"))
  }
  
  D_bar <- D - D[k]
  tau   <- D_bar[k - 1] / 2
  F <- 1 + rho * sum(z^2 / (D_bar - tau))
  
  # deciding sigma using F
  if (F > 0) {
    list(i = k, branch = "dk")
  } else {
    list(i = k - 1, branch = "dk_minus_1")
  }
}

dpr1eig_split <- function(D, z, i) {
  n <- length(D)
  sigma <- D[i]
  
  if (i > 1) {
    D1 <- D[1:(i - 1)] - sigma
    z1 <- z[1:(i - 1)]
  } else {
    D1 <- numeric(0)
    z1 <- numeric(0)
  }
  
  if (i < n) {
    D2 <- D[(i + 1):n] - sigma
    z2 <- z[(i + 1):n]
  } else {
    D2 <- numeric(0)
    z2 <- numeric(0)
  }
  
  zeta_i <- z[i]
  
  list(sigma = sigma, D1 = D1, z1 = z1, D2 = D2, z2 = z2, zeta_i = zeta_i)
}

dpr1eig_wb <- function(split, rho) {
  D1 <- split$D1
  z1 <- split$z1
  D2 <- split$D2
  z2 <- split$z2
  zeta_i <- split$zeta_i
  
  w1 <- -(z1 / D1) / zeta_i
  w2 <- -(z2 / D2) / zeta_i
  
  b <- (1 / rho + sum(z1^2 / D1) + sum(z2^2 / D2)) / zeta_i^2
  
  list(w1 = w1, w2 = w2, b = b)
}

build_secular_terms <- function(split, wb) {
  Delta <- c(1 / split$D1, 1 / split$D2)
  w     <- c(wb$w1, wb$w2)
  list(Delta = Delta, w = w, b = wb$b)
}

bracket_root <- function(Delta, w, b, direction = c("rightmost", "leftmost"),
                         offset = 1e-8, step = 1) {
  direction <- match.arg(direction)
  h <- function(nu) b - nu - sum(w^2 / (Delta - nu))
  
  if (direction == "rightmost") {
    anchor <- max(Delta)
    near <- anchor + offset
    while (h(near) <= 0) {
      offset <- offset / 2
      near <- anchor + offset
    }
    far <- near + step
    while (h(far) > 0) {
      step <- step * 2
      far <- near + step
    }
    c(lower = near, upper = far)
    
  } else {
    anchor <- min(Delta)
    near <- anchor - offset
    while (h(near) >= 0) {
      offset <- offset / 2
      near <- anchor - offset
    }
    far <- near - step
    while (h(far) < 0) {
      step <- step * 2
      far <- near - step
    }
    c(lower = far, upper = near)
  }
}

bisect_root <- function(Delta, w, b, bracket, tol = 1e-16, max_iter = 200) {
  h <- function(nu) b - nu - sum(w^2 / (Delta - nu))
  
  lo <- unname(bracket["lower"])
  hi <- unname(bracket["upper"])
  
  for (iter in seq_len(max_iter)) {
    mid <- (lo + hi) / 2
    h_mid <- h(mid)
    
    if (abs(h_mid) < tol || (hi - lo) < tol) {
      return(mid)
    }
    
    if (h_mid > 0) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  
  mid
}

dpr1eig_assemble <- function(split, nu) {
  sigma <- split$sigma
  D1 <- split$D1
  z1 <- split$z1
  D2 <- split$D2
  z2 <- split$z2
  zeta_i <- split$zeta_i
  
  mu <- 1 / nu
  x1 <- z1 / (D1 - mu)
  x2 <- z2 / (D2 - mu)
  # puts entries in the original index order
  x  <- c(x1, -zeta_i / mu, x2)
  
  v <- x / sqrt(sum(x^2))
  lambda <- mu + sigma
  
  list(lambda = lambda, v = v)
}

# diganostics

dpr1eig_precision_diagnostics <- function(split, rho) {
  D1 <- split$D1
  z1 <- split$z1
  D2 <- split$D2
  z2 <- split$z2
  zeta_i <- split$zeta_i
  
  n <- length(D1) + length(D2) + 1
  
  term1 <- rho * sum(z1^2 / D1)
  term2 <- rho * sum(z2^2 / D2)
  
  Kb <- abs(1 + term1 - term2) / abs(1 + term1 + term2)
  
  Kz <- (sum(abs(z1)) + sum(abs(z2))) / abs(zeta_i)
  
  kappa_nu <- min(
    (n + 4) * sqrt(n) * Kb,
    3 * sqrt(n) + (n + 4) * (1 + 2 * Kz)
  )
  
  list(Kb = Kb, Kz = Kz, kappa_nu = kappa_nu)
}

dpr1eig_shift_diagnostics <- function(terms, nu, direction) {
  other_direction <- if (direction == "rightmost") "leftmost" else "rightmost"
  
  other_bracket <- bracket_root(terms$Delta, terms$w, terms$b, other_direction)
  nu_other <- bisect_root(terms$Delta, terms$w, terms$b, other_bracket)
  
  K_nu <- max(abs(nu), abs(nu_other)) / abs(nu)
  
  list(nu_other = nu_other, K_nu = K_nu)
}

dpr1eig_cancellation_flag <- function(lambda, D, k, threshold = 1e-2) {
  if (k == 1) {
    gap <- abs(lambda - D[k])
  } else {
    gap <- min(abs(lambda - D[k]), abs(lambda - D[k - 1]))
  }
  
  flagged <- abs(lambda) < threshold * gap
  list(flagged = flagged, gap = gap)
}

# mid level orchestrated function
dpr1eig <- function(D, z, rho, k, diagnostics = FALSE) {
  if (length(D) == 1) {
    lambda <- D[1] + rho * z[1]^2
    result <- list(lambda = lambda, v = 1)
    
    if (diagnostics) {
      result$diagnostics <- list(
        Kb = NA, Kz = NA, kappa_nu = NA,
        K_nu = NA, cancellation_flagged = NA
      )
    }
    
    return(result)
  }
  
  idx   <- dpr1eig_shift_index(D, z, rho, k)
  split <- dpr1eig_split(D, z, idx$i)
  wb    <- dpr1eig_wb(split, rho)
  terms <- build_secular_terms(split, wb)
  
  direction <- if (idx$branch == "dk_minus_1") "leftmost" else "rightmost"
  
  bracket <- bracket_root(terms$Delta, terms$w, terms$b, direction)
  nu      <- bisect_root(terms$Delta, terms$w, terms$b, bracket)
  
  result <- dpr1eig_assemble(split, nu)
  
  if (diagnostics) {
    precision <- dpr1eig_precision_diagnostics(split, rho)
    shift     <- dpr1eig_shift_diagnostics(terms, nu, direction)
    cancel    <- dpr1eig_cancellation_flag(result$lambda, D, k)
    
    result$diagnostics <- list(
      Kb = precision$Kb,
      Kz = precision$Kz,
      kappa_nu = precision$kappa_nu,
      K_nu = shift$K_nu,
      cancellation_flagged = cancel$flagged
    )
  }
  
  result
}

lift_reduced_eigenvector <- function(v_reduced, group_indices, z, n) {
  v_full <- numeric(n)
  for (g in seq_along(group_indices)) {
    idx <- group_indices[[g]]
    if (length(idx) == 1) {
      v_full[idx] <- v_reduced[g]
    } else {
      z_group <- z[idx]
      v_full[idx] <- v_reduced[g] * (z_group / sqrt(sum(z_group^2)))
    }
  }
  v_full
}


# if poles are tied

find_tied_groups <- function(D) {
  split(seq_along(D), D)
}

generate_tied_eigenpairs <- function(D, z, indices) {
  d <- D[indices[1]]
  z_group <- z[indices]
  m <- length(indices)
  
  Q <- qr.Q(qr(matrix(z_group, ncol = 1)), complete = TRUE)
  basis <- Q[, 2:m, drop = FALSE]
  
  n <- length(D)
  eigenpairs <- vector("list", m - 1)
  for (j in seq_len(m - 1)) {
    v_full <- numeric(n)
    v_full[indices] <- basis[, j]
    eigenpairs[[j]] <- list(lambda = d, v = v_full)
  }
  eigenpairs
}

collapse_tied_poles <- function(D, z, tied) {
  n_groups <- length(tied)
  D_reduced <- numeric(n_groups)
  z_reduced <- numeric(n_groups)
  group_indices <- vector("list", n_groups)
  
  for (g in seq_len(n_groups)) {
    idx <- tied[[g]]
    D_reduced[g] <- D[idx[1]]
    z_reduced[g] <- sqrt(sum(z[idx]^2))
    group_indices[[g]] <- idx
  }
  
  order_idx <- order(D_reduced, decreasing = TRUE)
  
  list(
    D = D_reduced[order_idx],
    z = z_reduced[order_idx],
    group_indices = group_indices[order_idx]
  )
}

## top level orchestrator to output multiple eigenpairs
dpr1eig_full <- function(D, z, rho) {
  n <- length(D)
  tied <- find_tied_groups(D)
  
  closed_form <- list()
  for (g in tied) {
    if (length(g) > 1) {
      closed_form <- c(closed_form, generate_tied_eigenpairs(D, z, g))
    }
  }
  
  reduced <- collapse_tied_poles(D, z, tied)
  n_reduced <- length(reduced$D)
  
  reduced_pairs <- vector("list", n_reduced)
  for (k in seq_len(n_reduced)) {
    res <- dpr1eig(reduced$D, reduced$z, rho, k)
    v_full <- lift_reduced_eigenvector(res$v, reduced$group_indices, z, n)
    reduced_pairs[[k]] <- list(lambda = res$lambda, v = v_full)
  }
  
  all_pairs <- c(closed_form, reduced_pairs)
  lambdas <- sapply(all_pairs, function(p) p$lambda)
  all_pairs[order(lambdas, decreasing = TRUE)]
}