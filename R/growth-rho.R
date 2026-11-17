growth_rho_pattern <- function(pattern, occasions, factors) {
  n <- length(occasions)
  Q <- length(factors)
  if (!is.numeric(pattern) || any(is.nan(pattern)) || any(is.infinite(pattern))) {
    stop("Time scores must be numeric, with NA for scores to estimate.", call. = FALSE)
  }
  if (is.null(dim(pattern)) && length(pattern) == n) {
    pattern <- matrix(pattern, n, Q)
  }
  if (!identical(dim(pattern), c(n, Q))) {
    stop("Supply one score per occasion, or a matrix with one column per factor.", call. = FALSE)
  }
  dimnames(pattern) <- list(as.character(occasions), factors)
  free <- is.na(pattern)
  initial <- pattern
  zero_time <- integer(Q)
  names <- outer(rownames(pattern), factors, function(t, q) sprintf("rho[%s,%s]", t, q))

  for (q in seq_len(Q)) {
    zero <- which(pattern[, q] == 0)
    one <- which(pattern[, q] == 1)
    if (!length(zero) || !length(one)) {
      stop("Factor ", factors[q], " needs a fixed score of 0 and a fixed score of 1.",
           call. = FALSE)
    }
    zero_time[q] <- zero[1L]
    if (any(free[, q])) {
      t <- seq_len(n)
      anchors <- which(!free[, q])
      values <- stats::approx(anchors, pattern[anchors, q], xout = t, rule = 2)$y
      left <- t < anchors[1L]
      right <- t > tail(anchors, 1L)
      slope <- 1 / (one[1L] - zero[1L])
      values[left] <- values[left] + (t[left] - anchors[1L]) * slope
      values[right] <- values[right] + (t[right] - tail(anchors, 1L)) * slope
      initial[free[, q], q] <- values[free[, q]]
    }
  }
  list(pattern = pattern, free = free, initial = initial, names = names, zero_time = zero_time)
}

growth_rho_scores <- function(spec, par) {
  scores <- spec$pattern
  scores[spec$free] <- par[spec$index]
  scores
}

growth_rho_update_design <- function(X, Z, scores) {
  X[, , , 2L] <- rep(as.double(scores), each = dim(X)[1L])
  Z[, , , 2L] <- rep(as.double(scores), each = dim(Z)[1L])
  list(X = X, Z = Z)
}

growth_rho_native <- function(par, spec, X, Z) {
  design <- growth_rho_update_design(array(X, spec$dim_x), array(Z, spec$dim_z),
                                     growth_rho_scores(spec, par))
  keep <- setdiff(seq_along(par), spec$index)
  list(par = par[keep], X = as.double(design$X), Z = as.double(design$Z))
}
