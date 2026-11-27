#cross-sectional ordinal probit CFA
cfa_fit <- function(data, namvar = names(data), Q = NULL,
                    simple = NULL, loadmask = NULL, marker_item = NULL,
                    factor_names = NULL, methInteg = c("GH", "AGH", "QMC"),
                    nquadgh = 30L, nquadqmc = 16384L,
                    max_nodes = 200000L, max_patterns = 4096L,
                    full_patterns = TRUE, compute_derivatives = TRUE,
                    compute_hessian = FALSE, start = NULL, maxiter = 500L,
                    reltol = 1e-10, diff_step = 1e-5,
                    agh_maxiter = 100L, agh_tol = 1e-7) {
  methInteg <- match.arg(methInteg)
  check_count <- function(x, min = 1L) {
    if (length(x) != 1L || !is.finite(x) || x < min ||
        x != floor(x) || x > .Machine$integer.max)
      stop("invalid integer control")
    as.integer(x)
  }
  nquadgh <- check_count(nquadgh)
  nquadqmc <- check_count(nquadqmc)
  max_nodes <- check_count(max_nodes)
  max_patterns <- check_count(max_patterns)
  maxiter <- check_count(maxiter, 0L)
  agh_maxiter <- check_count(agh_maxiter)
  if (length(agh_tol) != 1L || !is.finite(agh_tol) || agh_tol <= 0)
    stop("invalid AGH tolerance")
  if (!is.finite(diff_step) || diff_step <= 0 ||
      !is.finite(reltol) || reltol <= 0) stop("invalid tolerance")
  data <- as.data.frame(data)
  if (!length(namvar) || anyDuplicated(namvar) ||
      !all(namvar %in% names(data)) || !nrow(data))
    stop("invalid item data")
  K <- length(namvar)
  N <- nrow(data)
  levels_item <- vector("list", K)
  Y <- matrix(NA_integer_, N, K, dimnames = list(NULL, namvar))
  for (k in seq_len(K)) {
    x <- data[[namvar[k]]]
    if (is.ordered(x)) {
      lev <- levels(x)
    } else if (is.numeric(x) || is.logical(x)) {
      x <- as.numeric(x)
      lev <- sort(unique(x[!is.na(x)]))
      if (length(lev) < 2L || !all(is.finite(lev)) ||
          !lev[1L] %in% c(0, 1) ||
          !all(diff(lev) == 1))
        stop("numeric items must have consecutive integer categories")
    } else {
      stop("use numeric responses or ordered factors")
    }
    if (length(lev) < 2L) stop("each item needs at least two categories")
    Y[, k] <- match(x, lev) - 1L
    if (length(unique(stats::na.omit(Y[, k]))) != length(lev))
      stop("each declared category must be observed")
    levels_item[[k]] <- lev
  }
  names(levels_item) <- namvar
  C <- lengths(levels_item)
  if (any(rowSums(!is.na(Y)) == 0L))
    warning("entirely missing subjects contribute zero log likelihood and score")

  #loading structure and marker identification
  if (!is.null(loadmask) && !is.null(simple))
    stop("supply either simple or loadmask")
  if (is.null(loadmask)) {
    if (is.null(simple)) {
      if (!is.null(Q) && Q != 1L)
        stop("for multiple factors supply simple or loadmask")
      simple <- rep(1L, K)
    }
    if (length(simple) != K || anyNA(simple) ||
        any(!is.finite(simple) | simple < 1 | simple != floor(simple)))
      stop("simple must contain one positive factor index per item")
    if (is.null(Q)) Q <- max(simple)
    Q <- check_count(Q)
    if (any(simple > Q)) stop("simple exceeds Q")
    loadmask <- outer(simple, seq_len(Q), "==")
  } else {
    loadmask <- as.matrix(loadmask)
    if (nrow(loadmask) != K || !ncol(loadmask) ||
        anyNA(loadmask) || any(!loadmask %in% c(0, 1)))
      stop("loadmask must be a K by Q zero-one matrix")
    if (!is.null(Q) && Q != ncol(loadmask)) stop("Q and loadmask disagree")
    Q <- ncol(loadmask)
    loadmask <- loadmask == 1
  }
  if (any(rowSums(loadmask) == 0L) || any(colSums(loadmask) == 0L))
    stop("every item and factor must have an active loading")
  if (is.null(factor_names)) factor_names <- paste0("F", seq_len(Q))
  if (length(factor_names) != Q || anyNA(factor_names) ||
      anyDuplicated(factor_names)) stop("invalid factor names")
  dimnames(loadmask) <- list(namvar, factor_names)
  if (is.null(marker_item)) {
    marker <- vapply(seq_len(Q), function(q) {
      z <- which(loadmask[, q] & rowSums(loadmask) == 1L)
      if (!length(z)) stop("each factor requires a pure marker item")
      z[1L]
    }, integer(1))
  } else {
    marker <- if (is.character(marker_item)) match(marker_item, namvar)
              else marker_item
    if (length(marker) != Q || anyNA(marker) ||
        any(!is.finite(marker) | marker < 1 | marker > K |
            marker != floor(marker)) || anyDuplicated(marker))
      stop("invalid marker items")
    marker <- as.integer(marker)
  }
  if (any(!loadmask[cbind(marker, seq_len(Q))]) ||
      any(rowSums(loadmask)[marker] != 1L))
    stop("marker items must load on their own factor only")
  free_load <- loadmask
  free_load[cbind(marker, seq_len(Q))] <- FALSE
  li <- which(free_load, arr.ind = TRUE)
  ni <- setdiff(seq_len(K), marker)
  ci <- which(lower.tri(matrix(0, Q, Q), diag = TRUE), arr.ind = TRUE)
  ti <- do.call(rbind, lapply(seq_len(K), function(k) {
    j <- if (C[k] > 2L) seq.int(2L, C[k] - 1L) else integer()
    cbind(rep.int(k, length(j)), j)
  }))
  colnames(ti) <- c("item", "threshold")
  nl <- nrow(li); nn <- length(ni); nc <- nrow(ci); nt <- nrow(ti)
  block <- function(n, offset) offset + seq_len(n)
  il <- block(nl, 0L); iv <- block(nn, nl)
  im <- block(Q, nl + nn); ic <- block(nc, nl + nn + Q)
  it <- block(nt, nl + nn + Q + nc)
  d <- nl + nn + Q + nc + nt
  indexed_names <- function(prefix, left, right) {
    if (!length(left)) return(character())
    paste0(prefix, "[", left, ",", right, "]")
  }
  par_names <- c(
    indexed_names("lambda", namvar[li[, 1]], factor_names[li[, 2]]),
    if (nn) paste0("nu[", namvar[ni], "]") else character(),
    paste0("mu[", factor_names, "]"),
    indexed_names("phi", factor_names[ci[, 1]], factor_names[ci[, 2]]),
    indexed_names("tau", namvar[ti[, 1]], ti[, 2])
  )
  if (d > prod(C) - 1) stop("more parameters than response-pattern degrees of freedom")
  diag_ci <- ci[, 1] == ci[, 2]

  #natural and unconstrained parameter scales
  decode <- function(theta, internal = FALSE) {
    if (length(theta) != d || any(!is.finite(theta)))
      stop("invalid parameter vector")
    lambda <- matrix(0, K, Q, dimnames = dimnames(loadmask))
    lambda[cbind(marker, seq_len(Q))] <- 1
    lambda[li] <- theta[il]
    nu <- setNames(numeric(K), namvar)
    nu[ni] <- theta[iv]
    mu <- setNames(theta[im], factor_names)
    Phi <- matrix(0, Q, Q, dimnames = list(factor_names, factor_names))
    if (internal) {
      ch <- theta[ic]
      ch[diag_ci] <- exp(ch[diag_ci])
      Phi[ci] <- ch
      Phi <- tcrossprod(Phi)
    } else {
      Phi[ci] <- theta[ic]
      Phi <- Phi + t(Phi) - diag(diag(Phi), Q, Q)
      dimnames(Phi) <- list(factor_names, factor_names)
    }
    if (any(!is.finite(Phi))) stop("nonfinite factor covariance")
    R <- chol(Phi)
    tau <- lapply(seq_len(K), function(k) {
      vals <- theta[it][ti[, 1] == k]
      if (internal) vals <- cumsum(exp(vals))
      if (length(vals) && (any(!is.finite(vals)) ||
                            any(diff(c(0, vals)) <= 0)))
        stop("thresholds must be strictly ordered above zero")
      c(-Inf, 0, vals, Inf)
    })
    names(tau) <- namvar
    list(lambda = lambda, nu = nu, mu = mu, Phi = Phi,
         tau = tau, residual_var = setNames(rep(1, K), namvar), R = R)
  }
  to_natural <- function(par) {
    z <- decode(par, internal = TRUE)
    setNames(c(z$lambda[li], z$nu[ni], z$mu, z$Phi[ci],
               unlist(lapply(z$tau, function(t) {
                 if (length(t) > 3L) t[seq.int(3L, length(t) - 1L)]
                 else numeric()
               }), use.names = FALSE)), par_names)
  }
  to_internal <- function(theta) {
    z <- decode(theta)
    par <- theta
    ch <- t(z$R)[ci]
    ch[diag_ci] <- log(ch[diag_ci])
    par[ic] <- ch
    for (k in seq_len(K)) {
      ii <- it[ti[, 1] == k]
      if (length(ii)) par[ii] <- log(diff(c(0, theta[ii])))
    }
    setNames(par, par_names)
  }
  if (is.null(start)) {
    lambda <- matrix(0, K, Q)
    lambda[loadmask] <- 0.7
    lambda[cbind(marker, seq_len(Q))] <- 1
    sd_y <- sqrt(1 + rowSums(lambda^2))
    p1 <- colMeans(Y > 0L, na.rm = TRUE)
    target <- stats::qnorm(p1) * sd_y
    mu <- target[marker]
    nu <- target - drop(lambda %*% mu)
    threshold <- if (nt) vapply(seq_len(nt), function(j) {
      k <- ti[j, 1]; t <- ti[j, 2]
      target[k] + sd_y[k] * stats::qnorm(mean(Y[, k] < t, na.rm = TRUE))
    }, numeric(1)) else numeric()
    start <- setNames(c(lambda[li], nu[ni], mu, diag(Q)[ci], threshold),
                      par_names)
  } else {
    if (!is.null(names(start))) {
      if (anyDuplicated(names(start)) || !setequal(names(start), par_names))
        stop("named start must match the natural parameter names")
      start <- start[par_names]
    }
    if (length(start) != d) stop("start has the wrong length")
  }
  par0 <- to_internal(start)

  #base integration rule
  if (methInteg %in% c("GH", "AGH")) {
    if (Q * log(nquadgh) > log(max_nodes))
      stop("GH/AGH exceeds max_nodes; reduce nquadgh or select QMC")
    rule <- fastGHQuad::gaussHermiteData(nquadgh)
    grid <- as.matrix(expand.grid(rep(list(seq_len(nquadgh)), Q)))
    nodes <- matrix(sqrt(2) * rule$x[grid], nrow(grid), Q)
    logw <- rowSums(matrix(log(rule$w[grid] / sqrt(pi)), nrow(grid), Q))
  } else {
    if (!requireNamespace("spacefillr", quietly = TRUE))
      stop("install spacefillr for QMC")
    if (nquadqmc > max_nodes) stop("QMC exceeds max_nodes")
    u <- spacefillr::generate_sobol_owen_set(nquadqmc, Q)
    u <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
    nodes <- matrix(stats::qnorm(u), nquadqmc, Q)
    logw <- rep(-log(nquadqmc), nquadqmc)
  }
  logw <- logw - (max(logw) + log(sum(exp(logw - max(logw)))))
  logdiff <- function(a, b) a + log(-expm1(pmin(b - a, 0)))
  logsumexp <- function(x) {
    m <- max(x)
    if (!is.finite(m)) return(m)
    m + log(sum(exp(x - m)))
  }
  validate_patterns <- function(pat) {
    pat <- as.matrix(pat)
    if (ncol(pat) != K || !is.numeric(pat)) stop("invalid pattern matrix")
    for (k in seq_len(K)) {
      v <- pat[, k]
      if (any(!is.na(v) & (!is.finite(v) | v < 0 | v >= C[k] | v != floor(v))))
        stop("patterns must use zero-based category codes")
    }
    pat
  }
  #stable ordinal category log probability
  category_logprob <- function(lo, hi) {
    if (all(lo == -Inf)) return(stats::pnorm(hi, log.p = TRUE))
    if (all(hi == Inf))
      return(stats::pnorm(lo, lower.tail = FALSE, log.p = TRUE))
    pos <- lo > 0
    out <- numeric(length(lo))
    out[!pos] <- logdiff(stats::pnorm(hi[!pos], log.p = TRUE),
                         stats::pnorm(lo[!pos], log.p = TRUE))
    out[pos] <- logdiff(stats::pnorm(lo[pos], lower.tail = FALSE, log.p = TRUE),
                        stats::pnorm(hi[pos], lower.tail = FALSE, log.p = TRUE))
    out
  }

  #adaptive quadrature in standard normal factor coordinates
  log_prob_agh <- function(z, pat) {
    base <- z$nu + drop(z$lambda %*% z$mu)
    B <- z$lambda %*% t(z$R)
    node_norm <- rowSums(nodes^2) / 2
    Iq <- diag(Q)
    vapply(seq_len(nrow(pat)), function(i) {
      obs <- which(!is.na(pat[i, ]))
      if (!length(obs)) return(0)
      agh_mode_calls <<- agh_mode_calls + 1
      b <- B[obs, , drop = FALSE]
      a <- base[obs]
      lo0 <- vapply(obs, function(k) z$tau[[k]][pat[i, k] + 1L], numeric(1))
      hi0 <- vapply(obs, function(k) z$tau[[k]][pat[i, k] + 2L], numeric(1))

      posterior <- function(x, deriv = FALSE) {
        m <- a + drop(b %*% x)
        lo <- lo0 - m
        hi <- hi0 - m
        lp <- vapply(seq_along(obs), function(j)
          category_logprob(lo[j], hi[j]), numeric(1))
        value <- sum(x^2) / 2 - sum(lp)
        if (!deriv) return(value)
        if (!is.finite(value)) stop("nonfinite AGH posterior")
        rl <- exp(stats::dnorm(lo, log = TRUE) - lp)
        rh <- exp(stats::dnorm(hi, log = TRUE) - lp)
        g <- rl - rh
        grad <- x - drop(crossprod(b, g))
        #curvature of the negative conditional log probability
        zl <- zh <- numeric(length(obs))
        jl <- is.finite(lo); jh <- is.finite(hi)
        zl[jl] <- lo[jl] * rl[jl]
        zh[jh] <- hi[jh] * rh[jh]
        curv <- g^2 - zl + zh
        H <- Iq + crossprod(b, b * curv)
        list(value = value, gradient = grad, hessian = H)
      }

      fn <- function(x) posterior(x)
      gr <- function(x) posterior(x, TRUE)$gradient
      #deterministic mode search and Newton refinement
      mode_fit <- stats::optim(rep(0, Q), fn, gr, method = "BFGS",
        control = list(maxit = agh_maxiter, reltol = 1e-13))
      mode <- mode_fit$par
      for (iter in seq_len(10L)) {
        p <- posterior(mode, TRUE)
        if (max(abs(p$gradient)) <= agh_tol * (1 + max(abs(mode)))) break
        direction <- tryCatch(solve(p$hessian, p$gradient),
                              error = function(e) NULL)
        if (is.null(direction) || any(!is.finite(direction))) break
        step <- 1
        accepted <- FALSE
        for (j in seq_len(20L)) {
          trial <- mode - step * direction
          ftrial <- fn(trial)
          if (is.finite(ftrial) && ftrial <= p$value + 1e-13 * (1 + abs(p$value))) {
            mode <- trial
            accepted <- TRUE
            break
          }
          step <- step / 2
        }
        if (!accepted) break
      }
      p <- posterior(mode, TRUE)
      if (any(!is.finite(p$gradient)) ||
          max(abs(p$gradient)) > agh_tol * (1 + max(abs(mode))))
        stop("AGH posterior mode did not converge")
      H <- (p$hessian + t(p$hessian)) / 2
      R <- tryCatch(chol(H), error = function(e) NULL)
      if (is.null(R)) {
        H <- stats::optimHess(mode, fn, gr,
                              control = list(ndeps = rep(1e-5, Q)))
        R <- chol((H + t(H)) / 2)
      }
      #if H = R'R, A'A = H^{-1} for A = R^{-T}
      A <- t(backsolve(R, Iq))
      adapted <- sweep(nodes %*% A, 2L, mode, "+")
      m <- sweep(adapted %*% t(b), 2L, a, "+")
      log_cond <- numeric(nrow(nodes))
      for (j in seq_along(obs)) {
        log_cond <- log_cond +
          category_logprob(lo0[j] - m[, j], hi0[j] - m[, j])
      }
      logsumexp(logw + log_cond - rowSums(adapted^2) / 2 +
                  node_norm) - sum(log(diag(R)))
    }, numeric(1))
  }

  log_prob <- function(theta, pat, internal = FALSE) {
    z <- decode(theta, internal)
    if (methInteg == "AGH") return(log_prob_agh(z, pat))
    eta <- sweep(nodes %*% z$R, 2L, z$mu, "+")
    m <- sweep(eta %*% t(z$lambda), 2L, z$nu, "+")
    if (any(!is.finite(m))) stop("nonfinite linear predictors")
    lp <- lapply(seq_len(K), function(k) {
      out <- matrix(0, nrow(nodes), C[k])
      for (j in seq_len(C[k])) {
        lo <- z$tau[[k]][j] - m[, k]
        hi <- z$tau[[k]][j + 1L] - m[, k]
        if (j == 1L) {
          out[, j] <- stats::pnorm(hi, log.p = TRUE)
        } else if (j == C[k]) {
          out[, j] <- stats::pnorm(lo, lower.tail = FALSE, log.p = TRUE)
        } else {
          pos <- lo > 0
          out[!pos, j] <- logdiff(
            stats::pnorm(hi[!pos], log.p = TRUE),
            stats::pnorm(lo[!pos], log.p = TRUE))
          out[pos, j] <- logdiff(
            stats::pnorm(lo[pos], lower.tail = FALSE, log.p = TRUE),
            stats::pnorm(hi[pos], lower.tail = FALSE, log.p = TRUE))
        }
      }
      out
    })
    vapply(seq_len(nrow(pat)), function(i) {
      v <- logw
      for (k in which(!is.na(pat[i, ])))
        v <- v + lp[[k]][, pat[i, k] + 1L]
      logsumexp(v)
    }, numeric(1))
  }
  pattern_key <- function(pat) apply(pat, 1L, paste, collapse = ":")
  key <- pattern_key(Y)
  obs_patterns <- Y[!duplicated(key), , drop = FALSE]
  id <- match(key, key[!duplicated(key)])
  counts <- tabulate(id, nrow(obs_patterns))
  #AGH fits each subject separately, including repeated response patterns
  individual_agh <- methInteg == "AGH"
  likelihood_data <- if (individual_agh) Y else obs_patterns
  likelihood_weights <- if (individual_agh) rep(1, N) else counts
  agh_mode_calls <- 0
  objective <- function(par) {
    val <- tryCatch(-sum(likelihood_weights *
                          log_prob(par, likelihood_data, TRUE)),
                    error = function(e) Inf)
    if (!is.finite(val)) 1e100 else val
  }
  if (objective(par0) >= 1e100) stop("invalid starting likelihood")
  opt <- if (maxiter == 0L) {
    list(par = par0, value = objective(par0), convergence = NA_integer_,
         message = "likelihood evaluation only")
  } else {
    stats::optim(par0, objective, method = "BFGS",
                 control = list(maxit = maxiter, reltol = reltol))
  }
  theta <- to_natural(opt$par)
  decoded <- decode(theta)
  if (maxiter > 0L && opt$convergence != 0L)
    warning("optimization did not converge")
  mode_calls_before <- agh_mode_calls
  fitted_logpi <- log_prob(theta, likelihood_data)
  final_adaptations <- agh_mode_calls - mode_calls_before
  individual_loglik <- if (individual_agh) fitted_logpi else fitted_logpi[id]
  obs_logpi <- if (individual_agh) fitted_logpi[!duplicated(key)] else fitted_logpi

  #derivatives on the natural parameter scale
  differentiate <- function(fun, x) {
    f0 <- fun(x)
    J <- matrix(NA_real_, length(f0), length(x),
                dimnames = list(NULL, par_names))
    for (j in seq_along(x)) {
      h <- diff_step * max(1, abs(x[j]))
      for (attempt in seq_len(18L)) {
        xp <- xm <- x
        xp[j] <- xp[j] + h; xm[j] <- xm[j] - h
        fp <- tryCatch(fun(xp), error = function(e) NA_real_)
        fm <- tryCatch(fun(xm), error = function(e) NA_real_)
        if (length(fp) == length(f0) && length(fm) == length(f0) &&
            all(is.finite(fp)) && all(is.finite(fm))) {
          J[, j] <- (fp - fm) / (2 * h)
          break
        }
        h <- h / 2
      }
    }
    if (anyNA(J)) stop("numerical derivative failed; check parameter boundaries")
    J
  }
  probability_function <- function(theta, patterns) {
    exp(log_prob(theta, validate_patterns(patterns)))
  }
  jacobian_function <- function(theta, patterns) {
    pat <- validate_patterns(patterns)
    lp <- log_prob(theta, pat)
    differentiate(function(t) log_prob(t, pat), theta) * exp(lp)
  }
  pat <- pi_hat <- J <- U <- U_obs <- A <- A_obs <- H <- NULL
  if (isTRUE(full_patterns)) {
    if (sum(log(C)) > log(max_patterns) + 1e-12) {
      warning("full response-pattern enumeration skipped; max_patterns exceeded")
    } else {
      pat <- as.matrix(expand.grid(lapply(C, function(c) seq.int(0L, c - 1L))))
      colnames(pat) <- namvar
      pi_hat <- exp(log_prob(theta, pat))
    }
  }
  if (isTRUE(compute_derivatives)) {
    U_obs <- differentiate(function(t) log_prob(t, likelihood_data), theta)
    if (!individual_agh) U_obs <- U_obs[id, , drop = FALSE]
    A_obs <- crossprod(U_obs) / N
    if (!is.null(pat)) {
      U <- differentiate(function(t) log_prob(t, pat), theta)
      J <- U * pi_hat
      A <- crossprod(U, U * pi_hat)
    }
    check_A <- if (is.null(A)) A_obs else A
    ev <- eigen(check_A, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 1e-8 * max(ev))
      warning("information is numerically rank deficient; check identification")
  }
  if (isTRUE(compute_hessian)) {
    H <- tryCatch(stats::optimHess(theta, function(t) {
      -sum(likelihood_weights * log_prob(t, likelihood_data))
    }, control = list(ndeps = rep(diff_step, d))),
    error = function(e) {
      warning("natural-scale Hessian failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(H)) dimnames(H) <- list(par_names, par_names)
  }
  structure(list(coefficients = theta, loglik = sum(individual_loglik),
    individual_loglik = individual_loglik,
    convergence = opt$convergence, message = opt$message,
    optimizer = opt, decode = decoded, nobs = N, npar = d,
    item_names = namvar, factor_names = factor_names, loadmask = loadmask,
    marker_item = namvar[marker], levels = levels_item,
    patterns = pat, probabilities = pi_hat, jacobian = J,
    pattern_scores = U, scores = U_obs, information = A,
    empirical_information = A_obs, hessian = H,
    observed_patterns = obs_patterns, observed_counts = counts,
    observed_pattern_id = id, observed_probabilities = exp(obs_logpi),
    probability_function = probability_function,
    jacobian_function = jacobian_function,
    integration = list(method = methInteg, dimensions = Q, nodes = nrow(nodes),
                       adaptive = methInteg == "AGH",
                       adaptation = if (individual_agh) "per individual" else NULL,
                       final_adaptations = if (individual_agh) final_adaptations else NULL)),
    class = "cross_sectional_probit_cfa")
}
