wu_stage_info <- function(stage) {
  stages <- c(
    baseline = 0L,
    threshold = 1L,
    threshold_loading = 2L,
    threshold_loading_intercept = 3L,
    strict = 4L
  )

  if (length(stage) != 1L || is.na(stage)) stop("stage non valido")

  if (is.character(stage)) {
    if (!stage %in% names(stages)) {
      stop("stage deve essere uno tra: ", paste(names(stages), collapse = ", "))
    }
    code <- unname(stages[[stage]])
    name <- stage
  } else {
    code <- as.integer(stage)
    if (!code %in% unname(stages)) stop("stage deve essere compreso tra 0 e 4")
    name <- names(stages)[match(code, stages)]
  }

  list(code = code, name = name)
}

wu_build_parameter_layout <- function(stage_code, n, nvar, Q, Kvec, namvar, loadmask, time_names, factor_names) {
  d <- n * Q
  ntre <- if (stage_code == 0L) n * sum(Kvec) else sum(Kvec)

  n_lambda <- switch(as.character(stage_code),
    `0` = n * nvar,
    `1` = n * nvar,
    `2` = nvar,
    `3` = nvar,
    `4` = nvar
  )

  n_gamma <- switch(as.character(stage_code),
    `0` = 0L,
    `1` = (n - 1L) * nvar,
    `2` = (n - 1L) * nvar,
    `3` = 0L,
    `4` = 0L
  )

  n_theta <- switch(as.character(stage_code),
    `0` = 0L,
    `1` = (n - 1L) * nvar,
    `2` = (n - 1L) * nvar,
    `3` = (n - 1L) * nvar,
    `4` = 0L
  )

  n_kappa <- if (stage_code >= 3L) (n - 1L) * Q else 0L
  n_corr <- d * (d - 1L) / 2L
  n_sd <- if (stage_code >= 2L) d - Q else 0L
  npm <- as.integer(ntre + n_lambda + n_gamma + n_theta + n_kappa + n_corr + n_sd)

  nm <- character(npm)
  pos <- 0L

  if (stage_code == 0L) {
    for (tt in seq_len(n)) for (k in seq_len(nvar)) {
      for (h in seq_len(Kvec[k])) {
        pos <- pos + 1L
        nm[pos] <- sprintf("tau_raw[%s,%s]_c%d", namvar[k], time_names[tt], h)
      }
    }
  } else {
    for (k in seq_len(nvar)) for (h in seq_len(Kvec[k])) {
      pos <- pos + 1L
      nm[pos] <- sprintf("tau_raw[%s]_c%d", namvar[k], h)
    }
  }

  if (stage_code == 0L) {
    for (tt in seq_len(n)) for (k in seq_len(nvar)) {
      q <- which(loadmask[k, ] == 1L)
      pos <- pos + 1L
      nm[pos] <- sprintf("logit_lambda[%s,%s,%s]", namvar[k], factor_names[q], time_names[tt])
    }
  } else if (stage_code == 1L) {
    for (tt in seq_len(n)) for (k in seq_len(nvar)) {
      q <- which(loadmask[k, ] == 1L)
      pos <- pos + 1L
      prefix <- if (tt == 1L) "logit_lambda" else "log_lambda"
      nm[pos] <- sprintf("%s[%s,%s,%s]", prefix, namvar[k], factor_names[q], time_names[tt])
    }
  } else {
    for (k in seq_len(nvar)) {
      q <- which(loadmask[k, ] == 1L)
      pos <- pos + 1L
      nm[pos] <- sprintf("logit_lambda[%s,%s]", namvar[k], factor_names[q])
    }
  }

  if (stage_code == 0L) {

  } else if (stage_code <= 2L) {
    for (tt in 2:n) for (k in seq_len(nvar)) {
      pos <- pos + 1L
      nm[pos] <- sprintf("gamma[%s,%s]", namvar[k], time_names[tt])
    }
  }

  if (stage_code == 0L) {

  } else if (stage_code <= 3L) {
    for (tt in 2:n) for (k in seq_len(nvar)) {
      pos <- pos + 1L
      nm[pos] <- sprintf("log_theta[%s,%s]", namvar[k], time_names[tt])
    }
  }

  if (stage_code >= 3L) {
    for (tt in 2:n) for (q in seq_len(Q)) {
      pos <- pos + 1L
      nm[pos] <- sprintf("kappa[%s,%s]", time_names[tt], factor_names[q])
    }
  }

  eta_names <- unlist(lapply(time_names, function(tt) paste0(tt, ":", factor_names)))
  for (row in 2:d) for (col in seq_len(row - 1L)) {
    pos <- pos + 1L
    nm[pos] <- sprintf("pcor_eta[%s,%s]", eta_names[row], eta_names[col])
  }

  if (stage_code >= 2L) {
    for (idx in (Q + 1L):d) {
      pos <- pos + 1L
      nm[pos] <- sprintf("log_sd_eta[%s]", eta_names[idx])
    }
  }

  if (pos != npm) stop("errore interno nel layout dei parametri")

  list(
    npm = npm,
    ntre = as.integer(ntre),
    names = nm,
    nblock = c(
      thresholds = ntre,
      loadings = n_lambda,
      intercepts = n_gamma,
      residuals = n_theta,
      factor_means = n_kappa,
      factor_correlations = n_corr,
      factor_sds = n_sd
    )
  )
}

wu_make_start <- function(stage_code, layout, n, nvar, Q, Kvec, monotone_tau, lambda_start = 0.7) {
  b <- numeric(layout$npm)
  pos <- 0L

  if (stage_code == 0L) {
    for (tt in seq_len(n)) for (k in seq_len(nvar)) {
      K <- Kvec[k]
      thr <- qnorm(seq_len(K) / (K + 1))
      raw <- if (monotone_tau) c(thr[1], log(diff(thr))) else thr
      b[pos + seq_len(K)] <- raw
      pos <- pos + K
    }
  } else {
    for (k in seq_len(nvar)) {
      K <- Kvec[k]
      thr <- qnorm(seq_len(K) / (K + 1))
      if (monotone_tau) {
        raw <- c(thr[1], log(diff(thr)))
      } else {
        raw <- thr
      }
      b[pos + seq_len(K)] <- raw
      pos <- pos + K
    }
  }

  if (stage_code == 0L) {
    idx <- pos + seq_len(n * nvar)
    b[idx] <- qlogis(lambda_start)
    pos <- pos + n * nvar
  } else if (stage_code == 1L) {
    for (tt in seq_len(n)) for (k in seq_len(nvar)) {
      pos <- pos + 1L
      b[pos] <- if (tt == 1L) qlogis(lambda_start) else log(lambda_start)
    }
  } else {
    idx <- pos + seq_len(nvar)
    b[idx] <- qlogis(lambda_start)
    pos <- pos + nvar
  }

  ng <- unname(layout$nblock["intercepts"])
  if (ng > 0L) pos <- pos + ng

  nt <- unname(layout$nblock["residuals"])
  if (nt > 0L) {
    idx <- pos + seq_len(nt)
    b[idx] <- log(1)
    pos <- pos + nt
  }

  nk <- unname(layout$nblock["factor_means"])
  if (nk > 0L) pos <- pos + nk

  nc <- unname(layout$nblock["factor_correlations"])
  if (nc > 0L) pos <- pos + nc

  ns <- unname(layout$nblock["factor_sds"])
  if (ns > 0L) pos <- pos + ns

  if (pos != layout$npm) stop("errore interno nei valori iniziali")
  names(b) <- layout$names
  b
}

wu_decode <- function(par, stage_code, n, nvar, Q, Kvec, namvar, loadmask, monotone_tau, time_names, factor_names) {
  pos <- 0L
  tau_time <- lapply(Kvec, function(K) matrix(NA_real_, K, n))
  names(tau_time) <- namvar

  if (stage_code == 0L) {
    for (tt in seq_len(n)) for (k in seq_len(nvar)) {
      K <- Kvec[k]
      u <- par[pos + seq_len(K)]
      pos <- pos + K
      if (monotone_tau) {
        thr <- numeric(K)
        thr[1] <- u[1]
        if (K > 1L) for (h in 2:K) thr[h] <- thr[h - 1L] + exp(u[h])
      } else {
        thr <- u
      }
      tau_time[[k]][, tt] <- thr
    }
  } else {
    for (k in seq_len(nvar)) {
      K <- Kvec[k]
      u <- par[pos + seq_len(K)]
      pos <- pos + K
      if (monotone_tau) {
        thr <- numeric(K)
        thr[1] <- u[1]
        if (K > 1L) for (h in 2:K) thr[h] <- thr[h - 1L] + exp(u[h])
      } else {
        thr <- u
      }
      tau_time[[k]] <- matrix(rep(thr, n), nrow = K, ncol = n)
    }
  }
  for (k in seq_len(nvar)) dimnames(tau_time[[k]]) <- list(paste0("c", seq_len(Kvec[k])), time_names)

  lambda_time <- array(0, c(nvar, Q, n), dimnames = list(namvar, factor_names, time_names))
  if (stage_code == 0L) {
    for (tt in seq_len(n)) for (k in seq_len(nvar)) {
      q <- which(loadmask[k, ] == 1L)
      pos <- pos + 1L
      lambda_time[k, q, tt] <- plogis(par[pos])
    }
  } else if (stage_code == 1L) {
    for (tt in seq_len(n)) for (k in seq_len(nvar)) {
      q <- which(loadmask[k, ] == 1L)
      pos <- pos + 1L
      lambda_time[k, q, tt] <- if (tt == 1L) plogis(par[pos]) else exp(par[pos])
    }
  } else {
    for (k in seq_len(nvar)) {
      q <- which(loadmask[k, ] == 1L)
      pos <- pos + 1L
      lam <- plogis(par[pos])
      lambda_time[k, q, ] <- lam
    }
  }

  gamma_time <- matrix(0, nvar, n, dimnames = list(namvar, time_names))
  if (stage_code == 0L) {

  } else if (stage_code <= 2L) {
    for (tt in 2:n) for (k in seq_len(nvar)) {
      pos <- pos + 1L
      gamma_time[k, tt] <- par[pos]
    }
  }

  theta_time <- matrix(NA_real_, nvar, n, dimnames = list(namvar, time_names))
  if (stage_code == 0L) {
    for (tt in seq_len(n)) for (k in seq_len(nvar)) {
      q <- which(loadmask[k, ] == 1L)
      theta_time[k, tt] <- 1 - lambda_time[k, q, tt]^2
    }
  } else {
    for (k in seq_len(nvar)) {
      q <- which(loadmask[k, ] == 1L)
      theta_time[k, 1] <- 1 - lambda_time[k, q, 1]^2
    }
    if (stage_code <= 3L) {
      for (tt in 2:n) for (k in seq_len(nvar)) {
        pos <- pos + 1L
        theta_time[k, tt] <- exp(par[pos])
      }
    } else {
      for (tt in 2:n) theta_time[, tt] <- theta_time[, 1]
    }
  }

  kappa <- matrix(0, n, Q, dimnames = list(time_names, factor_names))
  if (stage_code >= 3L) {
    for (tt in 2:n) for (q in seq_len(Q)) {
      pos <- pos + 1L
      kappa[tt, q] <- par[pos]
    }
  }

  d <- n * Q
  Lcorr <- matrix(0, d, d)
  Lcorr[1, 1] <- 1
  if (d > 1L) {
    for (row in 2:d) {
      prod <- 1
      for (col in seq_len(row - 1L)) {
        pos <- pos + 1L
        rho <- tanh(par[pos])
        Lcorr[row, col] <- rho * prod
        prod <- prod * sqrt(max(1e-14, 1 - rho^2))
      }
      Lcorr[row, row] <- prod
    }
  }

  sdev <- rep(1, d)
  if (stage_code >= 2L) {
    for (idx in (Q + 1L):d) {
      pos <- pos + 1L
      sdev[idx] <- exp(par[pos])
    }
  }

  Lphi <- diag(sdev, d, d) %*% Lcorr
  Phi <- Lphi %*% t(Lphi)
  eta_names <- unlist(lapply(time_names, function(tt) paste0(tt, ":", factor_names)))
  dimnames(Phi) <- dimnames(Lphi) <- list(eta_names, eta_names)

  if (pos != length(par)) stop("decode non allineato")

  sigma_diag_time <- matrix(NA_real_, nvar, n, dimnames = list(namvar, time_names))
  for (tt in seq_len(n)) for (k in seq_len(nvar)) {
    q <- which(loadmask[k, ] == 1L)
    idx_eta <- (tt - 1L) * Q + q
    sigma_diag_time[k, tt] <- lambda_time[k, q, tt]^2 * Phi[idx_eta, idx_eta] + theta_time[k, tt]
  }

  list(
    tau_time = tau_time,
    lambda_time = lambda_time,
    gamma_time = gamma_time,
    theta_time = theta_time,
    kappa = kappa,
    L_phi = Lphi,
    Phi = Phi,
    diag_Sigma_time = sigma_diag_time,
    diag_Sigma_t1 = sigma_diag_time[, 1]
  )
}

#fit one invariance stage
cfa_long_fit_wu2016 <- function(data,
                                namvar,
                                occ = "occ",
                                ID = "ID",
                                Q = NULL,
                                simple = NULL,
                                loadmask = NULL,
                                stage = c("baseline", "threshold", "threshold_loading", "threshold_loading_intercept", "strict"),
                                nquadgh = 7L,
                                update_gh = 0L,
                                maxiter = 500L,
                                monotone_tau = TRUE,
                                start = NULL,
                                parscale = NULL,
                                tol_loglik = 1e-7,
                                tol_par = 1e-3,
                                compute_hessian = FALSE,
                                verbose = FALSE) {

  stage <- match.arg(stage)
  st <- wu_stage_info(stage)
  stage_code <- st$code

  if (!native_has("loglik_long_tt_agh_full_equiv_wu")) {
    stop("simbolo loglik_long_tt_agh_full_equiv_wu non caricato")
  }

  maxiter <- check_count(maxiter, "maxiter", 0L)
  update_gh <- check_count(update_gh, "update_gh", 0L)
  packed <- prepare_long_data(data, namvar, occ, ID, min_categories = 3L, consecutive = TRUE)
  nvar <- length(namvar)
  nsub <- length(packed$subjects)
  n <- length(packed$times)
  if (n < 2L) stop("The invariance sequence requires at least two occasions.")
  Cvec <- packed$c
  Kvec <- Cvec - 1L
  if (stage_code == 1L && all(Kvec == 2L)) {
    warning("With only three-category items, threshold invariance is equivalent to the baseline model.")
  }

  loadmask0 <- loading_mask(namvar, Q, simple, loadmask)
  Q <- ncol(loadmask0)
  if (any(rowSums(loadmask0) != 1L)) {
    stop("The Wu and Estabrook models require one active loading per item.")
  }

  nmes0 <- packed$nmes
  Y0 <- as.integer(packed$y)
  nitem0 <- as.integer(packed$time)
  nobs <- length(Y0)

  time_names <- paste0("t", seq_len(n))
  factor_names <- paste0("F", seq_len(Q))
  layout <- wu_build_parameter_layout(stage_code, n, nvar, Q, Kvec, namvar, loadmask0, time_names, factor_names)

  gh <- make_gh_tensor(n * Q, nquadgh)
  nquad_re0 <- 1L
  nodes_re0 <- 0
  weights_re0 <- 1

  b0 <- if (is.null(start)) {
    wu_make_start(stage_code, layout, n, nvar, Q, Kvec, isTRUE(monotone_tau))
  } else {
    x <- as.double(start)
    if (length(x) != layout$npm) stop("start ha lunghezza diversa da npm")
    names(x) <- layout$names
    x
  }

  ps <- if (is.null(parscale)) rep(1, layout$npm) else as.double(parscale)
  if (length(ps) != layout$npm || any(!is.finite(ps)) || any(ps <= 0)) stop("parscale non valido")

  fit_id0 <- next_fit_id()

  call_ll <- function(par, update_grid = 1L) {
    native_call(
      "loglik_long_tt_agh_full_equiv_wu",
      b0 = as.double(par),
      Y0 = Y0,
      nvar0 = as.integer(nvar),
      nsub0 = as.integer(nsub),
      nobs0 = as.integer(nobs),
      n0 = as.integer(n),
      nmes0 = as.integer(nmes0),
      c0 = as.integer(Cvec),
      cm10 = as.integer(Kvec),
      ntre0 = as.integer(layout$ntre),
      nquad_eta0 = as.integer(gh$nquad),
      nitem0 = nitem0,
      nodes_eta0 = as.double(gh$nodes),
      weights_eta0 = as.double(gh$weights),
      nquad_re0 = as.integer(nquad_re0),
      nodes_re0 = as.double(nodes_re0),
      weights_re0 = as.double(weights_re0),
      npm0 = as.integer(layout$npm),
      re_item0 = 0L,
      Q0 = as.integer(Q),
      loadmask0 = as.integer(loadmask0),
      phi_struct0 = 0L,
      monotone_tau0 = as.integer(isTRUE(monotone_tau)),
      wu_stage0 = as.integer(stage_code),
      update_gh0 = as.integer(update_grid),
      fit_id0 = as.integer(fit_id0),
      vrais = double(1)
    )$vrais
  }

  fit_dims <- list(nsub = nsub, nvar = nvar, n = n, Q = Q, dim_eta = n * Q)
  fit_structure <- list(framework = "Wu & Estabrook 2016", methInteg = "AGH",
                        two_tier = TRUE, two_tier_type = "full_equiv",
                        re_item = FALSE, phi_struct = "FREE",
                        loadmask = loadmask0, nquad_eta = gh$nquad)
  condition <- c(
      baseline = "configural baseline: Condition 7 at each time (delta parametrization; nu_t=0, kappa_t=0, diag(Phi_tt)=I, diag(Sigma_y*_tt)=I)",
      threshold = "Wu-Estabrook Condition 15",
      threshold_loading = "Wu-Estabrook Condition 19, Proposition 7",
      threshold_loading_intercept = "Wu-Estabrook Condition 27, Proposition 11",
      strict = "full invariance, Condition 9 plus residual variance equality"
    )[[stage]]

  if (maxiter == 0L) {
    val <- call_ll(b0, 1L)
    if (!valid_objective(val)) stop("The initial invariance likelihood could not be evaluated.")
    dec <- wu_decode(b0, stage_code, n, nvar, Q, Kvec, namvar, loadmask0, isTRUE(monotone_tau), time_names, factor_names)
    return(new_long_fit(list(
      stage = stage,
      condition = condition,
      stage_code = stage_code,
      conv = NA_integer_,
      loglik = -val,
      negloglik = val,
      npar = layout$npm,
      par_opt = b0,
      decode = dec,
      dims = fit_dims,
      structure = fit_structure,
      layout = layout
    ), "Ordinal CFA with measurement invariance", packed))
  }

  if (update_gh <= 0L) {
    fn <- function(par) call_ll(par, 1L)
    opt <- optim(
      par = b0,
      fn = fn,
      method = "BFGS",
      control = list(maxit = as.integer(maxiter), parscale = ps, trace = if (verbose) 1 else 0)
    )
    best <- opt$par
    final_value <- call_ll(best, 1L)
  } else {
    best <- b0
    prev_best <- NULL
    prev_val <- NA_real_
    converged <- FALSE
    step <- as.integer(update_gh)
    outer_max <- max(1L, floor((sqrt(1 + 8 * maxiter / step) - 1) / 2))

    for (oo in seq_len(outer_max)) {
      val0 <- call_ll(best, 1L)
      if (!is.finite(val0) || val0 >= 1e8) stop("aggiornamento griglia AGH fallito")

      fn_fixed <- function(par) call_ll(par, 0L)
      oi <- optim(
        par = best,
        fn = fn_fixed,
        method = "BFGS",
        control = list(maxit = oo * step, parscale = ps, trace = if (verbose) 1 else 0)
      )
      best_new <- oi$par
      val_new <- call_ll(best_new, 1L)
      if (!is.finite(val_new) || val_new >= 1e8) stop("refresh AGH fallito")

      if (!is.null(prev_best) && is.finite(prev_val)) {
        rel_ll <- abs(val_new - prev_val) / max(1, abs(val_new), abs(prev_val))
        rel_par <- max(abs(best_new - prev_best) / pmax(1, abs(best_new), abs(prev_best)))
        if (rel_ll < tol_loglik && rel_par < tol_par) {
          converged <- TRUE
          best <- best_new
          prev_val <- val_new
          break
        }
      }

      prev_best <- best_new
      prev_val <- val_new
      best <- best_new
    }

    final_value <- call_ll(best, 1L)
    opt <- list(
      par = best,
      value = final_value,
      convergence = if (converged) 0L else 1L,
      message = if (converged) "converged" else "maximum iterations reached"
    )
  }

  names(best) <- layout$names
  dec <- wu_decode(best, stage_code, n, nvar, Q, Kvec, namvar, loadmask0, isTRUE(monotone_tau), time_names, factor_names)

  hess <- NULL
  vcov <- NULL
  se_opt <- rep(NA_real_, length(best))
  names(se_opt) <- names(best)
  if (isTRUE(compute_hessian)) {
    hess <- optimHess(best, function(par) call_ll(par, 1L))
    vcov <- tryCatch(solve(hess), error = function(e) NULL)
    if (!is.null(vcov)) {
      vv <- diag(vcov)
      vv[vv < 0] <- NA_real_
      se_opt <- sqrt(vv)
      names(se_opt) <- names(best)
    }
  }

  if (!valid_objective(final_value)) stop("The final invariance likelihood could not be evaluated.")
  new_long_fit(list(
    stage = stage,
    stage_code = stage_code,
    condition = condition,
    conv = opt$convergence,
    message = opt$message,
    loglik = -final_value,
    negloglik = final_value,
    npar = layout$npm,
    par_opt = best,
    hessian = hess,
    vcov_opt = vcov,
    se_opt = se_opt,
    decode = dec,
    layout = layout,
    dims = fit_dims,
    structure = fit_structure
  ), "Ordinal CFA with measurement invariance", packed)
}

#compare successive invariance stages
measurement_invariance_wu2016 <- function(data,
                                          namvar,
                                          occ = "occ",
                                          ID = "ID",
                                          Q = NULL,
                                          simple = NULL,
                                          loadmask = NULL,
                                          stages = c("baseline", "threshold", "threshold_loading", "threshold_loading_intercept", "strict"),
                                          ...) {
  valid <- c("baseline", "threshold", "threshold_loading", "threshold_loading_intercept", "strict")
  if (any(!stages %in% valid)) stop("stages non validi")

  dots <- list(...)
  if (!is.null(dots$start)) stop("nella sequenza completa non passare start, perch\u00e9 la dimensione cambia tra modelli")

  fits <- vector("list", length(stages))
  names(fits) <- stages

  for (j in seq_along(stages)) {
    if (isTRUE(dots$verbose)) cat("\nfit:", stages[j], "\n")
    fits[[j]] <- do.call(
      cfa_long_fit_wu2016,
      c(
        list(
          data = data,
          namvar = namvar,
          occ = occ,
          ID = ID,
          Q = Q,
          simple = simple,
          loadmask = loadmask,
          stage = stages[j]
        ),
        dots
      )
    )
  }

  tab <- data.frame(
    stage = stages,
    condition = vapply(fits, `[[`, character(1), "condition"),
    logLik = vapply(fits, `[[`, numeric(1), "loglik"),
    npar = vapply(fits, `[[`, numeric(1), "npar"),
    LR_vs_previous = NA_real_,
    df_vs_previous = NA_integer_,
    p_vs_previous = NA_real_,
    stringsAsFactors = FALSE
  )

  if (length(stages) > 1L) {
    for (j in 2:length(stages)) {
      LR <- 2 * (tab$logLik[j - 1L] - tab$logLik[j])
      df <- as.integer(tab$npar[j - 1L] - tab$npar[j])
      tab$LR_vs_previous[j] <- LR
      tab$df_vs_previous[j] <- df
      both_converged <- isTRUE(fits[[j - 1L]]$conv == 0L) && isTRUE(fits[[j]]$conv == 0L)
      if (df > 0L && LR >= -1e-8 && both_converged) {
        tab$p_vs_previous[j] <- pchisq(max(LR, 0), df = df, lower.tail = FALSE)
      }
    }
  }

  structure(list(table = tab, fits = fits), class = "wu2016_mi")
}

print.wu2016_mi <- function(x, ...) {
  print(x$table, row.names = FALSE)
  invisible(x)
}
