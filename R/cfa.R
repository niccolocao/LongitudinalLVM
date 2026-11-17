set_cfa_long_options <- function(marker0, Q0, nvar0, nsub0, n0, p_fix0 = 0L, X0 = NULL) {
  if (is.null(X0)) X0 <- 0
  out <- native_call("set_cfa_long_options", marker_item0 = as.integer(marker0), Q0 = as.integer(Q0),
            nvar0 = as.integer(nvar0), nsub0 = as.integer(nsub0), n0 = as.integer(n0),
            p_fix0 = as.integer(p_fix0), X0 = as.double(X0), ok0 = integer(1))
  if (out$ok0 != 1L) stop("Impossibile inizializzare marker/covariate CFA nel codice Fortran.")
  invisible(TRUE)
}

cfa_threshold_id_code <- function(threshold_id) {
  if (is.numeric(threshold_id) || is.integer(threshold_id)) {
    if (length(threshold_id) != 1L || is.na(threshold_id) ||
        !threshold_id %in% c(0L, 1L)) {
      stop("threshold_id numerico deve essere 0 oppure 1.")
    }
    return(as.integer(threshold_id))
  }

  threshold_id <- match.arg(
    as.character(threshold_id),
    c("tau2_tau3", "tau1_tau2")
  )

  if (identical(threshold_id, "tau2_tau3")) 0L else 1L
}

set_cfa_long_threshold_id <- function(threshold_id = "tau2_tau3") {
  threshold_id0 <- cfa_threshold_id_code(threshold_id)

  if (!native_has("set_cfa_long_threshold_id")) {
    stop("Simbolo set_cfa_long_threshold_id non caricato. Ricompila e ricarica la DLL.")
  }

  out <- native_call(
    "set_cfa_long_threshold_id",
    threshold_id0 = as.integer(threshold_id0),
    ok0 = integer(1)
  )

  if (out$ok0 != 1L) {
    stop("Impossibile impostare l'identificazione delle soglie CFA nel codice Fortran.")
  }

  invisible(threshold_id0)
}

cfa_mean_ident_code <- function(mean_ident) {
  if (is.numeric(mean_ident) || is.integer(mean_ident)) {
    if (length(mean_ident) != 1L || is.na(mean_ident) ||
        !mean_ident %in% c(0L, 1L)) {
      stop("mean_ident numerico deve essere 0 oppure 1.")
    }
    return(as.integer(mean_ident))
  }

  mean_ident <- match.arg(
    as.character(mean_ident),
    c("latent_baseline", "marker_intercept")
  )

  if (identical(mean_ident, "marker_intercept")) 1L else 0L
}

set_cfa_long_gamma_time_free <- function(gamma_time_free0) {
  mask <- as.matrix(gamma_time_free0)
  if (!native_has("set_cfa_long_gamma_time_free")) {
    stop("Simbolo set_cfa_long_gamma_time_free non caricato. Ricompila e ricarica la DLL.")
  }
  out <- native_call(
    "set_cfa_long_gamma_time_free",
    gamma_time_free0 = as.integer(mask),
    nvar0 = as.integer(nrow(mask)),
    n0 = as.integer(ncol(mask)),
    ok0 = integer(1)
  )
  if (out$ok0 != 1L) stop("Impossibile impostare le intercette time-specific della CFA.")
  invisible(mask)
}

set_cfa_long_mean_ident <- function(mean_ident = "latent_baseline") {
  mean_ident0 <- cfa_mean_ident_code(mean_ident)

  if (!native_has("set_cfa_long_mean_ident")) {
    stop("Simbolo set_cfa_long_mean_ident non caricato. Ricompila e ricarica la DLL.")
  }

  out <- native_call(
    "set_cfa_long_mean_ident",
    mean_ident0 = as.integer(mean_ident0),
    ok0 = integer(1)
  )

  if (out$ok0 != 1L) {
    stop("Impossibile impostare l'identificazione della media CFA nel codice Fortran.")
  }

  invisible(mean_ident0)
}

loglik_long <- function(b0,
                        Y0, nvar0, nsub0, nobs0, n0,
                        nmes0, c0, cm10, ntre0,
                        nquad0, nitem0,
                        nodes0, weights0,
                        npm0,
                        re_item0,
                        theta_time0 = 0L,
                        Q0 = 1L,
                        loadmask0,
                        phi_struct0,
                        monotone_tau0,
                        fit_id0) {

  native_call("loglik_long",
     b0          = as.double(b0),
     Y0          = as.integer(Y0),
     nvar0       = as.integer(nvar0),
     nsub0       = as.integer(nsub0),
     nobs0       = as.integer(nobs0),
     n0          = as.integer(n0),
     nmes0       = as.integer(nmes0),
     c0          = as.integer(c0),
     cm10        = as.integer(cm10),
     ntre0       = as.integer(ntre0),
     nquad0      = as.integer(nquad0),
     nitem0      = as.integer(nitem0),
     nodes0      = as.double(nodes0),
     weights0    = as.double(weights0),
     npm0        = as.integer(npm0),
     re_item0    = as.integer(re_item0),
     theta_time0 = as.integer(theta_time0),
     Q0          = as.integer(Q0),
     loadmask0   = as.integer(loadmask0),
     phi_struct0 = as.integer(phi_struct0),
     monotone_tau0 = as.integer(monotone_tau0),
     fit_id0     = as.integer(fit_id0),
     vrais       = double(1)
  )$vrais
}

loglik_long_agh <- function(b0,
                        Y0, nvar0, nsub0, nobs0, n0,
                        nmes0, c0, cm10, ntre0,
                        nquad0, nitem0,
                        nodes0, weights0,
                        npm0,
                        re_item0,
                        theta_time0 = 0L,
                        Q0 = 1L,
                        loadmask0,
                        phi_struct0,
                        monotone_tau0,
                        update_gh0,
                        fit_id0) {

  native_call("loglik_long_agh",
     b0         = as.double(b0),
     Y0         = as.integer(Y0),
     nvar0      = as.integer(nvar0),
     nsub0      = as.integer(nsub0),
     nobs0      = as.integer(nobs0),
     n0         = as.integer(n0),
     nmes0      = as.integer(nmes0),
     c0         = as.integer(c0),
     cm10       = as.integer(cm10),
     ntre0      = as.integer(ntre0),
     nquad0     = as.integer(nquad0),
     nitem0     = as.integer(nitem0),
     nodes0     = as.double(nodes0),
     weights0   = as.double(weights0),
     npm0       = as.integer(npm0),
     re_item0   = as.integer(re_item0),
     theta_time0= as.integer(theta_time0),
     Q0         = as.integer(Q0),
     loadmask0  = as.integer(loadmask0),
     phi_struct0 = as.integer(phi_struct0),
     monotone_tau0 = as.integer(monotone_tau0),
     update_gh0 = as.integer(update_gh0),
     fit_id0    = as.integer(fit_id0),
     vrais      = double(1)
  )$vrais
}

loglik_long_tt <- function(b0,
                           Y0, nvar0, nsub0, nobs0, n0,
                           nmes0, c0, cm10, ntre0,
                           nquad_eta0,
                           nodes_eta0, weights_eta0,
                           nquad_re0,
                           nodes_re0, weights_re0,
                           nitem0,
                           npm0,
                           re_item0,
                           theta_time0 = 0L,
                           Q0 = 1L,
                           loadmask0,
                           phi_struct0,
                           monotone_tau0,
                           fit_id0) {

  if (re_item0 == 0L) {
    nquad_re0   <- 1L
    nodes_re0   <- 0
    weights_re0 <- 1
  }

  native_call("loglik_long_tt",
     b0            = as.double(b0),
     Y0            = as.integer(Y0),
     nvar0         = as.integer(nvar0),
     nsub0         = as.integer(nsub0),
     nobs0         = as.integer(nobs0),
     n0            = as.integer(n0),
     nmes0         = as.integer(nmes0),
     c0            = as.integer(c0),
     cm10          = as.integer(cm10),
     ntre0         = as.integer(ntre0),
     nquad_eta0    = as.integer(nquad_eta0),
     nodes_eta0    = as.double(nodes_eta0),
     weights_eta0  = as.double(weights_eta0),
     nquad_re0     = as.integer(nquad_re0),
     nodes_re0     = as.double(nodes_re0),
     weights_re0   = as.double(weights_re0),
     nitem0        = as.integer(nitem0),
     npm0          = as.integer(npm0),
     re_item0      = as.integer(re_item0),
     theta_time0   = as.integer(theta_time0),
     Q0            = as.integer(Q0),
     loadmask0     = as.integer(loadmask0),
     phi_struct0   = as.integer(phi_struct0),
     monotone_tau0 = as.integer(monotone_tau0),
     fit_id0       = as.integer(fit_id0),
     vrais         = double(1)
  )$vrais
}

loglik_long_tt_agh <- function(b0,
                               Y0, nvar0, nsub0, nobs0, n0,
                               nmes0, c0, cm10, ntre0,
                               nquad_eta0,
                               nodes_eta0, weights_eta0,
                               nquad_re0,
                               nodes_re0, weights_re0,
                               nitem0,
                               npm0,
                               re_item0,
                               theta_time0 = 0L,
                               Q0 = 1L,
                               loadmask0,
                               phi_struct0,
                               monotone_tau0,
                               update_gh0,
                               fit_id0) {

  if (re_item0 == 0L) {
    nquad_re0   <- 1L
    nodes_re0   <- 0
    weights_re0 <- 1
  }

  native_call("loglik_long_tt_agh",
     b0            = as.double(b0),
     Y0            = as.integer(Y0),
     nvar0         = as.integer(nvar0),
     nsub0         = as.integer(nsub0),
     nobs0         = as.integer(nobs0),
     n0            = as.integer(n0),
     nmes0         = as.integer(nmes0),
     c0            = as.integer(c0),
     cm10          = as.integer(cm10),
     ntre0         = as.integer(ntre0),
     nquad_eta0    = as.integer(nquad_eta0),
     nitem0        = as.integer(nitem0),
     nodes_eta0    = as.double(nodes_eta0),
     weights_eta0  = as.double(weights_eta0),
     nquad_re0     = as.integer(nquad_re0),
     nodes_re0     = as.double(nodes_re0),
     weights_re0   = as.double(weights_re0),
     npm0          = as.integer(npm0),
     re_item0      = as.integer(re_item0),
     theta_time0   = as.integer(theta_time0),
     Q0            = as.integer(Q0),
     loadmask0     = as.integer(loadmask0),
     phi_struct0   = as.integer(phi_struct0),
     monotone_tau0 = as.integer(monotone_tau0),
     update_gh0    = as.integer(update_gh0),
     fit_id0       = as.integer(fit_id0),
     vrais         = double(1)
  )$vrais
}

loglik_long_tt_agh_full_equiv <- function(b0,
                               Y0, nvar0, nsub0, nobs0, n0,
                               nmes0, c0, cm10, ntre0,
                               nquad_eta0,
                               nodes_eta0, weights_eta0,
                               nquad_re0,
                               nodes_re0, weights_re0,
                               nitem0,
                               npm0,
                               re_item0,
                               theta_time0 = 0L,
                               meas_invar0 = c(1L, 1L, 1L, 1L),
                               Q0 = 1L,
                               loadmask0,
                               phi_struct0,
                               monotone_tau0,
                               update_gh0,
                               fit_id0) {

  if (re_item0 == 0L) {
    nquad_re0   <- 1L
    nodes_re0   <- 0
    weights_re0 <- 1
  }

  native_call("loglik_long_tt_agh_full_equiv",
     b0            = as.double(b0),
     Y0            = as.integer(Y0),
     nvar0         = as.integer(nvar0),
     nsub0         = as.integer(nsub0),
     nobs0         = as.integer(nobs0),
     n0            = as.integer(n0),
     nmes0         = as.integer(nmes0),
     c0            = as.integer(c0),
     cm10          = as.integer(cm10),
     ntre0         = as.integer(ntre0),
     nquad_eta0    = as.integer(nquad_eta0),
     nitem0        = as.integer(nitem0),
     nodes_eta0    = as.double(nodes_eta0),
     weights_eta0  = as.double(weights_eta0),
     nquad_re0     = as.integer(nquad_re0),
     nodes_re0     = as.double(nodes_re0),
     weights_re0   = as.double(weights_re0),
     npm0          = as.integer(npm0),
     re_item0      = as.integer(re_item0),
     theta_time0   = as.integer(theta_time0),
     meas_invar0   = as.integer(meas_invar0),
     Q0            = as.integer(Q0),
     loadmask0     = as.integer(loadmask0),
     phi_struct0   = as.integer(phi_struct0),
     monotone_tau0 = as.integer(monotone_tau0),
     update_gh0    = as.integer(update_gh0),
     fit_id0       = as.integer(fit_id0),
     vrais         = double(1)
  )$vrais
}

#fit an ordinal longitudinal cfa
cfa_long_fit_md <- function(data,
                            namvar,
                            occ = "occ",
                            ID  = "ID",
                            Q = NULL,
                            simple = NULL,
                            loadmask = NULL,
                            marker = NULL,
                            covariates = NULL,
                            post_covariates = NULL,
                            methInteg = c("QMC","GH","AGH"),
                            update_gh = 0,
                            nquadqmc = 2000,
                            nquadgh  = 8,
                            re_item  = TRUE,
                            two_tier = FALSE,
                            two_tier_type = c("node_specific", "full_equiv"),
                            theta_time = FALSE,
                            meas_invar = NULL,
                            intercepts_free_items = NULL,
                            monotone_tau = TRUE,

                            phi_struct = "FREE",

                            maxiter = 500,
                            start = NULL,
                            verbose = FALSE,
                            compute_hessian = FALSE,
                            analytic_score = TRUE,
                            hess_eps = 1e-4,
                            profile = TRUE,
                            profile_print = TRUE,
                            nquadqmc_re = NULL,
                            re_tier = c("GH", "QMC"),
                            parscale = NULL,
                            tol_loglik = 1e-7,
                            tol_par = 1e-3,
                            threshold_id = c("tau1_tau2", "tau2_tau3"),
                            mean_ident = c("latent_baseline", "marker_intercept")) {

  maxiter <- check_count(maxiter, "maxiter", 0L)
  update_gh <- check_count(update_gh, "update_gh", 0L)
  re_item <- check_flag(re_item, "re_item")
  two_tier <- check_flag(two_tier, "two_tier")
  methInteg  <- match.arg(methInteg)
  if (methInteg == "QMC" && two_tier) {
    stop("QMC requires two_tier = FALSE.", call. = FALSE)
  }
  phi_struct <- match.arg(phi_struct, "FREE")
  two_tier_type <- match.arg(two_tier_type)
  re_tier <- match.arg(re_tier)
  threshold_id <- match.arg(threshold_id)
  threshold_id0 <- cfa_threshold_id_code(threshold_id)
  mean_ident <- match.arg(mean_ident)
  mean_ident0 <- cfa_mean_ident_code(mean_ident)

  if (length(tol_loglik) != 1L ||
      !is.finite(tol_loglik) ||
      tol_loglik <= 0) {
    stop("tol_loglik deve essere uno scalare finito > 0.")
  }

  if (length(tol_par) != 1L ||
      !is.finite(tol_par) ||
      tol_par <= 0) {
    stop("tol_par deve essere uno scalare finito > 0.")
  }

  if (length(hess_eps) != 1L ||
      !is.finite(hess_eps) ||
      hess_eps <= 0) {
    stop("hess_eps deve essere uno scalare finito > 0.")
  }

  eval_ll <- function(par, use_agh, update_gh_eval) {
    t0_fn <- proc.time()[["elapsed"]]

    val <- loglik_long_call(
      par = par,
      use_agh = use_agh,
      update_gh_eval = update_gh_eval
    )

    if (isTRUE(profile)) {
      dt <- proc.time()[["elapsed"]] - t0_fn
      nfn <<- nfn + 1L
      fn_elapsed <<- c(fn_elapsed, dt)
    }

    val
  }

make_two_tier_nodes <- function(dim_eta, re_item0, nquadgh) {
  if (dim_eta >= 9L) {
    warning(
      "griglia tensoriale su eta con dim_eta = ", dim_eta,
      ": nquad_eta = nquadgh^dim_eta pu\u00f2 essere molto grande."
    )
  }

  eta_rule <- make_gh_tensor(dim_eta, nquadgh)

  if (re_item0 == 1L) {
    re_rule <- make_gh_1d(nquadgh)
  } else {
    re_rule <- list(
      nquad = 1L,
      nodes = 0,
      weights = 1
    )
  }

  list(
    nquad_eta0 = eta_rule$nquad,
    nodes_eta0 = eta_rule$nodes,
    weights_eta0 = eta_rule$weights,
    nquad_re0 = re_rule$nquad,
    nodes_re0 = re_rule$nodes,
    weights_re0 = re_rule$weights
  )
}

tau_se_vector <- function(u_all, V_all, cm1, cm3, namvar, monotone_tau0,
                          n = 1L, time_names = NULL, tau_invar0 = 1L,
                          threshold_id0 = 0L) {
  out_names <- character(0)
  out_se <- numeric(0)

  if (is.null(time_names)) {
    time_names <- paste0("t", seq_len(n))
  }

  if (length(threshold_id0) != 1L ||
      is.na(threshold_id0) ||
      !threshold_id0 %in% c(0L, 1L)) {
    stop("threshold_id0 deve essere 0 oppure 1.")
  }

  add_item <- function(k, u, V, tt = NULL) {
    pref <- if (is.null(tt)) {
      sprintf("tau[%s]", namvar[k])
    } else {
      sprintf("tau[%s,%s]", namvar[k], time_names[tt])
    }

    nf <- cm3[k]

    if (threshold_id0 == 0L) {
      if (is.null(V)) {
        se1 <- NA_real_
      } else if (monotone_tau0 == 0L) {
        se1 <- if (is.finite(V[1, 1]) && V[1, 1] >= 0) {
          sqrt(V[1, 1])
        } else {
          NA_real_
        }
      } else {
        g <- numeric(nf)
        g[1] <- -exp(u[1])
        vv <- as.numeric(t(g) %*% V %*% g)
        se1 <- if (is.finite(vv) && vv >= 0) sqrt(vv) else NA_real_
      }

      out_names <<- c(out_names, paste0(pref, "_c1"))
      out_se <<- c(out_se, se1)

      out_names <<- c(out_names, paste0(pref, "_c2"))
      out_se <<- c(out_se, 0)

      out_names <<- c(out_names, paste0(pref, "_c3"))
      out_se <<- c(out_se, 0)

      if (cm1[k] > 3L) {
        for (cc in 4:cm1[k]) {
          h <- cc - 2L

          if (is.null(V)) {
            se <- NA_real_
          } else if (monotone_tau0 == 0L) {
            vhh <- V[h, h]
            se <- if (is.finite(vhh) && vhh >= 0) sqrt(vhh) else NA_real_
          } else {
            g <- numeric(nf)
            g[2:h] <- exp(u[2:h])
            vv <- as.numeric(t(g) %*% V %*% g)
            se <- if (is.finite(vv) && vv >= 0) sqrt(vv) else NA_real_
          }

          out_names <<- c(out_names, paste0(pref, "_c", cc))
          out_se <<- c(out_se, se)
        }
      }
    } else {
      out_names <<- c(
        out_names,
        paste0(pref, "_c1"),
        paste0(pref, "_c2")
      )
      out_se <<- c(out_se, 0, 0)

      if (cm1[k] > 2L) {
        for (cc in 3:cm1[k]) {
          h <- cc - 2L

          if (is.null(V)) {
            se <- NA_real_
          } else if (monotone_tau0 == 0L) {
            vhh <- V[h, h]
            se <- if (is.finite(vhh) && vhh >= 0) sqrt(vhh) else NA_real_
          } else {
            g <- numeric(nf)
            g[seq_len(h)] <- exp(u[seq_len(h)])
            vv <- as.numeric(t(g) %*% V %*% g)
            se <- if (is.finite(vv) && vv >= 0) sqrt(vv) else NA_real_
          }

          out_names <<- c(out_names, paste0(pref, "_c", cc))
          out_se <<- c(out_se, se)
        }
      }
    }
  }

  jj <- 0L

  if (tau_invar0 == 1L) {
    for (k in seq_along(namvar)) {
      nf <- cm3[k]
      idx <- (jj + 1L):(jj + nf)
      u <- u_all[idx]
      V <- if (!is.null(V_all)) V_all[idx, idx, drop = FALSE] else NULL
      add_item(k, u, V)
      jj <- jj + nf
    }
  } else {
    for (tt in seq_len(n)) {
      for (k in seq_along(namvar)) {
        nf <- cm3[k]
        idx <- (jj + 1L):(jj + nf)
        u <- u_all[idx]
        V <- if (!is.null(V_all)) V_all[idx, idx, drop = FALSE] else NULL
        add_item(k, u, V, tt)
        jj <- jj + nf
      }
    }
  }

  names(out_se) <- out_names
  out_se
}

  nvar <- length(namvar)
  covariates <- if (is.null(covariates)) character(0) else as.character(covariates)
  post_covariates <- if (is.null(post_covariates)) character(0) else as.character(post_covariates)
  if (anyDuplicated(covariates)) stop("covariates contiene nomi duplicati")
  if (anyDuplicated(post_covariates)) stop("post_covariates contiene nomi duplicati")
  if (length(covariates) && any(!covariates %in% names(data))) {
    stop("Covariate non trovate in data: ", paste(setdiff(covariates, names(data)), collapse = ", "))
  }
  if (length(post_covariates) && any(!post_covariates %in% covariates)) {
    stop("post_covariates deve essere un sottoinsieme di covariates: ",
         paste(setdiff(post_covariates, covariates), collapse = ", "))
  }
  packed <- prepare_long_data(data, namvar, occ, ID, covariates, min_categories = 4L)
  dat0 <- packed$data
  IDs <- packed$subjects
  n <- length(packed$times)
  nsub <- length(IDs)
  c0 <- as.integer(packed$c)
  cm10 <- c0 - 1L
  nmes0 <- packed$nmes
  Y0 <- as.integer(packed$y)
  nitem0 <- as.integer(packed$time)
  nobs <- length(Y0)

  loadmask0 <- loading_mask(namvar, Q, simple, loadmask)
  Q <- ncol(loadmask0)

  resolve_marker <- function(marker, loadmask0, namvar, Q) {
    if (is.null(marker)) {
      idx <- vapply(seq_len(Q), function(q) which(loadmask0[, q] == 1L)[1], integer(1))
    } else if (is.character(marker)) {
      if (length(marker) != Q) stop("marker deve avere lunghezza Q")
      idx <- match(marker, namvar)
      if (anyNA(idx)) stop("Marker non trovato in namvar")
    } else {
      idx <- as.integer(marker)
      if (length(idx) != Q || anyNA(idx) || any(idx < 1L | idx > nvar)) stop("marker numerico non valido")
    }
    for (q in seq_len(Q)) if (loadmask0[idx[q], q] != 1L) stop("Il marker del fattore ", q, " non ha loading attivo.")
    idx
  }
  marker_idx <- resolve_marker(marker, loadmask0, namvar, Q)
  names(marker_idx) <- paste0("F", seq_len(Q))

  marker_mean_ident <- identical(mean_ident0, 1L)
  if (marker_mean_ident &&
      !(isTRUE(two_tier) &&
        identical(methInteg, "AGH") &&
        identical(two_tier_type, "full_equiv"))) {
    stop("mean_ident='marker_intercept' \u00e8 implementato solo con two_tier=TRUE, methInteg='AGH', two_tier_type='full_equiv'")
  }

  gamma_free_item <- rep(TRUE, nvar)
  if (marker_mean_ident) {
    gamma_free_item[unique(marker_idx)] <- FALSE
  }

  kappa_times <- if (marker_mean_ident) {
    seq_len(n)
  } else if (n > 1L) {
    2:n
  } else {
    integer(0L)
  }

  p_fix <- length(covariates)
  if (p_fix > 0L) {
    for (v in covariates) if (!is.numeric(dat0[[v]])) stop("La covariata ", v, " deve essere numerica.")

    cov_time_mask <- matrix(
      1,
      nrow = n,
      ncol = p_fix,
      dimnames = list(paste0("t", seq_len(n)), covariates)
    )
    if (length(post_covariates)) cov_time_mask[1L, post_covariates] <- 0

    X0 <- array(0, dim = c(nsub, n, Q, p_fix),
                dimnames = list(as.character(IDs), paste0("t", seq_len(n)), paste0("F", seq_len(Q)), covariates))
    for (ii in seq_len(nsub)) {
      di <- dat0[dat0[[ID]] == IDs[ii], c("..t", covariates), drop = FALSE]
      if (anyDuplicated(di$..t)) stop("Pi\u00f9 di una riga per ID-tempo per ID=", IDs[ii])
      for (rr in seq_len(nrow(di))) {
        tt <- di$..t[rr]
        xv <- as.numeric(di[rr, covariates, drop = TRUE])
        active <- cov_time_mask[tt, ] == 1
        if (any(active & !is.finite(xv))) {
          bad <- covariates[active & !is.finite(xv)]
          stop("Covariate mancanti/non finite per ID=", IDs[ii], ", tempo=", tt, ": ", paste(bad, collapse = ", "))
        }
        xv[!active] <- 0
        for (q in seq_len(Q)) X0[ii, tt, q, ] <- xv
      }
    }
  } else {
    cov_time_mask <- matrix(integer(0), nrow = n, ncol = 0L)
    X0 <- 0
  }

  if (is.null(meas_invar)) {
    meas_invar0 <- c(
      thresholds = 1L,
      intercepts = 1L,
      loadings = 1L,
      residuals = if (isTRUE(theta_time)) 0L else 1L
    )
  } else {
    meas_invar0 <- as.integer(meas_invar)
    if (length(meas_invar0) != 4L) {
      stop("meas_invar deve avere lunghezza 4: thresholds, intercepts, loadings, residuals")
    }
    if (any(!meas_invar0 %in% c(0L, 1L))) {
      stop("meas_invar deve contenere solo 0/1 oppure TRUE/FALSE")
    }
    names(meas_invar0) <- c("thresholds", "intercepts", "loadings", "residuals")
  }

  if (any(meas_invar0[c("thresholds", "intercepts", "loadings")] == 0L) &&
      !(isTRUE(two_tier) && methInteg == "AGH" && identical(two_tier_type, "full_equiv"))) {
    stop("thresholds/intercepts/loadings time-specific sono implementati solo con two_tier=TRUE, methInteg='AGH', two_tier_type='full_equiv'")
  }

  if (is.null(intercepts_free_items)) intercepts_free_items <- character(0)
  if (!is.character(intercepts_free_items) || anyNA(intercepts_free_items) ||
      anyDuplicated(intercepts_free_items)) {
    stop("intercepts_free_items deve essere un vettore di nomi di item distinti.")
  }
  unknown_free <- setdiff(intercepts_free_items, namvar)
  if (length(unknown_free)) {
    stop("Item sconosciuti in intercepts_free_items: ", paste(unknown_free, collapse = ", "))
  }
  if (length(intercepts_free_items)) {
    if (!(isTRUE(two_tier) && methInteg == "AGH" && identical(two_tier_type, "full_equiv"))) {
      stop("intercepts_free_items \u00e8 supportato solo con AGH two_tier full_equiv.")
    }
    if (meas_invar0[["intercepts"]] != 1L) {
      stop("Con intercepts_free_items usa meas_invar['intercepts'] = 1L.")
    }
    if (any(intercepts_free_items %in% namvar[marker_idx]) && marker_mean_ident) {
      stop("Non \u00e8 possibile liberare le intercette del marker con mean_ident='marker_intercept'.")
    }
  }
  gamma_time_free0 <- matrix(0L, nrow = nvar, ncol = n)
  if (length(intercepts_free_items) && n > 1L) {
    gamma_time_free0[match(intercepts_free_items, namvar), 2:n] <- 1L
  }

  theta_time0 <- 1L - meas_invar0[["residuals"]]
  re_item0    <- if (isTRUE(re_item)) 1L else 0L
  monotone_tau0 <- if (isTRUE(monotone_tau)) 1L else 0L
  update_gh0 <- update_gh

  dim_eta_int <- n * Q
  dimInt <- n * Q + (if (re_item) nvar else 0L)
nquad0 <- NA_integer_
nodes0 <- NULL
weights0 <- NULL

nquad_eta0 <- NA_integer_
nodes_eta0 <- NULL
weights_eta0 <- NULL

nquad_re0 <- NA_integer_
nodes_re0 <- NULL
weights_re0 <- NULL

if (!isTRUE(two_tier)) {

  if (methInteg %in% c("GH", "AGH")) {
    full_rule <- make_gh_tensor(dimInt, nquadgh)

    nquad0   <- full_rule$nquad
    nodes0   <- full_rule$nodes
    weights0 <- full_rule$weights

  } else {
    full_rule <- make_qmc(dimInt, nquadqmc)

    nquad0   <- full_rule$nquad
    nodes0   <- full_rule$nodes
    weights0 <- full_rule$weights
  }

} else {

  quad_obj <- make_two_tier_nodes(dim_eta_int, re_item0, nquadgh)

  nquad_eta0   <- quad_obj$nquad_eta0
  nodes_eta0   <- quad_obj$nodes_eta0
  weights_eta0 <- quad_obj$weights_eta0

  nquad_re0    <- quad_obj$nquad_re0
  nodes_re0    <- quad_obj$nodes_re0
  weights_re0  <- quad_obj$weights_re0
}

if (isTRUE(two_tier) && !isTRUE(re_item)) {
  message("two_tier=TRUE ma re_item=FALSE: integrazione solo sui fattori eta.")
}

  cm1 <- cm10
  cm3 <- pmax(cm1 - 2L, 0L)

  ntre_base <- sum(cm3)
  ntre0 <- as.integer(if (meas_invar0[["thresholds"]] == 1L) ntre_base else n * ntre_base)

  n_lambda_base <- sum(loadmask0) - Q
  n_lambda_free <- if (meas_invar0[["loadings"]] == 1L) n_lambda_base else n * n_lambda_base

  n_gamma_base <- sum(gamma_free_item)
  n_gamma_est <- if (meas_invar0[["intercepts"]] == 1L) {
    n_gamma_base + sum(gamma_time_free0)
  } else {
    n * n_gamma_base
  }

  n_theta_base <- sum(cm1 > 1L)
  n_theta_est <- if (meas_invar0[["residuals"]] == 1L) n_theta_base else n * n_theta_base

  n_kappa_est <- length(kappa_times) * Q
  dim_eta <- n*Q

  n_phi_par <- dim_eta * (dim_eta + 1L) / 2L

  npm0 <- as.integer(ntre0 +
                     n_lambda_free +
                     n_gamma_est +
                     n_theta_est +
                     (if (re_item) nvar else 0L) +
                     n_kappa_est +
                     Q*p_fix +
                     n_phi_par)

  make_start <- function() {
    b <- numeric(npm0)
    pos <- 0L

    tau_start_values <- function(nf) {
      if (nf <= 0L) return(numeric(0))

      if (threshold_id0 == 0L) {
        lower <- -0.5
        upper <- if (nf > 1L) {
          seq(1.5, 1.5 + 0.75 * (nf - 2L), length.out = nf - 1L)
        } else {
          numeric(0)
        }

        if (monotone_tau0 == 0L) {
          c(lower, upper)
        } else {
          c(
            log(-lower),
            if (length(upper)) log(diff(c(1, upper))) else numeric(0)
          )
        }
      } else {
        upper <- seq(
          from = 1.5,
          by = 0.75,
          length.out = nf
        )

        if (monotone_tau0 == 0L) {
          upper
        } else {
          log(diff(c(1, upper)))
        }
      }
    }

    if (meas_invar0[["thresholds"]] == 1L) {
      for (k in 1:nvar) {
        nf <- cm3[k]
        if (nf > 0L) {
          b[pos + seq_len(nf)] <- tau_start_values(nf)
          pos <- pos + nf
        }
      }
    } else {
      for (tt in 1:n) {
        for (k in 1:nvar) {
          nf <- cm3[k]
          if (nf > 0L) {
            b[pos + seq_len(nf)] <- tau_start_values(nf)
            pos <- pos + nf
          }
        }
      }
    }

      if (meas_invar0[["loadings"]] == 1L) {
        for (q in 1:Q) {
          ref_item <- marker_idx[q]
          for (k in 1:nvar) {
            if (loadmask0[k,q]==0L) next
            if (k == ref_item) next
            pos <- pos + 1L
            b[pos] <- 0.7
          }
        }
      } else {
        for (tt in 1:n) {
          for (q in 1:Q) {
            ref_item <- marker_idx[q]
            for (k in 1:nvar) {
              if (loadmask0[k,q]==0L) next
              if (k == ref_item) next
              pos <- pos + 1L
              b[pos] <- 0.7
            }
          }
        }
      }

      if (meas_invar0[["intercepts"]] == 1L) {
        for (k in 1:nvar) {
          if (!gamma_free_item[k]) next
          pos <- pos + 1L
          b[pos] <- 0
        }
        for (tt in if (n > 1L) 2:n else integer(0)) for (k in 1:nvar) {
          if (gamma_time_free0[k, tt] != 1L) next
          pos <- pos + 1L
          b[pos] <- 0
        }
      } else {
        for (tt in 1:n) for (k in 1:nvar) {
          if (!gamma_free_item[k]) next
          pos <- pos + 1L
          b[pos] <- 0
        }
      }

      if (meas_invar0[["residuals"]] == 1L) {
        for (k in 1:nvar) {
          if (cm1[k] == 1L) next
          pos <- pos + 1L
          b[pos] <- log(1.0)
        }
      } else {
        for (k in 1:nvar) {
          if (cm1[k] == 1L) next
          for (tt in 1:n) { pos <- pos + 1L; b[pos] <- log(1.0) }
        }
      }

      if (re_item0==1L) {
        for (k in 1:nvar) { pos <- pos + 1L; b[pos] <- log(1.0) }
      }

      for (tt in kappa_times) for (q in 1:Q) {
        pos <- pos + 1L
        b[pos] <- 0
      }

      if (p_fix > 0L) for (q in 1:Q) for (r in 1:p_fix) { pos <- pos + 1L; b[pos] <- 0 }

      for (d in 1:dim_eta) { pos <- pos + 1L; b[pos] <- log(1.0) }
      for (col in 1:dim_eta) if (col < dim_eta) {
        for (row in (col+1):dim_eta) { pos <- pos + 1L; b[pos] <- 0 }
      }


      if (pos != npm0) stop("Bug start: pos != npm0 (pos=", pos, ", npm0=", npm0, ")")
      b
    }

    make_parscale_adaptive <- function(bref) {
      ps <- numeric(npm0)
      pos <- 0L

      clamp_scale <- function(x, smin, smax = 5) {
        pmin(pmax(abs(x), smin), smax)
      }

      if (meas_invar0[["thresholds"]] == 1L) {
        for (k in 1:nvar) {
          nf <- cm3[k]
          if (nf > 0L) {
            idx <- pos + seq_len(nf)
            ps[idx] <- clamp_scale(bref[idx], smin = 0.5, smax = 3)
            pos <- pos + nf
          }
        }
      } else {
        for (tt in 1:n) {
          for (k in 1:nvar) {
            nf <- cm3[k]
            if (nf > 0L) {
              idx <- pos + seq_len(nf)
              ps[idx] <- clamp_scale(bref[idx], smin = 0.5, smax = 3)
              pos <- pos + nf
            }
          }
        }
      }

      if (meas_invar0[["loadings"]] == 1L) {
        for (q in 1:Q) {
          ref_item <- marker_idx[q]
          for (k in 1:nvar) {
            if (loadmask0[k, q] == 0L) next
            if (k == ref_item) next
            pos <- pos + 1L
            ps[pos] <- clamp_scale(bref[pos], smin = 0.3, smax = 2)
          }
        }
      } else {
        for (tt in 1:n) {
          for (q in 1:Q) {
            ref_item <- marker_idx[q]
            for (k in 1:nvar) {
              if (loadmask0[k, q] == 0L) next
              if (k == ref_item) next
              pos <- pos + 1L
              ps[pos] <- clamp_scale(bref[pos], smin = 0.3, smax = 2)
            }
          }
        }
      }

      if (meas_invar0[["intercepts"]] == 1L) {
        for (k in 1:nvar) {
          if (!gamma_free_item[k]) next
          pos <- pos + 1L
          ps[pos] <- clamp_scale(bref[pos], smin = 0.5, smax = 3)
        }
        for (tt in if (n > 1L) 2:n else integer(0)) for (k in 1:nvar) {
          if (gamma_time_free0[k, tt] != 1L) next
          pos <- pos + 1L
          ps[pos] <- clamp_scale(bref[pos], smin = 0.5, smax = 3)
        }
      } else {
        for (tt in 1:n) for (k in 1:nvar) {
          if (!gamma_free_item[k]) next
          pos <- pos + 1L
          ps[pos] <- clamp_scale(bref[pos], smin = 0.5, smax = 3)
        }
      }

      if (meas_invar0[["residuals"]] == 1L) {
        for (k in 1:nvar) {
          if (cm1[k] == 1L) next
          pos <- pos + 1L
          ps[pos] <- clamp_scale(bref[pos], smin = 0.3, smax = 2)
        }
      } else {
        for (k in 1:nvar) {
          if (cm1[k] == 1L) next
          for (tt in 1:n) {
            pos <- pos + 1L
            ps[pos] <- clamp_scale(bref[pos], smin = 0.3, smax = 2)
          }
        }
      }

      if (re_item0 == 1L) {
        for (k in 1:nvar) {
          pos <- pos + 1L
          ps[pos] <- clamp_scale(bref[pos], smin = 0.3, smax = 2)
        }
      }

      for (tt in kappa_times) for (q in 1:Q) {
        pos <- pos + 1L
        ps[pos] <- clamp_scale(bref[pos], smin = 0.25, smax = 2)
      }

      if (p_fix > 0L) for (q in 1:Q) for (r in 1:p_fix) {
        pos <- pos + 1L
        ps[pos] <- clamp_scale(bref[pos], smin = 0.25, smax = 2)
      }

      for (d in 1:dim_eta) {
        pos <- pos + 1L
        ps[pos] <- clamp_scale(bref[pos], smin = 0.3, smax = 2)
      }

      for (col in 1:dim_eta) if (col < dim_eta) {
        for (row in (col + 1):dim_eta) {
          pos <- pos + 1L
          ps[pos] <- clamp_scale(bref[pos], smin = 0.2, smax = 2)
        }
      }


      if (pos != npm0) stop("Bug parscale adaptive: pos != npm0")
      ps
    }

    b0 <- if (is.null(start)) make_start() else {
      start <- as.double(start)
      if (length(start) != npm0) stop("start length != npm0")
      start
    }

    ps0 <- if (is.null(parscale)) {
      make_parscale_adaptive(b0)
    } else {
      parscale <- as.double(parscale)
      if (length(parscale) != npm0) stop("parscale length != npm0")
      if (any(!is.finite(parscale)) || any(parscale <= 0)) {
        stop("parscale deve avere valori finiti > 0")
      }
      parscale
    }

  run_bfgs_two_criteria <- function(par_start,
                                    fn,
                                    maxiter_total,
                                    parscale,
                                    label) {
    outer_max <- max(
      1L,
      floor((sqrt(1 + 8 * maxiter_total) - 1) / 2)
    )

    best <- par_start
    prev_best <- NULL
    prev_value <- NA_real_
    last_value <- NA_real_
    converged <- FALSE

    for (oo in seq_len(outer_max)) {
      opt_inner <- optim(
        par = best,
        fn = fn,
        method = "BFGS",
        control = list(
          maxit = oo,
          trace = if (verbose) 1 else 0,
          parscale = parscale
        )
      )

      best <- opt_inner$par
      last_value <- fn(best)

      if (!is.null(prev_best) &&
          is.finite(prev_value) &&
          is.finite(last_value)) {
        rel_loglik <- abs(last_value - prev_value) /
          max(
            1,
            abs(last_value),
            abs(prev_value)
          )

        rel_par <- max(
          abs(best - prev_best) /
            pmax(
              1,
              abs(best),
              abs(prev_best)
            )
        )

        if (isTRUE(verbose)) {
          cat("outer", oo, "relative likelihood change:", rel_loglik, "\n")
          cat("outer", oo, "relative parameter change:", rel_par, "\n")
        }

        if (rel_loglik < tol_loglik &&
            rel_par < tol_par) {
          converged <- TRUE
          break
        }
      }

      prev_best <- best
      prev_value <- last_value
    }

    list(
      par = best,
      value = last_value,
      convergence = if (isTRUE(converged)) 0L else 1L,
      message = if (isTRUE(converged)) {
        paste0(label, " converged")
      } else {
        paste0(label, " reached maximum iterations")
      },
      counts = setNames(
        c(nfn, NA_integer_),
        c("function", "gradient")
      )
    )
  }

  use_agh <- (methInteg == "AGH")

  fit_id0 <- next_fit_id()

  Y0_c          <- as.integer(Y0)
  nvar0_c       <- as.integer(nvar)
  nsub0_c       <- as.integer(nsub)
  nobs0_c       <- as.integer(nobs)
  n0_c          <- as.integer(n)
  nmes0_c       <- as.integer(nmes0)
  c0_c          <- as.integer(c0)
  cm10_c        <- as.integer(cm10)
  ntre0_c       <- as.integer(ntre0)
  nitem0_c      <- as.integer(nitem0)
  npm0_c        <- as.integer(npm0)
  re_item0_c    <- as.integer(re_item0)
  theta_time0_c <- as.integer(theta_time0)
  meas_invar0_c <- as.integer(meas_invar0)
  Q0_c          <- as.integer(Q)
  loadmask0_c   <- as.integer(loadmask0)
  phi_struct0_c <- 0L
  monotone_tau0_c <- as.integer(monotone_tau0)
  threshold_id0_c <- as.integer(threshold_id0)
  mean_ident0_c <- as.integer(mean_ident0)
  update_gh0_c <- as.integer(update_gh0)
  fit_id0_c     <- as.integer(fit_id0)
  marker0_c     <- as.integer(marker_idx)
  p_fix0_c      <- as.integer(p_fix)
  X0_c          <- as.double(X0)

  set_cfa_long_options(marker0_c, Q0_c, nvar0_c, nsub0_c, n0_c, p_fix0_c, X0_c)
  set_cfa_long_threshold_id(threshold_id0_c)
  set_cfa_long_mean_ident(mean_ident0_c)
  set_cfa_long_gamma_time_free(gamma_time_free0)

  if (!isTRUE(two_tier)) {
    nquad0_c   <- as.integer(nquad0)
    nodes0_c   <- as.double(nodes0)
    weights0_c <- as.double(weights0)
  } else {
    nquad_eta0_c   <- as.integer(nquad_eta0)
    nodes_eta0_c   <- as.double(nodes_eta0)
    weights_eta0_c <- as.double(weights_eta0)

    nquad_re0_c    <- as.integer(nquad_re0)
    nodes_re0_c    <- as.double(nodes_re0)
    weights_re0_c  <- as.double(weights_re0)
  }

   check_symbol <- function(symb) {
  if (!native_has(symb)) {
    stop("Native routine ", symb, " is unavailable.")
  }
}

call_fortran <- function(symb, ...) {
  check_symbol(symb)
  native_call(symb, ..., vrais = double(1))$vrais
}

loglik_long_call <- function(par, use_agh = FALSE, update_gh_eval = 1L) {

  if (!isTRUE(two_tier) && !isTRUE(use_agh)) {
    return(call_fortran(
      "loglik_long",
      b0            = as.double(par),
      Y0            = Y0_c,
      nvar0         = nvar0_c,
      nsub0         = nsub0_c,
      nobs0         = nobs0_c,
      n0            = n0_c,
      nmes0         = nmes0_c,
      c0            = c0_c,
      cm10          = cm10_c,
      ntre0         = ntre0_c,
      nquad0        = nquad0_c,
      nitem0        = nitem0_c,
      nodes0        = nodes0_c,
      weights0      = weights0_c,
      npm0          = npm0_c,
      re_item0      = re_item0_c,
      theta_time0   = theta_time0_c,
      Q0            = Q0_c,
      loadmask0     = loadmask0_c,
      phi_struct0   = phi_struct0_c,
      monotone_tau0 = monotone_tau0_c,
      fit_id0       = fit_id0_c
    ))
  }

  if (!isTRUE(two_tier) && isTRUE(use_agh)) {
    return(call_fortran(
      "loglik_long_agh",
      b0            = as.double(par),
      Y0            = Y0_c,
      nvar0         = nvar0_c,
      nsub0         = nsub0_c,
      nobs0         = nobs0_c,
      n0            = n0_c,
      nmes0         = nmes0_c,
      c0            = c0_c,
      cm10          = cm10_c,
      ntre0         = ntre0_c,
      nquad0        = nquad0_c,
      nitem0        = nitem0_c,
      nodes0        = nodes0_c,
      weights0      = weights0_c,
      npm0          = npm0_c,
      re_item0      = re_item0_c,
      theta_time0   = theta_time0_c,
      Q0            = Q0_c,
      loadmask0     = loadmask0_c,
      phi_struct0   = phi_struct0_c,
      monotone_tau0 = monotone_tau0_c,
      update_gh0    = as.integer(update_gh_eval),
      fit_id0       = fit_id0_c
    ))
  }

if (isTRUE(two_tier) && !isTRUE(use_agh)) {
  return(call_fortran(
    "loglik_long_tt",
    b0            = as.double(par),
    Y0            = Y0_c,
    nvar0         = nvar0_c,
    nsub0         = nsub0_c,
    nobs0         = nobs0_c,
    n0            = n0_c,
    nmes0         = nmes0_c,
    c0            = c0_c,
    cm10          = cm10_c,
    ntre0         = ntre0_c,
    nquad_eta0    = nquad_eta0_c,
    nodes_eta0    = nodes_eta0_c,
    weights_eta0  = weights_eta0_c,
    nquad_re0     = nquad_re0_c,
    nodes_re0     = nodes_re0_c,
    weights_re0   = weights_re0_c,
    nitem0        = nitem0_c,
    npm0          = npm0_c,
    re_item0      = re_item0_c,
    theta_time0   = theta_time0_c,
    Q0            = Q0_c,
    loadmask0     = loadmask0_c,
    phi_struct0   = phi_struct0_c,
    monotone_tau0 = monotone_tau0_c,
    fit_id0       = fit_id0_c
  ))
}
  tt_agh_symbol <- if (identical(two_tier_type, "full_equiv")) {
    "loglik_long_tt_agh_full_equiv"
  } else {
    "loglik_long_tt_agh"
  }

  if (identical(tt_agh_symbol, "loglik_long_tt_agh_full_equiv")) {
    return(call_fortran(
      tt_agh_symbol,
      b0            = as.double(par),
      Y0            = Y0_c,
      nvar0         = nvar0_c,
      nsub0         = nsub0_c,
      nobs0         = nobs0_c,
      n0            = n0_c,
      nmes0         = nmes0_c,
      c0            = c0_c,
      cm10          = cm10_c,
      ntre0         = ntre0_c,
      nquad_eta0    = nquad_eta0_c,
      nitem0        = nitem0_c,
      nodes_eta0    = nodes_eta0_c,
      weights_eta0  = weights_eta0_c,
      nquad_re0     = nquad_re0_c,
      nodes_re0     = nodes_re0_c,
      weights_re0   = weights_re0_c,
      npm0          = npm0_c,
      re_item0      = re_item0_c,
      theta_time0   = theta_time0_c,
      meas_invar0   = meas_invar0_c,
      Q0            = Q0_c,
      loadmask0     = loadmask0_c,
      phi_struct0   = phi_struct0_c,
      monotone_tau0 = monotone_tau0_c,
      update_gh0    = as.integer(update_gh_eval),
      fit_id0       = fit_id0_c
    ))
  }

  call_fortran(
    tt_agh_symbol,
    b0            = as.double(par),
    Y0            = Y0_c,
    nvar0         = nvar0_c,
    nsub0         = nsub0_c,
    nobs0         = nobs0_c,
    n0            = n0_c,
    nmes0         = nmes0_c,
    c0            = c0_c,
    cm10          = cm10_c,
    ntre0         = ntre0_c,
    nquad_eta0    = nquad_eta0_c,
    nitem0        = nitem0_c,
    nodes_eta0    = nodes_eta0_c,
    weights_eta0  = weights_eta0_c,
    nquad_re0     = nquad_re0_c,
    nodes_re0     = nodes_re0_c,
    weights_re0   = weights_re0_c,
    npm0          = npm0_c,
    re_item0      = re_item0_c,
    theta_time0   = theta_time0_c,
    Q0            = Q0_c,
    loadmask0     = loadmask0_c,
    phi_struct0   = phi_struct0_c,
    monotone_tau0 = monotone_tau0_c,
    update_gh0    = as.integer(update_gh_eval),
    fit_id0       = fit_id0_c
  )
}

  score_symbol <- "loglik_long_tt_agh_full_equiv_score_q1"
  score_compatible <- isTRUE(analytic_score) &&
    isTRUE(use_agh) &&
    isTRUE(two_tier) &&
    identical(two_tier_type, "full_equiv") &&
    Q == 1L &&
    meas_invar0[["thresholds"]] == 1L &&
    meas_invar0[["loadings"]] == 1L &&
    !is.null(update_gh) &&
    update_gh > 0L

  use_analytic_score <- score_compatible && native_has(score_symbol)

  if (isTRUE(analytic_score) && !score_compatible && isTRUE(verbose)) {
    cat("analytic score not used: requires Q=1, AGH full_equiv, FREE Phi, invariant thresholds/loadings, update_gh>0\n")
  }

  if (score_compatible && !native_has(score_symbol)) {
    warning("analytic score requested but symbol ", score_symbol,
            " is not loaded; falling back to numerical BFGS gradient")
  }

  call_fortran_score_fixed <- function(par) {
    out <- native_call(
      score_symbol,
      b0            = as.double(par),
      Y0            = Y0_c,
      nvar0         = nvar0_c,
      nsub0         = nsub0_c,
      nobs0         = nobs0_c,
      n0            = n0_c,
      nmes0         = nmes0_c,
      c0            = c0_c,
      cm10          = cm10_c,
      ntre0         = ntre0_c,
      nquad_eta0    = nquad_eta0_c,
      nitem0        = nitem0_c,
      nodes_eta0    = nodes_eta0_c,
      weights_eta0  = weights_eta0_c,
      nquad_re0     = nquad_re0_c,
      nodes_re0     = nodes_re0_c,
      weights_re0   = weights_re0_c,
      npm0          = npm0_c,
      re_item0      = re_item0_c,
      theta_time0   = theta_time0_c,
      meas_invar0   = meas_invar0_c,
      Q0            = Q0_c,
      loadmask0     = loadmask0_c,
      phi_struct0   = phi_struct0_c,
      monotone_tau0 = monotone_tau0_c,
      update_gh0    = 0L,
      fit_id0       = fit_id0_c,
      vrais         = double(1),
      score0        = double(npm0),
      ok0           = integer(1)
    )

    list(
      value = out$vrais,
      gradient = out$score0,
      ok = identical(out$ok0, 1L) &&
        is.finite(out$vrais) &&
        all(is.finite(out$score0))
    )
  }

  nfn <- 0L
  ngrad <- 0L
  fn_elapsed <- numeric(0)

  fn <- function(par) {
  t0_fn <- proc.time()[["elapsed"]]

  val <- loglik_long_call(
    par = par,
    use_agh = use_agh
  )

  if (isTRUE(profile)) {
    dt <- proc.time()[["elapsed"]] - t0_fn
    nfn <<- nfn + 1L
    fn_elapsed <<- c(fn_elapsed, dt)
  }

  val
}

if (maxiter == 0L) {
  t0_eval <- proc.time()

  val0 <- loglik_long_call(
    par = b0,
    use_agh = use_agh,
    update_gh_eval = 1L
  )

  t_eval <- proc.time() - t0_eval

  if (!valid_objective(val0)) stop("The initial CFA likelihood could not be evaluated.")
  return(new_long_fit(list(
    conv = NA_integer_,
    message = "single likelihood evaluation only",
    loglik = -val0,
    negloglik = val0,
    par_opt = b0,
    parscale = ps0,
    timing = list(
      eval_elapsed = unname(t_eval["elapsed"]),
      eval_user = unname(t_eval["user.self"]),
      eval_system = unname(t_eval["sys.self"])
    ),
    dims = list(
      nsub = nsub,
      nvar = nvar,
      n = n,
      Q = Q,
      dim_eta = dim_eta,
      p_fix = p_fix
    ),
    structure = list(
      loadmask = loadmask0,
      marker = setNames(namvar[marker_idx], paste0("F", seq_len(Q))),
      covariates = covariates,
      post_covariates = post_covariates,
      cov_time_mask = cov_time_mask,
      methInteg = methInteg,
      re_item = re_item0,
      two_tier = isTRUE(two_tier),
      two_tier_type = if (isTRUE(two_tier)) two_tier_type else NA_character_,
      analytic_score = isTRUE(use_analytic_score),
      threshold_id = threshold_id,
      threshold_id_code = threshold_id0,
      mean_ident = mean_ident,
      mean_ident_code = mean_ident0,
      threshold_fixed = if (threshold_id0 == 0L) {
        c(tau2 = 0, tau3 = 1)
      } else {
        c(tau1 = 0, tau2 = 1)
      },
      theta_time = theta_time0,
      meas_invar = meas_invar0,
      intercepts_free_items = intercepts_free_items,
      gamma_time_free = gamma_time_free0,
      phi_struct = phi_struct,
      nquad = if (!isTRUE(two_tier)) nquad0 else NULL,
      nquad_eta = if (isTRUE(two_tier)) nquad_eta0 else NULL,
      nquad_re = if (isTRUE(two_tier)) nquad_re0 else NULL,
      dimInt = if (!isTRUE(two_tier)) dimInt else NULL,
      dim_eta_int = if (isTRUE(two_tier)) dim_eta_int else NULL
    )
  ), "Longitudinal ordinal CFA", packed))
}
agh_fast <- isTRUE(use_agh) && !is.null(update_gh) && update_gh > 0L

if (!isTRUE(use_agh)) {

  if (verbose) {
    cat("\nstandard optimization with", methInteg, "\n")
  }

  fn <- function(par) {
    t0_fn <- proc.time()[["elapsed"]]

    val <- loglik_long_call(
      par = par,
      use_agh = FALSE
    )

    if (isTRUE(profile)) {
      dt <- proc.time()[["elapsed"]] - t0_fn
      nfn <<- nfn + 1L
      fn_elapsed <<- c(fn_elapsed, dt)
    }

    val
  }

  t0 <- proc.time()

  opt <- optim(
  par = b0,
  fn = fn,
  method = "BFGS",
  control = list(
    maxit = maxiter,
    trace = if (verbose) 1 else 0,
    parscale = ps0
  )
)

  t1 <- proc.time() - t0
  best <- opt$par
  nfn_opt <- nfn
  fn_hess <- fn

} else if (!agh_fast) {

  if (verbose) {
    cat("\nfull AGH optimization: updating modes and adaptive nodes at every likelihood evaluation\n")
  }

  fn <- function(par) {
    t0_fn <- proc.time()[["elapsed"]]

    val <- loglik_long_call(
      par = par,
      use_agh = TRUE,
      update_gh_eval = 1L
    )

    if (isTRUE(profile)) {
      dt <- proc.time()[["elapsed"]] - t0_fn
      nfn <<- nfn + 1L
      fn_elapsed <<- c(fn_elapsed, dt)
    }

    val
  }

  t0 <- proc.time()

  opt <- run_bfgs_two_criteria(
    par_start = b0,
    fn = fn,
    maxiter_total = maxiter,
    parscale = ps0,
    label = "full AGH optimization"
  )

  t1 <- proc.time() - t0
  best <- opt$par
  nfn_opt <- nfn
  fn_hess <- fn

} else {

  if (verbose) {
    cat("\nfast AGH optimization: outer-inner scheme\n")
  }

  agh_step <- as.integer(update_gh)
  agh_outer <- max(
    1L,
    floor((sqrt(1 + 8 * maxiter / agh_step) - 1) / 2)
  )

  best <- b0

  agh_trace <- data.frame(
    outer = integer(),
    inner_maxit = integer(),
    negloglik_grid_start = numeric(),
    negloglik_fixed_end = numeric(),
    negloglik_refreshed = numeric(),
    loglik_refreshed = numeric(),
    convergence = integer(),
    nfn = integer()
  )

  t0 <- proc.time()
  prev_refreshed <- NA_real_
  prev_best <- NULL
  agh_converged <- FALSE

  agh_failed <- FALSE
  agh_fail_message <- NULL

  for (oo in seq_len(agh_outer)) {

    inner_maxit <- oo * agh_step

    if (verbose) {
      cat("\nAGH outer:", oo, "\n")
      cat("updating AGH grid...\n")
    }

    val_grid <- eval_ll(best, TRUE, 1L)

    if (!is.finite(val_grid) || val_grid >= 1e8) {
      agh_failed <- TRUE
      agh_fail_message <- paste0("AGH grid update failed at outer ", oo)
      break
    }

    if (verbose) {
      cat("initial fixed-grid negloglik:", val_grid, "\n")
      cat("optimizing fixed AGH grid for maxit =", inner_maxit, "\n")
    }

    score_cache <- new.env(parent = emptyenv())
    score_cache$par <- NULL
    score_cache$result <- NULL

    numeric_fixed_gradient <- function(par) {
      g <- numeric(length(par))
      step <- 1e-6 * pmax(1, abs(par))
      for (jj in seq_along(par)) {
        p1 <- par
        p2 <- par
        p1[jj] <- p1[jj] + step[jj]
        p2[jj] <- p2[jj] - step[jj]
        f1 <- loglik_long_call(p1, use_agh = TRUE, update_gh_eval = 0L)
        f2 <- loglik_long_call(p2, use_agh = TRUE, update_gh_eval = 0L)
        g[jj] <- (f1 - f2) / (2 * step[jj])
      }
      g
    }

    eval_fixed_score <- function(par) {
      p <- as.double(par)
      if (!is.null(score_cache$par) && identical(p, score_cache$par)) {
        return(score_cache$result)
      }

      ans <- call_fortran_score_fixed(p)
      if (!isTRUE(ans$ok)) {
        ans$value <- loglik_long_call(
          p,
          use_agh = TRUE,
          update_gh_eval = 0L
        )
        ans$gradient <- NULL
      }

      score_cache$par <- p
      score_cache$result <- ans
      ans
    }

    fn_fixed <- function(par) {
      t0_fn <- proc.time()[["elapsed"]]

      val <- if (isTRUE(use_analytic_score)) {
        eval_fixed_score(par)$value
      } else {
        loglik_long_call(
          par = par,
          use_agh = TRUE,
          update_gh_eval = 0L
        )
      }

      if (isTRUE(profile)) {
        dt <- proc.time()[["elapsed"]] - t0_fn
        nfn <<- nfn + 1L
        fn_elapsed <<- c(fn_elapsed, dt)
      }

      val
    }

    gr_fixed <- NULL
    if (isTRUE(use_analytic_score)) {
      gr_fixed <- function(par) {
        ngrad <<- ngrad + 1L
        ans <- eval_fixed_score(par)
        if (!is.null(ans$gradient)) return(ans$gradient)
        numeric_fixed_gradient(par)
      }
    }

    if (isTRUE(use_analytic_score) && isTRUE(verbose)) {
      cat("using analytic fixed-grid AGH score in BFGS\n")
    }

    opt_inner <- tryCatch(
      optim(
        par = best,
        fn = fn_fixed,
        gr = gr_fixed,
        method = "BFGS",
        control = list(
          maxit = inner_maxit,
          trace = if (verbose) 1 else 0,
          parscale = ps0
        )
      ),
      error = function(e) e
    )

    if (inherits(opt_inner, "error")) {
      agh_failed <- TRUE
      agh_fail_message <- paste0("AGH fixed-grid optimization failed at outer ", oo, ": ", opt_inner$message)
      break
    }

    best <- opt_inner$par

    val_fixed <- eval_ll(best, TRUE, 0L)
    val_refreshed <- eval_ll(best, TRUE, 1L)

    if (!is.finite(val_refreshed) || val_refreshed >= 1e8) {
      agh_failed <- TRUE
      agh_fail_message <- paste0("final AGH grid refresh failed at outer ", oo)
      break
    }

    agh_trace <- rbind(
      agh_trace,
      data.frame(
        outer = oo,
        inner_maxit = inner_maxit,
        negloglik_grid_start = val_grid,
        negloglik_fixed_end = val_fixed,
        negloglik_refreshed = val_refreshed,
        loglik_refreshed = -val_refreshed,
        convergence = opt_inner$convergence,
        nfn = nfn
      )
    )

    if (verbose) {
      cat("outer", oo, "inner maxit:", inner_maxit, "\n")
      cat("outer", oo, "fixed-grid negloglik:", val_fixed, "\n")
      cat("outer", oo, "refreshed AGH negloglik:", val_refreshed, "\n")
      cat("outer", oo, "refreshed AGH loglik:", -val_refreshed, "\n")
    }

    if (is.finite(prev_refreshed) &&
        is.finite(val_refreshed) &&
        !is.null(prev_best)) {
      rel_loglik <- abs(val_refreshed - prev_refreshed) /
        max(
          1,
          abs(val_refreshed),
          abs(prev_refreshed)
        )

      rel_par <- max(
        abs(best - prev_best) /
          pmax(
            1,
            abs(best),
            abs(prev_best)
          )
      )

      if (isTRUE(verbose)) {
        cat("outer", oo, "relative likelihood change:", rel_loglik, "\n")
        cat("outer", oo, "relative parameter change:", rel_par, "\n")
      }

      if (rel_loglik < tol_loglik &&
          rel_par < tol_par) {
        agh_converged <- TRUE
        break
      }
    }

    prev_refreshed <- val_refreshed
    prev_best <- best
  }

 if (isTRUE(agh_failed) || nrow(agh_trace) == 0L) {
  stop(if (!is.null(agh_fail_message)) agh_fail_message else "AGH fast failed before producing a valid iteration.")
}

final_value <- tail(agh_trace$negloglik_refreshed, 1)

  t1 <- proc.time() - t0

  opt <- list(
    par = best,
    value = final_value,
    convergence = if (isTRUE(agh_converged)) 0L else 1L,
    message = if (isTRUE(agh_converged)) {
      "AGH fast outer-inner optimization converged"
    } else {
      "AGH fast outer-inner optimization reached maximum iterations"
    },
    counts = setNames(
      c(nfn, if (isTRUE(use_analytic_score)) ngrad else NA_integer_),
      c("function", "gradient")
    )
  )

  nfn_opt <- nfn

  fn_hess <- function(par) {
    loglik_long_call(
      par = par,
      use_agh = TRUE,
      update_gh_eval = 1L
    )
  }
}

best <- opt$par
nfn_opt <- nfn

  hess <- NULL
  hessian_raw <- NULL
  hessian_step <- NULL
  score_at_optimum <- NULL
  hessian_eigenvalues <- NULL
  hessian_positive_definite <- NA
  hessian_method <- "none"
  ngrad_hess <- 0L
  t3 <- proc.time() - proc.time()
  nfn_hess <- 0L

  if (isTRUE(compute_hessian)) {
    t2 <- proc.time()

    if (isTRUE(use_analytic_score)) {
      hessian_method <- "analytic_fixed"

      val_final_refresh <- loglik_long_call(
        par = best,
        use_agh = TRUE,
        update_gh_eval = 1L
      )
      if (!is.finite(val_final_refresh) || val_final_refresh >= 1e8) {
        stop("refresh finale AGH fallito prima della Hessiana analytic_fixed")
      }
      opt$value <- val_final_refresh

      e0 <- call_fortran_score_fixed(best)
      if (!isTRUE(e0$ok)) {
        stop("score analitico fallito alla soluzione ottima")
      }
      score_at_optimum <- e0$gradient
      ngrad_hess <- ngrad_hess + 1L

      hessian_step <- hess_eps * (abs(best) + 1)
      hessian_raw <- matrix(NA_real_, npm0, npm0)

      for (jj in seq_len(npm0)) {
        pp <- best
        pm <- best
        pp[jj] <- pp[jj] + hessian_step[jj]
        pm[jj] <- pm[jj] - hessian_step[jj]

        sp <- call_fortran_score_fixed(pp)
        sm <- call_fortran_score_fixed(pm)
        ngrad_hess <- ngrad_hess + 2L

        if (!isTRUE(sp$ok) || !isTRUE(sm$ok)) {
          stop("score analitico fallito durante la Hessiana al parametro ", jj)
        }

        hessian_raw[, jj] <-
          (sp$gradient - sm$gradient) / (2 * hessian_step[jj])
      }

      hess <- (hessian_raw + t(hessian_raw)) / 2
      hessian_eigenvalues <- eigen(
        hess, symmetric = TRUE, only.values = TRUE
      )$values
      hessian_positive_definite <-
        all(is.finite(hessian_eigenvalues)) &&
        min(hessian_eigenvalues) > 0

      if (!isTRUE(hessian_positive_definite)) {
        warning(
          "Hessiana analytic_fixed non positiva definita; SE/Wald non affidabili. min eigenvalue=",
          format(min(hessian_eigenvalues), scientific = TRUE)
        )
      }
    } else {
      hessian_method <- "numeric_refresh"
      hess <- optimHess(par = best, fn = fn_hess)
      hess <- (hess + t(hess)) / 2
      hessian_eigenvalues <- eigen(
        hess, symmetric = TRUE, only.values = TRUE
      )$values
      hessian_positive_definite <-
        all(is.finite(hessian_eigenvalues)) &&
        min(hessian_eigenvalues) > 0
      nfn_hess <- nfn - nfn_opt
    }

    t3 <- proc.time() - t2
  }
  timing <- list(
    opt_elapsed  = unname(t1["elapsed"]),
    opt_user     = unname(t1["user.self"]),
    opt_system   = unname(t1["sys.self"]),
    hess_elapsed = unname(t3["elapsed"]),
    hess_user    = unname(t3["user.self"]),
    hess_system  = unname(t3["sys.self"]),
    total_elapsed = unname(t1["elapsed"]) + unname(t3["elapsed"])
  )

  profile_info <- NULL
  if (isTRUE(profile)) {
    profile_info <- list(
      npar = length(b0),
      nfn_total = nfn,
      nfn_opt = nfn_opt,
      nfn_hess = nfn_hess,
      ngrad_hess = ngrad_hess,
      hessian_method = hessian_method,
      fn_mean = mean(fn_elapsed),
      fn_median = median(fn_elapsed),
      fn_min = min(fn_elapsed),
      fn_max = max(fn_elapsed),
      fn_sum = sum(fn_elapsed),
      parscale_summary = c(
      min = min(ps0),
      median = median(ps0),
      mean = mean(ps0),
      max = max(ps0)
    )
    )

        if (isTRUE(profile_print)) {
      cat("\n")
      cat("profiling summary\n")
      cat("nsub:", nsub, "\n")
      cat("npar:", profile_info$npar, "\n")
      cat("two_tier:", isTRUE(two_tier), "\n")
      if (isTRUE(two_tier)) cat("two_tier_type:", two_tier_type, "\n")
      if (!isTRUE(two_tier)) {
        cat("dimInt:", dimInt, "\n")
        cat("nquad:", nquad0, "\n")
      } else {
        cat("dim_eta_int:", dim_eta_int, "\n")
        cat("nquad_eta:", nquad_eta0, "\n")
        cat("nquad_re:", nquad_re0, "\n")
      }
      cat("nfn_total:", profile_info$nfn_total, "\n")
      cat("nfn_opt:", profile_info$nfn_opt, "\n")
      cat("nfn_hess:", profile_info$nfn_hess, "\n")
      cat("ngrad_hess:", profile_info$ngrad_hess, "\n")
      cat("hessian_method:", profile_info$hessian_method, "\n")
      cat("fn_mean:", sprintf("%.6f", profile_info$fn_mean), "\n")
      cat("fn_median:", sprintf("%.6f", profile_info$fn_median), "\n")
      cat("fn_min:", sprintf("%.6f", profile_info$fn_min), "\n")
      cat("fn_max:", sprintf("%.6f", profile_info$fn_max), "\n")
      cat("fn_sum:", sprintf("%.6f", profile_info$fn_sum), "\n")
      cat("optim_elapsed:", sprintf("%.6f", timing$opt_elapsed), "\n")
      cat("hess_elapsed:", sprintf("%.6f", timing$hess_elapsed), "\n")
      cat("total_elapsed:", sprintf("%.6f", timing$total_elapsed), "\n")
      cat("\n")
      flush.console()
    }
  }
  inv_sym <- function(M) {
    out <- tryCatch(solve(M), error = function(e) NULL)
    if (is.null(out)) return(NULL)
    out <- (out + t(out)) / 2
    out
  }

  vcov_opt <- if (!is.null(hess)) {
    if (identical(hessian_method, "analytic_fixed") &&
        !isTRUE(hessian_positive_definite)) {
      NULL
    } else {
      inv_sym(hess)
    }
  } else {
    NULL
  }

  factor_names <- paste0("F", seq_len(Q))
  time_names   <- paste0("t", seq_len(n))
  eta_labels   <- unlist(lapply(time_names, function(tt) paste0(tt, ":", factor_names)))

  nm_opt <- character(npm0)
  pos_nm <- 0L

  if (meas_invar0[["thresholds"]] == 1L) {
    for (k in seq_len(nvar)) {
      if (cm3[k] > 0L) for (jj in seq_len(cm3[k])) {
        pos_nm <- pos_nm + 1L
        cc_free <- if (threshold_id0 == 0L) {
          if (jj == 1L) 1L else jj + 2L
        } else {
          jj + 2L
        }
        nm_opt[pos_nm] <- sprintf("tau_free[%s]_c%d", namvar[k], cc_free)
      }
    }
  } else {
    for (tt in seq_len(n)) {
      for (k in seq_len(nvar)) {
        if (cm3[k] > 0L) for (jj in seq_len(cm3[k])) {
          pos_nm <- pos_nm + 1L
          cc_free <- if (threshold_id0 == 0L) {
            if (jj == 1L) 1L else jj + 2L
          } else {
            jj + 2L
          }
          nm_opt[pos_nm] <- sprintf("tau_free[%s,%s]_c%d", namvar[k], time_names[tt], cc_free)
        }
      }
    }
  }

  if (meas_invar0[["loadings"]] == 1L) {
    for (qq in seq_len(Q)) {
      ref_item <- marker_idx[qq]
      for (k in seq_len(nvar)) {
        if (loadmask0[k, qq] == 0L) next
        if (k == ref_item) next
        pos_nm <- pos_nm + 1L
        nm_opt[pos_nm] <- sprintf("lambda_free[%s,%s]", namvar[k], factor_names[qq])
      }
    }
  } else {
    for (tt in seq_len(n)) {
      for (qq in seq_len(Q)) {
        ref_item <- marker_idx[qq]
        for (k in seq_len(nvar)) {
          if (loadmask0[k, qq] == 0L) next
          if (k == ref_item) next
          pos_nm <- pos_nm + 1L
          nm_opt[pos_nm] <- sprintf("lambda_free[%s,%s,%s]", namvar[k], factor_names[qq], time_names[tt])
        }
      }
    }
  }

  if (meas_invar0[["intercepts"]] == 1L) {
    for (k in seq_len(nvar)) {
      if (!gamma_free_item[k]) next
      pos_nm <- pos_nm + 1L
      nm_opt[pos_nm] <- sprintf("gamma[%s]", namvar[k])
    }
    for (tt in if (n > 1L) 2:n else integer(0)) for (k in seq_len(nvar)) {
      if (gamma_time_free0[k, tt] != 1L) next
      pos_nm <- pos_nm + 1L
      nm_opt[pos_nm] <- sprintf("gamma[%s,%s]", namvar[k], time_names[tt])
    }
  } else {
    for (tt in seq_len(n)) for (k in seq_len(nvar)) {
      if (!gamma_free_item[k]) next
      pos_nm <- pos_nm + 1L
      nm_opt[pos_nm] <- sprintf("gamma[%s,%s]", namvar[k], time_names[tt])
    }
  }

  if (meas_invar0[["residuals"]] == 1L) {
    for (k in seq_len(nvar)) {
      if (cm1[k] == 1L) next
      pos_nm <- pos_nm + 1L
      nm_opt[pos_nm] <- sprintf("log_theta[%s]", namvar[k])
    }
  } else {
    for (k in seq_len(nvar)) {
      if (cm1[k] == 1L) next
      for (tt in seq_len(n)) {
        pos_nm <- pos_nm + 1L
        nm_opt[pos_nm] <- sprintf("log_theta[%s,%s]", namvar[k], time_names[tt])
      }
    }
  }

  if (re_item0 == 1L) {
    for (k in seq_len(nvar)) {
      pos_nm <- pos_nm + 1L
      nm_opt[pos_nm] <- sprintf("log_sigmak[%s]", namvar[k])
    }
  }

  for (tt in kappa_times) for (qq in seq_len(Q)) {
    pos_nm <- pos_nm + 1L
    nm_opt[pos_nm] <- sprintf("kappa[%s,%s]", time_names[tt], factor_names[qq])
  }

  if (p_fix > 0L) {
    for (qq in seq_len(Q)) for (r in seq_len(p_fix)) {
      pos_nm <- pos_nm + 1L
      nm_opt[pos_nm] <- sprintf("beta[%s,%s]", factor_names[qq], covariates[r])
    }
  }

  for (d in seq_len(dim_eta)) {
    pos_nm <- pos_nm + 1L
    nm_opt[pos_nm] <- sprintf("log_Lphi[%s]", eta_labels[d])
  }
  for (col in 1:dim_eta) if (col < dim_eta) for (row in (col + 1L):dim_eta) {
    pos_nm <- pos_nm + 1L
    nm_opt[pos_nm] <- sprintf("Lphi[%s,%s]", eta_labels[row], eta_labels[col])
  }


  if (pos_nm != npm0) stop("Naming: pos != npm0 (pos=", pos_nm, ", npm0=", npm0, "). Blocchi disallineati.")

  names(best) <- nm_opt
  par_opt <- best

  pos <- 0L

  ta <- if (ntre0 > 0L) best[1:ntre0] else numeric(0)
  pos <- ntre0

  tau_time <- vector("list", nvar)
  names(tau_time) <- namvar
  jj <- 0L

  build_tau <- function(k, u) {
    m <- cm1[k]
    nf <- cm3[k]
    thr <- numeric(m)

    if (threshold_id0 == 0L) {
      thr[2] <- 0
      thr[3] <- 1

      if (monotone_tau0 == 0L) {
        thr[1] <- u[1]
        if (m > 3L) thr[4:m] <- u[2:nf]
      } else {
        thr[1] <- -exp(u[1])
        if (m > 3L) thr[4:m] <- 1 + cumsum(exp(u[2:nf]))
      }
    } else {
      thr[1] <- 0
      thr[2] <- 1

      if (monotone_tau0 == 0L) {
        thr[3:m] <- u[seq_len(nf)]
      } else {
        thr[3:m] <- 1 + cumsum(exp(u[seq_len(nf)]))
      }
    }

    thr
  }

  if (meas_invar0[["thresholds"]] == 1L) {
    for (k in 1:nvar) {
      nf <- cm3[k]
      u <- if (nf > 0L) ta[(jj + 1L):(jj + nf)] else numeric(0)
      thr <- build_tau(k, u)
      tau_time[[k]] <- matrix(rep(thr, n), nrow = length(thr), ncol = n,
                              dimnames = list(paste0("c", seq_along(thr)), time_names))
      jj <- jj + nf
    }
  } else {
    tmp <- vector("list", nvar)
    for (k in 1:nvar) tmp[[k]] <- matrix(NA_real_, nrow = cm1[k], ncol = n,
                                        dimnames = list(paste0("c", seq_len(cm1[k])), time_names))
    for (tt in 1:n) {
      for (k in 1:nvar) {
        nf <- cm3[k]
        u <- if (nf > 0L) ta[(jj + 1L):(jj + nf)] else numeric(0)
        tmp[[k]][, tt] <- build_tau(k, u)
        jj <- jj + nf
      }
    }
    tau_time <- tmp
    names(tau_time) <- namvar
  }

  tau_full <- lapply(tau_time, function(x) x[, 1])
  names(tau_full) <- namvar

  lambda_time <- array(0, dim = c(nvar, Q, n), dimnames = list(namvar, factor_names, time_names))
  if (meas_invar0[["loadings"]] == 1L) {
    lambda_base <- matrix(0, nvar, Q, dimnames = list(namvar, factor_names))
    for (qq in 1:Q) {
      ref_item <- marker_idx[qq]
      lambda_base[ref_item, qq] <- 1
      for (k in 1:nvar) {
        if (loadmask0[k, qq] == 0L) next
        if (k == ref_item) next
        pos <- pos + 1L
        lambda_base[k, qq] <- best[pos]
      }
    }
    for (tt in 1:n) lambda_time[, , tt] <- lambda_base
  } else {
    for (tt in 1:n) {
      for (qq in 1:Q) {
        ref_item <- marker_idx[qq]
        lambda_time[ref_item, qq, tt] <- 1
        for (k in 1:nvar) {
          if (loadmask0[k, qq] == 0L) next
          if (k == ref_item) next
          pos <- pos + 1L
          lambda_time[k, qq, tt] <- best[pos]
        }
      }
    }
  }
  lambda_hat <- matrix(
  lambda_time[, , 1, drop = FALSE],
  nrow = nvar,
  ncol = Q,
  dimnames = list(namvar, factor_names)
)

  gamma_time <- matrix(0, nvar, n, dimnames = list(namvar, time_names))
  gamma_hat <- setNames(numeric(nvar), namvar)
  if (meas_invar0[["intercepts"]] == 1L) {
    for (k in seq_len(nvar)) {
      if (!gamma_free_item[k]) next
      pos <- pos + 1L
      gamma_hat[k] <- best[pos]
    }
    for (tt in 1:n) gamma_time[, tt] <- gamma_hat
    for (tt in if (n > 1L) 2:n else integer(0)) for (k in 1:nvar) {
      if (gamma_time_free0[k, tt] != 1L) next
      pos <- pos + 1L
      gamma_time[k, tt] <- best[pos]
    }
  } else {
    for (tt in 1:n) {
      for (k in 1:nvar) {
        if (!gamma_free_item[k]) next
        pos <- pos + 1L
        gamma_time[k, tt] <- best[pos]
      }
    }
    gamma_hat <- gamma_time[, 1]
    names(gamma_hat) <- namvar
  }

  theta_hat <- matrix(NA_real_, nvar, n, dimnames = list(namvar, time_names))
  if (meas_invar0[["residuals"]] == 1L) {
    for (k in 1:nvar) {
      if (cm1[k] == 1L) {
        theta_hat[k, ] <- 1
      } else {
        pos <- pos + 1L
        theta_hat[k, ] <- exp(best[pos])
      }
    }
  } else {
    for (k in 1:nvar) {
      if (cm1[k] == 1L) {
        theta_hat[k, ] <- 1
      } else {
        for (tt in 1:n) {
          pos <- pos + 1L
          theta_hat[k, tt] <- exp(best[pos])
        }
      }
    }
  }

  sigmak_hat <- rep(0, nvar); names(sigmak_hat) <- namvar
  if (re_item0 == 1L) {
    for (k in 1:nvar) {
      pos <- pos + 1L
      sigmak_hat[k] <- exp(best[pos])
    }
  }

  kappa_hat <- matrix(0, n, Q, dimnames = list(time_names, factor_names))
  for (tt in kappa_times) for (qq in 1:Q) {
    pos <- pos + 1L
    kappa_hat[tt, qq] <- best[pos]
  }

  beta_hat <- matrix(0, nrow = Q, ncol = p_fix, dimnames = list(factor_names, covariates))
  if (p_fix > 0L) {
    for (qq in seq_len(Q)) for (r in seq_len(p_fix)) {
      pos <- pos + 1L
      beta_hat[qq, r] <- best[pos]
    }
  }

  Lphi <- matrix(0, dim_eta, dim_eta, dimnames = list(eta_labels, eta_labels))

  diag(Lphi) <- exp(best[(pos + 1L):(pos + dim_eta)])
  pos <- pos + dim_eta
  for (col in 1:dim_eta) if (col < dim_eta) for (row in (col + 1L):dim_eta) {
    pos <- pos + 1L
    Lphi[row, col] <- best[pos]
  }


  Phi <- Lphi %*% t(Lphi)

  if (pos != npm0) stop("Decode: pos != npm0 (pos=", pos, ", npm0=", npm0, ")")

  par_interest <- par_opt

  idx_log_theta <- grep("^log_theta\\[", names(par_interest))
  if (length(idx_log_theta)) {
    par_interest[idx_log_theta] <- exp(par_interest[idx_log_theta])
    names(par_interest)[idx_log_theta] <- sub("^log_theta\\[", "theta[", names(par_interest)[idx_log_theta])
  }

  idx_log_sig <- grep("^log_sigmak\\[", names(par_interest))
  if (length(idx_log_sig)) {
    par_interest[idx_log_sig] <- exp(par_interest[idx_log_sig])
    names(par_interest)[idx_log_sig] <- sub("^log_sigmak\\[", "Var_itemRE[", names(par_interest)[idx_log_sig])
  }

  idx_Lphi <- grep("^(log_Lphi\\[|Lphi\\[)", names(par_interest))
  npar_Lphi <- dim_eta * (dim_eta + 1L) / 2L
  if (length(idx_Lphi) != npar_Lphi) stop("Blocco L_phi non ha dimensione attesa.")

  tmp <- numeric(npar_Lphi)
  tmp_names <- character(npar_Lphi)
  tpos <- 0L

  for (d in 1:dim_eta) {
    tpos <- tpos + 1L
    tmp[tpos] <- Phi[d, d]
    tmp_names[tpos] <- sprintf("Var_eta[%s]", eta_labels[d])
  }
  for (col in 1:dim_eta) if (col < dim_eta) for (row in (col + 1L):dim_eta) {
    tpos <- tpos + 1L
    tmp[tpos] <- Phi[row, col]
    tmp_names[tpos] <- sprintf("Cov_eta[%s,%s]", eta_labels[row], eta_labels[col])
  }

  par_interest[idx_Lphi] <- tmp
  names(par_interest)[idx_Lphi] <- tmp_names


  par_nat <- numeric(0)
  nm_nat  <- character(0)

  if (meas_invar0[["thresholds"]] == 1L) {
    for (k in seq_len(nvar)) {
      tk <- tau_time[[k]][, 1]
      for (cc in seq_along(tk)) {
        nm_nat  <- c(nm_nat, sprintf("tau[%s]_c%d", namvar[k], cc))
        par_nat <- c(par_nat, tk[cc])
      }
    }
  } else {
    for (tt in seq_len(n)) for (k in seq_len(nvar)) {
      tk <- tau_time[[k]][, tt]
      for (cc in seq_along(tk)) {
        nm_nat  <- c(nm_nat, sprintf("tau[%s,%s]_c%d", namvar[k], time_names[tt], cc))
        par_nat <- c(par_nat, tk[cc])
      }
    }
  }

  if (meas_invar0[["loadings"]] == 1L) {
    for (qq in seq_len(Q)) {
      for (k in seq_len(nvar)) {
        if (loadmask0[k, qq] == 0L) next
        nm_nat  <- c(nm_nat, sprintf("lambda[%s,%s]", namvar[k], factor_names[qq]))
        par_nat <- c(par_nat, lambda_hat[k, qq])
      }
    }
  } else {
    for (tt in seq_len(n)) for (qq in seq_len(Q)) {
      for (k in seq_len(nvar)) {
        if (loadmask0[k, qq] == 0L) next
        nm_nat  <- c(nm_nat, sprintf("lambda[%s,%s,%s]", namvar[k], factor_names[qq], time_names[tt]))
        par_nat <- c(par_nat, lambda_time[k, qq, tt])
      }
    }
  }

  if (meas_invar0[["intercepts"]] == 1L) {
    for (k in seq_len(nvar)) {
      nm_nat  <- c(nm_nat, sprintf("gamma[%s]", namvar[k]))
      par_nat <- c(par_nat, gamma_hat[k])
    }
    for (tt in if (n > 1L) 2:n else integer(0)) for (k in seq_len(nvar)) {
      if (gamma_time_free0[k, tt] != 1L) next
      nm_nat  <- c(nm_nat, sprintf("gamma[%s,%s]", namvar[k], time_names[tt]))
      par_nat <- c(par_nat, gamma_time[k, tt])
    }
  } else {
    for (tt in seq_len(n)) for (k in seq_len(nvar)) {
      nm_nat  <- c(nm_nat, sprintf("gamma[%s,%s]", namvar[k], time_names[tt]))
      par_nat <- c(par_nat, gamma_time[k, tt])
    }
  }

  if (meas_invar0[["residuals"]] == 1L) {
    for (k in seq_len(nvar)) {
      nm_nat  <- c(nm_nat, sprintf("theta[%s]", namvar[k]))
      par_nat <- c(par_nat, theta_hat[k, 1])
    }
  } else {
    for (k in seq_len(nvar)) for (tt in seq_len(n)) {
      nm_nat  <- c(nm_nat, sprintf("theta[%s,%s]", namvar[k], time_names[tt]))
      par_nat <- c(par_nat, theta_hat[k, tt])
    }
  }

  if (re_item0 == 1L) {
    for (k in seq_len(nvar)) {
      nm_nat  <- c(nm_nat, sprintf("Var_itemRE[%s]", namvar[k]))
      par_nat <- c(par_nat, sigmak_hat[k])
    }
  }

  for (tt in seq_len(n)) for (qq in seq_len(Q)) {
    nm_nat  <- c(nm_nat, sprintf("kappa[%s,%s]", time_names[tt], factor_names[qq]))
    par_nat <- c(par_nat, kappa_hat[tt, qq])
  }

  if (p_fix > 0L) {
    for (qq in seq_len(Q)) for (r in seq_len(p_fix)) {
      nm_nat <- c(nm_nat, sprintf("beta[%s,%s]", factor_names[qq], covariates[r]))
      par_nat <- c(par_nat, beta_hat[qq, r])
    }
  }

  for (d in seq_len(dim_eta)) {
    nm_nat  <- c(nm_nat, sprintf("Var_eta[%s]", eta_labels[d]))
    par_nat <- c(par_nat, Phi[d, d])
  }
  for (col in 1:dim_eta) if (col < dim_eta) for (row in (col + 1L):dim_eta) {
    nm_nat  <- c(nm_nat, sprintf("Cov_eta[%s,%s]", eta_labels[row], eta_labels[col]))
    par_nat <- c(par_nat, Phi[row, col])
  }



  names(par_nat) <- nm_nat

  se_nat_all <- rep(NA_real_, length(par_nat))
  names(se_nat_all) <- names(par_nat)

  tau_se <- tau_se_vector(
                          u_all = if (ntre0>0) best[1:ntre0] else numeric(0),
                          V_all = if (!is.null(vcov_opt) && ntre0>0) vcov_opt[1:ntre0,1:ntre0,drop=FALSE] else NULL,
                          cm1 = cm1, cm3 = cm3, namvar = namvar,
                          monotone_tau0 = monotone_tau0,
                          n = n, time_names = time_names,
                          tau_invar0 = meas_invar0[["thresholds"]],
                          threshold_id0 = threshold_id0
                        )

  se_nat_all[names(tau_se)] <- tau_se

  if (meas_invar0[["loadings"]] == 1L) {
    if (!is.null(vcov_opt)) {
      for (qq in seq_len(Q)) {
        ref_item <- marker_idx[qq]
        for (k in seq_len(nvar)) {
          if (loadmask0[k, qq] == 0L) next
          nmL <- sprintf("lambda[%s,%s]", namvar[k], factor_names[qq])
          if (k == ref_item) {
            se_nat_all[nmL] <- 0
          } else {
            nm_free <- sprintf("lambda_free[%s,%s]", namvar[k], factor_names[qq])
            j <- match(nm_free, names(par_opt))
            if (!is.na(j)) {
              vjj <- vcov_opt[j, j]
              if (is.finite(vjj) && vjj >= 0) se_nat_all[nmL] <- sqrt(vjj)
            }
          }
        }
      }
    } else {
      for (qq in seq_len(Q)) {
        ref_item <- marker_idx[qq]
        nmL <- sprintf("lambda[%s,%s]", namvar[ref_item], factor_names[qq])
        se_nat_all[nmL] <- 0
      }
    }
  } else {
    if (!is.null(vcov_opt)) {
      for (tt in seq_len(n)) for (qq in seq_len(Q)) {
        ref_item <- marker_idx[qq]
        for (k in seq_len(nvar)) {
          if (loadmask0[k, qq] == 0L) next
          nmL <- sprintf("lambda[%s,%s,%s]", namvar[k], factor_names[qq], time_names[tt])
          if (k == ref_item) {
            se_nat_all[nmL] <- 0
          } else {
            nm_free <- sprintf("lambda_free[%s,%s,%s]", namvar[k], factor_names[qq], time_names[tt])
            j <- match(nm_free, names(par_opt))
            if (!is.na(j)) {
              vjj <- vcov_opt[j, j]
              if (is.finite(vjj) && vjj >= 0) se_nat_all[nmL] <- sqrt(vjj)
            }
          }
        }
      }
    } else {
      for (tt in seq_len(n)) for (qq in seq_len(Q)) {
        ref_item <- marker_idx[qq]
        nmL <- sprintf("lambda[%s,%s,%s]", namvar[ref_item], factor_names[qq], time_names[tt])
        se_nat_all[nmL] <- 0
      }
    }
  }

  decode <- list(
    threshold_id = threshold_id,
    threshold_id_code = threshold_id0,
    mean_ident = mean_ident,
    mean_ident_code = mean_ident0,
    tau_list = tau_full,
    tau_time = tau_time,
    lambda   = lambda_hat,
    lambda_time = lambda_time,
    gamma    = gamma_hat,
    gamma_time = gamma_time,
    theta    = theta_hat,
    sigmak   = sigmak_hat,
    kappa    = kappa_hat,
    beta     = beta_hat,
    L_phi    = Lphi,
    Phi      = Phi
  )


  se_opt <- rep(NA_real_, npm0); names(se_opt) <- names(par_opt)
  if (!is.null(vcov_opt)) {
    d <- diag(vcov_opt)
    d[d < 0] <- NA_real_
    se_opt <- sqrt(d)
  }

  se_nat <- rep(NA_real_, length(par_interest)); names(se_nat) <- names(par_interest)

  if (!is.null(vcov_opt)) {

    se_exp <- function(name_nat, name_log) {
      j <- match(name_log, names(par_opt))
      if (is.na(j)) return(NA_real_)
      est <- unname(par_interest[name_nat])
      vjj <- vcov_opt[j, j]
      if (!is.finite(vjj) || vjj < 0) return(NA_real_)
      abs(est) * sqrt(vjj)
    }

    for (nm in names(par_interest)) {
      j <- match(nm, names(par_opt))
      if (!is.na(j)) {
        vjj <- vcov_opt[j, j]
        if (is.finite(vjj) && vjj >= 0) se_nat[nm] <- sqrt(vjj)
      }
    }

    for (nm in names(par_interest)) {
      if (grepl("^theta\\[", nm)) {
        nm_log <- sub("^theta\\[", "log_theta[", nm)
        se_nat[nm] <- se_exp(nm, nm_log)
      }
      if (grepl("^Var_itemRE\\[", nm)) {
        nm_log <- sub("^Var_itemRE\\[", "log_sigmak[", nm)
        se_nat[nm] <- se_exp(nm, nm_log)
      }
      if (grepl("^Var_eta1\\[", nm)) {
        nm_log <- sub("^Var_eta1\\[", "log_Sigma1[", nm)
        se_nat[nm] <- se_exp(nm, nm_log)
      }
    }

    nm_Lblock <- character(0)
    nm_Lblock <- c(nm_Lblock, sprintf("log_Lphi[%s]", eta_labels))
    for (col in 1:dim_eta) if (col < dim_eta) for (row in (col + 1L):dim_eta) {
      nm_Lblock <- c(nm_Lblock, sprintf("Lphi[%s,%s]", eta_labels[row], eta_labels[col]))
    }
    idx_Lblock <- match(nm_Lblock, names(par_opt))
    if (anyNA(idx_Lblock)) stop("Delta method: non riesco a ricostruire gli indici del blocco Lphi in par_opt.")

    V_L <- vcov_opt[idx_Lblock, idx_Lblock, drop = FALSE]

    pos_map <- matrix(NA_integer_, dim_eta, dim_eta)
    for (a in 1:dim_eta) pos_map[a, a] <- a
    pp <- dim_eta
    for (col in 1:dim_eta) if (col < dim_eta) for (row in (col + 1L):dim_eta) {
      pp <- pp + 1L
      pos_map[row, col] <- pp
    }

    grad_phi_rs <- function(r, s) {
      g <- numeric(length(nm_Lblock))
      m <- min(r, s)

      if (r == s) {
        for (k in 1:m) {
          dL <- 2 * Lphi[r, k]
          p <- pos_map[r, k]
          if (is.na(p)) next
          if (r == k) dL <- dL * Lphi[r, r]
          g[p] <- g[p] + dL
        }
      } else {
        for (k in 1:m) {

          d_rk <- Lphi[s, k]
          p_rk <- pos_map[r, k]
          if (!is.na(p_rk)) {
            if (r == k) d_rk <- d_rk * Lphi[r, r]
            g[p_rk] <- g[p_rk] + d_rk
          }

          d_sk <- Lphi[r, k]
          p_sk <- pos_map[s, k]
          if (!is.na(p_sk)) {
            if (s == k) d_sk <- d_sk * Lphi[s, s]
            g[p_sk] <- g[p_sk] + d_sk
          }
        }
      }
      g
    }

    for (d in 1:dim_eta) {
      nmv <- sprintf("Var_eta[%s]", eta_labels[d])
      g <- grad_phi_rs(d, d)
      vv <- as.numeric(t(g) %*% V_L %*% g)
      if (is.finite(vv) && vv >= 0) se_nat[nmv] <- sqrt(vv)
    }

    for (col in 1:dim_eta) if (col < dim_eta) for (row in (col + 1L):dim_eta) {
      nmc <- sprintf("Cov_eta[%s,%s]", eta_labels[row], eta_labels[col])
      g <- grad_phi_rs(row, col)
      vv <- as.numeric(t(g) %*% V_L %*% g)
      if (is.finite(vv) && vv >= 0) se_nat[nmc] <- sqrt(vv)
    }

  }

  zcrit <- qnorm(0.975)
  est_nat <- as.numeric(par_interest)
  z_nat <- est_nat / se_nat
  p_nat <- 2 * pnorm(-abs(z_nat))
  ci_low <- est_nat - zcrit * se_nat
  ci_high <- est_nat + zcrit * se_nat

  wald_nat <- data.frame(
    param   = names(par_interest),
    estimate= est_nat,
    se      = as.numeric(se_nat),
    z       = as.numeric(z_nat),
    p       = as.numeric(p_nat),
    reject_0.05 = as.logical(p_nat < 0.05),
    ci_low  = as.numeric(ci_low),
    ci_high = as.numeric(ci_high),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

common <- intersect(names(se_nat_all), names(se_nat))
se_nat_all[common] <- se_nat[common]

if (marker_mean_ident) {
  if (meas_invar0[["intercepts"]] == 1L) {
    for (qq in seq_len(Q)) {
      nm0 <- sprintf("gamma[%s]", namvar[marker_idx[qq]])
      if (nm0 %in% names(se_nat_all)) se_nat_all[nm0] <- 0
    }
  } else {
    for (tt in seq_len(n)) for (qq in seq_len(Q)) {
      nm0 <- sprintf("gamma[%s,%s]", namvar[marker_idx[qq]], time_names[tt])
      if (nm0 %in% names(se_nat_all)) se_nat_all[nm0] <- 0
    }
  }
} else {
  for (qq in factor_names) {
    nm0 <- sprintf("kappa[%s,%s]", time_names[1], qq)
    if (nm0 %in% names(se_nat_all)) se_nat_all[nm0] <- 0
  }
}

bin_items <- namvar[cm1 == 1L]
if (length(bin_items)) {
  if (!theta_time0) {
    for (it in bin_items) {
      nm0 <- sprintf("theta[%s]", it)
      if (nm0 %in% names(se_nat_all)) se_nat_all[nm0] <- 0
    }
  } else {
    for (it in bin_items) for (tt in time_names) {
      nm0 <- sprintf("theta[%s,%s]", it, tt)
      if (nm0 %in% names(se_nat_all)) se_nat_all[nm0] <- 0
    }
  }
}

  zcrit <- qnorm(0.975)

  est_all <- as.numeric(par_nat)
  se_all  <- as.numeric(se_nat_all[names(par_nat)])

  z_all <- ifelse(is.finite(se_all) & se_all > 0, est_all / se_all, NA_real_)
  p_all <- ifelse(is.finite(z_all), 2 * pnorm(-abs(z_all)), NA_real_)

  ci_low_all  <- ifelse(is.finite(se_all) & se_all > 0, est_all - zcrit * se_all, NA_real_)
  ci_high_all <- ifelse(is.finite(se_all) & se_all > 0, est_all + zcrit * se_all, NA_real_)

  wald_all <- data.frame(
    param   = names(par_nat),
    estimate= est_all,
    se      = se_all,
    z       = z_all,
    p       = p_all,
    reject_0.05 = as.logical(!is.na(p_all) & p_all < 0.05),
    ci_low  = ci_low_all,
    ci_high = ci_high_all,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  if (!valid_objective(opt$value)) stop("The final CFA likelihood could not be evaluated.")
  new_long_fit(list(
    conv    = opt$convergence,
    message = opt$message,
    loglik  = -opt$value,
    timing  = timing,
    parscale = ps0,
    par_opt = par_opt,
    par = par_nat,
    par_natural = par_nat,
    hessian = hess,
    hessian_raw = hessian_raw,
    hessian_method = hessian_method,
    hessian_step = hessian_step,
    score_at_optimum = score_at_optimum,
    hessian_eigenvalues = hessian_eigenvalues,
    hessian_positive_definite = hessian_positive_definite,
    se_opt  = se_opt,
    vcov_opt = vcov_opt,
    se  = se_nat_all,
    wald_nat = wald_nat,
    wald_all = wald_all,
    dims    = list(nsub = nsub, nvar = nvar, n = n, Q = Q, dim_eta = dim_eta, p_fix = p_fix),
        structure = list(
      loadmask = loadmask0,
      marker = setNames(namvar[marker_idx], paste0("F", seq_len(Q))),
      covariates = covariates,
      post_covariates = post_covariates,
      cov_time_mask = cov_time_mask,
      methInteg = methInteg,
      re_item = re_item0,
      two_tier = isTRUE(two_tier),
      two_tier_type = if (isTRUE(two_tier)) two_tier_type else NA_character_,
      analytic_score = isTRUE(use_analytic_score),
      threshold_id = threshold_id,
      threshold_id_code = threshold_id0,
      mean_ident = mean_ident,
      mean_ident_code = mean_ident0,
      threshold_fixed = if (threshold_id0 == 0L) {
        c(tau2 = 0, tau3 = 1)
      } else {
        c(tau1 = 0, tau2 = 1)
      },
      theta_time = theta_time0,
      meas_invar = meas_invar0,
      intercepts_free_items = intercepts_free_items,
      gamma_time_free = gamma_time_free0,
      phi_struct = phi_struct,
      nquad = if (!isTRUE(two_tier)) nquad0 else NULL,
      nquad_eta = if (isTRUE(two_tier)) nquad_eta0 else NULL,
      nquad_re = if (isTRUE(two_tier)) nquad_re0 else NULL,
      dimInt = if (!isTRUE(two_tier)) dimInt else NULL,
      dim_eta_int = if (isTRUE(two_tier)) dim_eta_int else NULL
    ),
    decode = decode,
        profile = profile_info
  ), "Longitudinal ordinal CFA", packed)
}
