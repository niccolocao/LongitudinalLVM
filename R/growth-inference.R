loglik_growth_tt_i_vec <- function(b0,
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
                                   delta_free0 = 1L,
                                   monotone_tau0 = 0L,
                                   Q0 = 1L,
                                   loadmask0 = NULL,
                             p_fix0 = NULL,
                             p_ran0 = NULL,
                             X0 = NULL,
                             Z0 = NULL,
                             delta_mask0 = NULL,
                             delta_x_col0 = NULL,
                             delta_z_col0 = NULL,
                             subject_specific_design0 = NULL,
                             prepared_design0 = FALSE) {
  if (is.null(loadmask0)) {
    loadmask0 <- matrix(1L, nrow = as.integer(nvar0), ncol = as.integer(Q0))
  }
  loadmask0 <- matrix(as.integer(loadmask0), nrow = as.integer(nvar0), ncol = as.integer(Q0))

  design0 <- .growth_prepare_lowlevel_call(
    nsub0, n0, Q0, p_fix0, p_ran0, X0, Z0,
    delta_mask0, delta_x_col0, delta_z_col0,
    subject_specific_design0,
    prepared_design0 = prepared_design0
  )

  if (re_item0 == 0L) {
    nquad_re0 <- 1L
    nodes_re0 <- 0
    weights_re0 <- 1
  }

  native_call("loglik_growth_tt_i_vec",
     b0 = as.double(b0),
     Y0 = as.integer(Y0),
     nvar0 = as.integer(nvar0),
     nsub0 = as.integer(nsub0),
     nobs0 = as.integer(nobs0),
     n0 = as.integer(n0),
     nmes0 = as.integer(nmes0),
     c0 = as.integer(c0),
     cm10 = as.integer(cm10),
     ntre0 = as.integer(ntre0),
     nquad_eta0 = as.integer(nquad_eta0),
     nodes_eta0 = as.double(nodes_eta0),
     weights_eta0 = as.double(weights_eta0),
     nquad_re0 = as.integer(nquad_re0),
     nodes_re0 = as.double(nodes_re0),
     weights_re0 = as.double(weights_re0),
     nitem0 = as.integer(nitem0),
     npm0 = as.integer(npm0),
     re_item0 = as.integer(re_item0),
     theta_time0 = as.integer(theta_time0),
     delta_free0 = as.integer(delta_free0),
     monotone_tau0 = as.integer(monotone_tau0),
     Q0 = as.integer(Q0),
     loadmask0 = as.integer(loadmask0),
     p_fix0 = design0$p_fix0,
     p_ran0 = design0$p_ran0,
     X0 = design0$X0,
     Z0 = design0$Z0,
     delta_mask0 = design0$delta_mask0,
     delta_x_col0 = design0$delta_x_col0,
     delta_z_col0 = design0$delta_z_col0,
     loglik_i = double(as.integer(nsub0)))$loglik_i
}

loglik_growth_tt_agh_i_vec <- function(b0,
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
                                       delta_free0 = 1L,
                                       monotone_tau0 = 0L,
                                       update_gh0 = 1L,
                                       fit_id0 = 1L,
                                       Q0 = 1L,
                                       loadmask0 = NULL,
                             p_fix0 = NULL,
                             p_ran0 = NULL,
                             X0 = NULL,
                             Z0 = NULL,
                             delta_mask0 = NULL,
                             delta_x_col0 = NULL,
                             delta_z_col0 = NULL,
                             subject_specific_design0 = NULL,
                             prepared_design0 = FALSE) {
  if (is.null(loadmask0)) {
    loadmask0 <- matrix(1L, nrow = as.integer(nvar0), ncol = as.integer(Q0))
  }
  loadmask0 <- matrix(as.integer(loadmask0), nrow = as.integer(nvar0), ncol = as.integer(Q0))

  design0 <- .growth_prepare_lowlevel_call(
    nsub0, n0, Q0, p_fix0, p_ran0, X0, Z0,
    delta_mask0, delta_x_col0, delta_z_col0,
    subject_specific_design0,
    prepared_design0 = prepared_design0
  )

  if (re_item0 == 0L) {
    nquad_re0 <- 1L
    nodes_re0 <- 0
    weights_re0 <- 1
  }

  native_call("loglik_growth_tt_agh_i_vec",
     b0 = as.double(b0),
     Y0 = as.integer(Y0),
     nvar0 = as.integer(nvar0),
     nsub0 = as.integer(nsub0),
     nobs0 = as.integer(nobs0),
     n0 = as.integer(n0),
     nmes0 = as.integer(nmes0),
     c0 = as.integer(c0),
     cm10 = as.integer(cm10),
     ntre0 = as.integer(ntre0),
     nquad_eta0 = as.integer(nquad_eta0),
     nodes_eta0 = as.double(nodes_eta0),
     weights_eta0 = as.double(weights_eta0),
     nquad_re0 = as.integer(nquad_re0),
     nodes_re0 = as.double(nodes_re0),
     weights_re0 = as.double(weights_re0),
     nitem0 = as.integer(nitem0),
     npm0 = as.integer(npm0),
     re_item0 = as.integer(re_item0),
     theta_time0 = as.integer(theta_time0),
     delta_free0 = as.integer(delta_free0),
     monotone_tau0 = as.integer(monotone_tau0),
     update_gh0 = as.integer(update_gh0),
     fit_id0 = as.integer(fit_id0),
     Q0 = as.integer(Q0),
     loadmask0 = as.integer(loadmask0),
     p_fix0 = design0$p_fix0,
     p_ran0 = design0$p_ran0,
     X0 = design0$X0,
     Z0 = design0$Z0,
     delta_mask0 = design0$delta_mask0,
     delta_x_col0 = design0$delta_x_col0,
     delta_z_col0 = design0$delta_z_col0,
     loglik_i = double(as.integer(nsub0)))$loglik_i
}

loglik_growth_tt_agh_full_equiv_i_vec <- function(b0,
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
                                                   delta_free0 = 1L,
                                                   monotone_tau0 = 0L,
                                                   update_gh0 = 1L,
                                                   fit_id0 = 1L,
                                                   Q0 = 1L,
                                                   loadmask0 = NULL,
                             p_fix0 = NULL,
                             p_ran0 = NULL,
                             X0 = NULL,
                             Z0 = NULL,
                             delta_mask0 = NULL,
                             delta_x_col0 = NULL,
                             delta_z_col0 = NULL,
                             subject_specific_design0 = NULL,
                             prepared_design0 = FALSE,
                             joint_rule0 = 0L,
                             joint_grid0 = NULL) {
  if (is.null(loadmask0)) {
    loadmask0 <- matrix(1L, nrow = as.integer(nvar0), ncol = as.integer(Q0))
  }
  loadmask0 <- matrix(as.integer(loadmask0), nrow = as.integer(nvar0), ncol = as.integer(Q0))

  design0 <- .growth_prepare_lowlevel_call(
    nsub0, n0, Q0, p_fix0, p_ran0, X0, Z0,
    delta_mask0, delta_x_col0, delta_z_col0,
    subject_specific_design0,
    prepared_design0 = prepared_design0
  )

  if (re_item0 == 0L) {
    nquad_re0 <- 1L
    nodes_re0 <- 0
    weights_re0 <- 1
  }

  if (joint_rule0 != 2L) {
    joint_grid0 <- list(modes = 0, roots = 0, logdet = 0)
  } else {
    dim_joint <- as.double(n0) * Q0 + nvar0 * re_item0
    sizes <- c(modes = dim_joint * nsub0, roots = dim_joint^2 * nsub0, logdet = nsub0)
    if (is.null(joint_grid0)) {
      joint_grid0 <- lapply(sizes, double)
    } else if (!all(lengths(joint_grid0[names(sizes)]) == sizes)) {
      stop("The saved joint AGH grid does not match the fitted dimensions.", call. = FALSE)
    }
  }
  out <- native_call("loglik_growth_tt_agh_full_equiv_i_vec",
     b0 = as.double(b0),
     Y0 = as.integer(Y0),
     nvar0 = as.integer(nvar0),
     nsub0 = as.integer(nsub0),
     nobs0 = as.integer(nobs0),
     n0 = as.integer(n0),
     nmes0 = as.integer(nmes0),
     c0 = as.integer(c0),
     cm10 = as.integer(cm10),
     ntre0 = as.integer(ntre0),
     nquad_eta0 = as.integer(nquad_eta0),
     nodes_eta0 = as.double(nodes_eta0),
     weights_eta0 = as.double(weights_eta0),
     nquad_re0 = as.integer(nquad_re0),
     nodes_re0 = as.double(nodes_re0),
     weights_re0 = as.double(weights_re0),
     nitem0 = as.integer(nitem0),
     npm0 = as.integer(npm0),
     re_item0 = as.integer(re_item0),
     theta_time0 = as.integer(theta_time0),
     delta_free0 = as.integer(delta_free0),
     monotone_tau0 = as.integer(monotone_tau0),
     update_gh0 = as.integer(update_gh0),
     fit_id0 = as.integer(fit_id0),
     Q0 = as.integer(Q0),
     loadmask0 = as.integer(loadmask0),
     p_fix0 = design0$p_fix0,
     p_ran0 = design0$p_ran0,
     X0 = design0$X0,
     Z0 = design0$Z0,
     delta_mask0 = design0$delta_mask0,
     delta_x_col0 = design0$delta_x_col0,
     delta_z_col0 = design0$delta_z_col0,
     joint_rule0 = as.integer(joint_rule0),
     joint_modes0 = as.double(joint_grid0$modes),
     joint_roots0 = as.double(joint_grid0$roots),
     joint_logdet0 = as.double(joint_grid0$logdet),
     loglik_i = double(as.integer(nsub0)))
  result <- out$loglik_i
  if (joint_rule0 == 2L && update_gh0 == 1L) {
    attr(result, "joint_grid") <- list(type = "AGH_JOINT", modes = out$joint_modes0,
                                      roots = out$joint_roots0, logdet = out$joint_logdet0)
  }
  result
}

build_growth_tt_c_inputs <- function(data,
                                     namvar,
                                     occ = "occ",
                                     ID = "ID",
                                     fit = NULL,
                                     loadmask = NULL,
                                     factor_names = NULL,
                                     methInteg = NULL,
                                     nquadgh = NULL,
                                     nquadqmc = 2000,
                                     re_item = NULL,
                                     theta_time = NULL,
                                     delta_free = NULL,
                                     monotone_tau = NULL,
                                     simple_growth = TRUE,
                                     same_design = TRUE,
                                     subject_specific_design = NULL,
                                     time_scores = NULL,
                                     fixed_covariates = NULL,
                                     time_interactions = NULL,
                                     X = NULL,
                                     Z = NULL,
                                     delta_mask = NULL,
                                     fixed_names = NULL,
                                     random_names = NULL) {

  if (!is.null(fit)) {
    if (!is.null(fit$c_inputs)) return(fit$c_inputs)
    stop(
      "fit$c_inputs mancante: non ricostruisco automaticamente X/Z, perch\u00e9 ",
      "potrei cambiare il modello stimato. Rifare il fit con la versione corrente ",
      "oppure chiamare build_growth_tt_c_inputs() senza fit e specificare il design."
    )
  }

  if (!all(c(ID, occ, namvar) %in% names(data))) {
    stop("variabili non trovate in data.")
  }

  if (is.null(methInteg)) methInteg <- "AGH"
  methInteg <- match.arg(
    methInteg,
    c("GH", "QMC", "AGH", "AGH_FULL_EQUIV", "AGH_RE_NESTED", "AGH_RE_FULL_EQUIV")
  )
  if (methInteg %in% c("AGH_RE_NESTED", "AGH_RE_FULL_EQUIV")) {
    stop("build_growth_tt_c_inputs non ricostruisce le griglie u/zeta; usa c_inputs salvato da growth_tt_fit().")
  }

  if (is.null(nquadgh)) nquadgh <- 5L
  nquadgh <- as.integer(nquadgh)
  if (is.null(re_item)) re_item <- TRUE
  if (is.null(theta_time)) theta_time <- FALSE
  if (is.null(delta_free)) delta_free <- isTRUE(simple_growth)
  if (is.null(monotone_tau)) monotone_tau <- FALSE
  if (!isTRUE(simple_growth) && isTRUE(delta_free)) {
    stop("delta_free \u00e8 disponibile soltanto con simple_growth=TRUE.")
  }

  nvar <- length(namvar)
  if (is.null(loadmask)) {
    loadmask <- matrix(1L, nrow = nvar, ncol = 1L)
  } else {
    loadmask <- as.matrix(loadmask)
    storage.mode(loadmask) <- "integer"
  }
  if (nrow(loadmask) != nvar) stop("nrow(loadmask) != length(namvar).")
  if (!all(loadmask %in% c(0L, 1L))) stop("loadmask deve contenere solo 0/1.")
  if (any(rowSums(loadmask) < 1L)) stop("ogni item deve caricare su almeno un fattore.")
  if (any(colSums(loadmask) < 1L)) stop("ogni fattore deve avere almeno un item.")

  Q <- as.integer(ncol(loadmask))
  if (is.null(factor_names)) factor_names <- paste0("F", seq_len(Q))
  if (length(factor_names) != Q) stop("length(factor_names) != Q.")

  data_sorted <- data[order(data[[ID]], data[[occ]]), , drop = FALSE]
  dat0 <- data_sorted[, c(ID, namvar, occ), drop = FALSE]
  occ_vals <- sort(unique(dat0[[occ]]))
  dat0$..t <- as.integer(match(dat0[[occ]], occ_vals))
  n <- as.integer(max(dat0$..t, na.rm = TRUE))
  if (!is.finite(n) || n < 2L) stop("servono almeno due occasioni.")

  IDs <- sort(unique(dat0[[ID]]))
  nsub <- length(IDs)

  design <- .growth_build_design(
    data_sorted = data_sorted,
    ID = ID,
    occ = occ,
    IDs = IDs,
    occ_vals = occ_vals,
    Q = Q,
    simple_growth = simple_growth,
    same_design = same_design,
    subject_specific_design = subject_specific_design,
    time_scores = time_scores,
    fixed_covariates = fixed_covariates,
    time_interactions = time_interactions,
    X = X,
    Z = Z,
    delta_mask = delta_mask,
    fixed_names = fixed_names,
    random_names = random_names
  )

  if (isTRUE(simple_growth) && any(!design$delta_mask %in% c(0, 1))) {
    stop("delta_mask deve contenere soltanto 0 e 1.")
  }
  if (isTRUE(simple_growth) && isTRUE(delta_free)) {
    active_delta <- vapply(
      seq_len(Q),
      function(q) sum(design$delta_mask[, , q] != 0),
      numeric(1L)
    )
    if (any(active_delta == 0L)) {
      stop("delta_free richiede almeno un loading selezionato in delta_mask per ogni fattore.")
    }
  }

  design0 <- .growth_lowlevel_design(
    nsub0 = nsub,
    n0 = n,
    Q0 = Q,
    p_fix0 = design$p_fix,
    p_ran0 = design$p_ran,
    X0 = design$X,
    Z0 = design$Z,
    delta_mask0 = design$delta_mask,
    delta_x_col0 = design$delta_x_col,
    delta_z_col0 = design$delta_z_col,
    subject_specific_design0 = design$subject_specific_design
  )

  c0 <- integer(nvar)
  cat_offset <- integer(nvar)
  for (k in seq_len(nvar)) {
    vals <- dat0[[namvar[k]]]
    vals <- vals[!is.na(vals)]
    if (!length(vals)) stop("item ", namvar[k], " ha solo NA.")
    mn <- min(vals)
    mx <- max(vals)
    if (mn == 0L) {
      c0[k] <- mx + 1L
      cat_offset[k] <- 1L
    } else if (mn == 1L) {
      c0[k] <- mx
      cat_offset[k] <- 0L
    } else {
      stop("item ", namvar[k], " ha categorie che non partono da 0 o 1.")
    }
    if (c0[k] < 2L) stop("item ", namvar[k], " ha meno di due categorie osservate.")
  }

  c0 <- as.integer(c0)
  cm10 <- as.integer(c0 - 1L)
  cm1 <- cm10
  cm3 <- pmax(cm1 - 2L, 0L)
  ntre0 <- as.integer(sum(cm3))

  nmes0 <- matrix(0L, nrow = nsub, ncol = nvar)
  Y_list <- vector("list", nsub * nvar)
  t_list <- vector("list", nsub * nvar)
  idx <- 0L

  for (ii in seq_len(nsub)) {
    di <- dat0[dat0[[ID]] == IDs[ii], c("..t", namvar), drop = FALSE]
    for (k in seq_len(nvar)) {
      idx <- idx + 1L
      y <- di[[namvar[k]]]
      tt <- di[["..t"]]
      ok <- which(!is.na(y))
      nmes0[ii, k] <- length(ok)
      if (length(ok)) {
        yv <- y[ok]
        Y_list[[idx]] <- as.integer(if (cat_offset[k] == 1L) yv + 1L else yv)
        t_list[[idx]] <- as.integer(tt[ok])
      }
    }
  }

  Y0 <- as.integer(unlist(Y_list))
  nitem0 <- as.integer(unlist(t_list))
  nobs <- length(Y0)

  re_item0 <- if (isTRUE(re_item)) 1L else 0L
  theta_time0 <- if (isTRUE(theta_time)) 1L else 0L
  delta_free0 <- if (isTRUE(delta_free)) 1L else 0L
  monotone_tau0 <- if (isTRUE(monotone_tau)) 1L else 0L

  eta_rule <- if (methInteg %in% c("GH", "AGH", "AGH_FULL_EQUIV")) {
    make_gh_tensor(dim_eta, nquadgh)
  } else {
    make_qmc(dim_eta, nquadqmc)
  }
  re_rule <- if (re_item0 == 1L) {
    make_gh_1d(nquadgh)
  } else {
    list(nquad = 1L, nodes = 0, weights = 1)
  }

  n_theta_est <- if (theta_time0 == 0L) sum(cm1 > 1L) else sum(cm1 > 1L) * n
  n_lambda_free <- as.integer(sum(colSums(loadmask) - 1L))
  n_nu_free <- as.integer(nvar - Q)
  n_growth_est <- as.integer(Q * p_fix + if (delta_free0 == 1L) Q else 0L)
  n_Lu_est <- as.integer(dim_u * (dim_u + 1L) / 2L)
  n_zeta_est <- as.integer(n * Q)
  n_item_re_est <- if (re_item0 == 1L) nvar else 0L

  npm0 <- as.integer(
    ntre0 + n_lambda_free + n_nu_free + n_theta_est +
      n_growth_est + n_Lu_est + n_zeta_est + n_item_re_est
  )

  list(
    Y0_c = as.integer(Y0),
    nvar0_c = as.integer(nvar),
    nsub0_c = as.integer(nsub),
    nobs0_c = as.integer(nobs),
    n0_c = as.integer(n),
    nmes0_c = as.integer(nmes0),
    c0_c = as.integer(c0),
    cm10_c = as.integer(cm10),
    ntre0_c = as.integer(ntre0),
    nquad_eta0_c = as.integer(eta_rule$nquad),
    nodes_eta0_c = as.double(eta_rule$nodes),
    weights_eta0_c = as.double(eta_rule$weights),
    nquad_re0_c = as.integer(re_rule$nquad),
    nodes_re0_c = as.double(re_rule$nodes),
    weights_re0_c = as.double(re_rule$weights),
    nitem0_c = as.integer(nitem0),
    npm0_c = as.integer(npm0),
    re_item0_c = as.integer(re_item0),
    theta_time0_c = as.integer(theta_time0),
    delta_free0_c = as.integer(delta_free0),
    monotone_tau0_c = as.integer(monotone_tau0),
    Q0_c = as.integer(Q),
    loadmask0_c = as.integer(loadmask),
    p_fix0_c = design0$p_fix0,
    p_ran0_c = design0$p_ran0,
    subject_specific_design0_c = design0$subject_specific_design0,
    X0_c = design0$X0,
    Z0_c = design0$Z0,
    delta_mask0_c = design0$delta_mask0,
    delta_x_col0_c = design0$delta_x_col0,
    delta_z_col0_c = design0$delta_z_col0,
    IDs = IDs,
    namvar = namvar,
    factor_names = factor_names,
    methInteg = methInteg,
    dims = list(
      nsub = nsub,
      nvar = nvar,
      n = n,
      Q = Q,
      dim_eta = dim_eta,
      dim_u = dim_u,
      p_fix = p_fix,
      p_ran = p_ran,
      npm0 = npm0
    )
  )
}

#evaluate subject likelihoods and derivatives
growth_tt_individual_derivatives <- function(fit,
                                             data = NULL,
                                             namvar = fit$c_inputs$namvar,
                                             occ = "occ",
                                             ID = "ID",
                                             nquadgh = NULL,
                                             nquadqmc = 2000,
                                             eps_score = 1e-5,
                                             eps_hess = 1e-4,
                                             par_index = NULL,
                                             compute_score = TRUE,
                                             compute_hessian = TRUE,
                                             agh_fixed_grid = TRUE,
                                             progress = FALSE,
                                             attach_to_fit = FALSE) {
  if (is.null(fit$par_opt)) stop("fit$par_opt mancante.")
  if (!is.null(fit$structure$methInteg) &&
      fit$structure$methInteg %in% c("AGH_RE_NESTED", "AGH_RE_FULL_EQUIV")) {
    stop("Le derivate individuali per AGH_RE_NESTED e AGH_RE_FULL_EQUIV richiedono le corrispondenti routine Fortran i_vec, non ancora incluse.")
  }
  par_hat <- as.double(fit$par_opt)
  par_names <- names(fit$par_opt)
  if (is.null(par_names)) par_names <- paste0("par", seq_along(par_hat))
  names(par_hat) <- par_names
  npm0 <- length(par_hat)

  if (is.null(par_index)) par_index <- seq_len(npm0)
  par_index <- as.integer(par_index)
  if (any(par_index < 1L | par_index > npm0)) stop("par_index fuori range.")
  Psub <- length(par_index)

  if (!is.null(fit$c_inputs)) {
    ci <- fit$c_inputs
    if (!is.null(data)) {
      packed <- prepare_long_data(data, namvar, occ, ID)
      if (!identical(namvar, ci$namvar) || !identical(packed$subjects, ci$IDs) ||
          !identical(as.integer(packed$y), ci$Y0_c) ||
          !identical(as.integer(packed$nmes), ci$nmes0_c) ||
          !identical(as.integer(packed$time), ci$nitem0_c)) {
        stop("The supplied data differ from the fitted subject and item responses.")
      }
    }
    methInteg <- if (!is.null(fit$structure$methInteg)) fit$structure$methInteg else ci$methInteg
    fit_id0_c <- as.integer(ci$fit_id0_c)
    if (methInteg %in% c("AGH", "AGH_FULL_EQUIV")) {
      if (is.null(fit$agh_grid)) {
        stop("fit$agh_grid mancante: la griglia AGH alla convergenza deve essere salvata nel fit.")
      }
      if (!identical(fit$agh_grid$type, "AGH_JOINT")) growth_tt_set_agh_grid(fit$agh_grid)
    }
  } else {
    if (is.null(data) || is.null(namvar)) {
      stop("Supply the original data and item names for a fit without stored inputs.")
    }
    ci <- build_growth_tt_c_inputs(
      data = data,
      namvar = namvar,
      occ = occ,
      ID = ID,
      fit = fit,
      nquadgh = nquadgh,
      nquadqmc = nquadqmc
    )
    methInteg <- ci$methInteg
    fit_id0_c <- next_fit_id()
  }

  joint_rule <- if (is.null(ci$joint_rule_c)) 0L else ci$joint_rule_c
  if (ci$npm0_c != npm0) {
    stop("npm0 != length(fit$par_opt): ", ci$npm0_c, " vs ", npm0)
  }
  subject_specific_design0_ci <- if (is.null(ci$subject_specific_design0_c)) {
    NULL
  } else {
    as.logical(ci$subject_specific_design0_c)
  }

  check_symbol <- function(symb) {
    if (!native_has(symb)) stop("simbolo ", symb, " non caricato. ricompila fort.so con la patch individuale.")
  }

  ll_i_vec <- function(par, update_gh_eval = 0L) {
    par <- as.double(par)
    if (length(par) != npm0) stop("length(par) != npm0.")
    X0 <- ci$X0_c
    Z0 <- ci$Z0_c
    if (!is.null(ci$rho_spec)) {
      native <- growth_rho_native(par, ci$rho_spec, X0, Z0)
      par <- native$par
      X0 <- native$X
      Z0 <- native$Z
    }
    npm_native <- length(par)

    if (!is.null(ci$growth_options)) {
      go <- ci$growth_options
      .growth_set_fortran_options(
        identification = go$identification,
        threshold_identification = go$threshold_identification,
        growth_shape = go$growth_shape,
        intercept_mask = go$intercept_mask,
        gamma_dep = go$gamma_dep,
        x_slope_col = go$x_slope_col,
        nvar = ci$dims$nvar,
        n = ci$dims$n,
        Q = ci$dims$Q,
        p_fix = ci$dims$p_fix,
        rho_fixed = go$rho_fixed,
        marker_index = go$marker_index,
        zero_time = go$zero_time
      )
    }

    if (methInteg == "AGH" && joint_rule == 0L) {
      check_symbol("loglik_growth_tt_agh_i_vec")
      out <- loglik_growth_tt_agh_i_vec(
        par,
        ci$Y0_c, ci$nvar0_c, ci$nsub0_c, ci$nobs0_c, ci$n0_c,
        ci$nmes0_c, ci$c0_c, ci$cm10_c, ci$ntre0_c,
        ci$nquad_eta0_c, ci$nodes_eta0_c, ci$weights_eta0_c,
        ci$nquad_re0_c, ci$nodes_re0_c, ci$weights_re0_c,
        ci$nitem0_c, npm_native, ci$re_item0_c, ci$theta_time0_c,
        ci$delta_free0_c, ci$monotone_tau0_c,
        as.integer(update_gh_eval), fit_id0_c,
        ci$Q0_c, ci$loadmask0_c,
        ci$p_fix0_c, ci$p_ran0_c, X0, Z0, ci$delta_mask0_c,
        ci$delta_x_col0_c, ci$delta_z_col0_c,
        subject_specific_design0_ci, TRUE
      )
    } else if (methInteg == "AGH_FULL_EQUIV" || joint_rule != 0L) {
      check_symbol("loglik_growth_tt_agh_full_equiv_i_vec")
      out <- loglik_growth_tt_agh_full_equiv_i_vec(
        par,
        ci$Y0_c, ci$nvar0_c, ci$nsub0_c, ci$nobs0_c, ci$n0_c,
        ci$nmes0_c, ci$c0_c, ci$cm10_c, ci$ntre0_c,
        ci$nquad_eta0_c, ci$nodes_eta0_c, ci$weights_eta0_c,
        ci$nquad_re0_c, ci$nodes_re0_c, ci$weights_re0_c,
        ci$nitem0_c, npm_native, ci$re_item0_c, ci$theta_time0_c,
        ci$delta_free0_c, ci$monotone_tau0_c,
        as.integer(update_gh_eval), fit_id0_c,
        ci$Q0_c, ci$loadmask0_c,
        ci$p_fix0_c, ci$p_ran0_c, X0, Z0, ci$delta_mask0_c,
        ci$delta_x_col0_c, ci$delta_z_col0_c,
        subject_specific_design0_ci, TRUE,
        joint_rule0 = joint_rule, joint_grid0 = fit$agh_grid
      )
    } else {
      check_symbol("loglik_growth_tt_i_vec")
      out <- loglik_growth_tt_i_vec(
        par,
        ci$Y0_c, ci$nvar0_c, ci$nsub0_c, ci$nobs0_c, ci$n0_c,
        ci$nmes0_c, ci$c0_c, ci$cm10_c, ci$ntre0_c,
        ci$nquad_eta0_c, ci$nodes_eta0_c, ci$weights_eta0_c,
        ci$nquad_re0_c, ci$nodes_re0_c, ci$weights_re0_c,
        ci$nitem0_c, npm_native, ci$re_item0_c, ci$theta_time0_c,
        ci$delta_free0_c, ci$monotone_tau0_c,
        ci$Q0_c, ci$loadmask0_c,
        ci$p_fix0_c, ci$p_ran0_c, X0, Z0, ci$delta_mask0_c,
        ci$delta_x_col0_c, ci$delta_z_col0_c,
        subject_specific_design0_ci, TRUE
      )
    }

    if (any(!is.finite(out) | out <= -1e8)) {
      stop("A subject likelihood could not be evaluated.")
    }
    names(out) <- as.character(ci$IDs)
    out
  }

  if (methInteg %in% c("AGH", "AGH_FULL_EQUIV")) {
    if (!isTRUE(agh_fixed_grid)) {
      stop("Per AGH score/Hessiana usa agh_fixed_grid=TRUE: la griglia deve restare quella salvata alla convergenza.")
    }

    if (isTRUE(progress)) {
      cat("restore AGH grid saved at convergence; update_gh0 = 0 for derivatives\n")
    }

    if (joint_rule == 0L) growth_tt_set_agh_grid(fit$agh_grid)
    loglik_i_refresh <- NULL
  }

  deriv_update <- if (methInteg %in% c("AGH", "AGH_FULL_EQUIV") &&
                         isTRUE(agh_fixed_grid)) 0L else 1L
  loglik_i <- ll_i_vec(par_hat, update_gh_eval = deriv_update)
  nsub <- length(loglik_i)

  h_score <- eps_score * (abs(par_hat) + 1)
  h_hess <- eps_hess * (abs(par_hat) + 1)

  score_ind <- NULL
  hess_ind <- NULL

  plus_score <- vector("list", npm0)
  minus_score <- vector("list", npm0)

  if (isTRUE(compute_score)) {
    score_ind <- matrix(NA_real_, nrow = nsub, ncol = Psub,
                        dimnames = list(as.character(ci$IDs), par_names[par_index]))

    for (aa in seq_along(par_index)) {
      p <- par_index[aa]
      if (isTRUE(progress)) cat("score parameter", aa, "/", Psub, ":", par_names[p], "\n")
      pp <- par_hat
      pm <- par_hat
      pp[p] <- pp[p] + h_score[p]
      pm[p] <- pm[p] - h_score[p]
      lp <- ll_i_vec(pp, update_gh_eval = deriv_update)
      lm <- ll_i_vec(pm, update_gh_eval = deriv_update)
      plus_score[[p]] <- lp
      minus_score[[p]] <- lm
      score_ind[, aa] <- (lp - lm) / (2 * h_score[p])
    }
  }

  if (isTRUE(compute_hessian)) {
    hess_ind <- array(NA_real_, dim = c(nsub, Psub, Psub),
                      dimnames = list(as.character(ci$IDs), par_names[par_index], par_names[par_index]))

    plus_hess <- vector("list", npm0)
    minus_hess <- vector("list", npm0)

    for (aa in seq_along(par_index)) {
      p <- par_index[aa]
      if (isTRUE(progress)) cat("hessian diagonal", aa, "/", Psub, ":", par_names[p], "\n")
      pp <- par_hat
      pm <- par_hat
      pp[p] <- pp[p] + h_hess[p]
      pm[p] <- pm[p] - h_hess[p]
      lp <- ll_i_vec(pp, update_gh_eval = deriv_update)
      lm <- ll_i_vec(pm, update_gh_eval = deriv_update)
      plus_hess[[p]] <- lp
      minus_hess[[p]] <- lm
      hess_ind[, aa, aa] <- (lp - 2 * loglik_i + lm) / (h_hess[p]^2)
    }

    if (Psub >= 2L) {
      for (aa in seq_len(Psub - 1L)) {
        p <- par_index[aa]
        for (bb in (aa + 1L):Psub) {
          q <- par_index[bb]
          if (isTRUE(progress)) cat("hessian offdiag", aa, bb, "/", Psub, ":", par_names[p], "x", par_names[q], "\n")
          p1 <- par_hat
          p2 <- par_hat
          p3 <- par_hat
          p4 <- par_hat
          p1[p] <- p1[p] + h_hess[p]; p1[q] <- p1[q] + h_hess[q]
          p2[p] <- p2[p] + h_hess[p]; p2[q] <- p2[q] - h_hess[q]
          p3[p] <- p3[p] - h_hess[p]; p3[q] <- p3[q] + h_hess[q]
          p4[p] <- p4[p] - h_hess[p]; p4[q] <- p4[q] - h_hess[q]
          fpp <- ll_i_vec(p1, update_gh_eval = deriv_update)
          fpm <- ll_i_vec(p2, update_gh_eval = deriv_update)
          fmp <- ll_i_vec(p3, update_gh_eval = deriv_update)
          fmm <- ll_i_vec(p4, update_gh_eval = deriv_update)
          val <- (fpp - fpm - fmp + fmm) / (4 * h_hess[p] * h_hess[q])
          hess_ind[, aa, bb] <- val
          hess_ind[, bb, aa] <- val
        }
      }
    }
  }

  loglik_sum <- sum(loglik_i)
  fit_loglik <- if (!is.null(fit$loglik)) fit$loglik else NA_real_

  score_total <- if (!is.null(score_ind)) colSums(score_ind) else NULL
  hess_total_loglik <- if (!is.null(hess_ind)) apply(hess_ind, c(2, 3), sum) else NULL
  hess_total_negloglik <- if (!is.null(hess_total_loglik)) -hess_total_loglik else NULL

  hess_symmetry_max <- NA_real_
  if (!is.null(hess_ind)) {
    hess_symmetry_max <- max(abs(hess_ind - aperm(hess_ind, c(1, 3, 2))), na.rm = TRUE)
  }

  hess_compare <- NULL
  if (!is.null(fit$hessian) && !is.null(hess_total_negloglik) &&
      all(dim(fit$hessian) == dim(hess_total_negloglik)) && Psub == npm0) {
    D <- hess_total_negloglik - fit$hessian
    hess_compare <- list(
      max_abs_diff = max(abs(D), na.rm = TRUE),
      rel_fro_diff = sqrt(sum(D^2, na.rm = TRUE)) / max(1, sqrt(sum(fit$hessian^2, na.rm = TRUE)))
    )
  }

  eig_neg_min <- NA_real_
  if (!is.null(hess_total_negloglik) && nrow(hess_total_negloglik) == ncol(hess_total_negloglik)) {
    eig_neg_min <- tryCatch(min(eigen((hess_total_negloglik + t(hess_total_negloglik)) / 2,
                                      symmetric = TRUE, only.values = TRUE)$values),
                            error = function(e) NA_real_)
  }

  sandwich <- NULL
  if (!is.null(score_ind) && !is.null(hess_total_negloglik) && Psub == npm0) {
    A <- (hess_total_negloglik + t(hess_total_negloglik)) / 2
    B <- crossprod(score_ind)
    Ainv <- tryCatch(solve(A), error = function(e) NULL)
    if (!is.null(Ainv)) {
      V <- Ainv %*% B %*% Ainv
      V <- (V + t(V)) / 2
      sandwich <- list(
        bread = A,
        meat = B,
        vcov = V,
        se = sqrt(pmax(diag(V), 0))
      )
      names(sandwich$se) <- colnames(score_ind)
    }
  }

  checks <- list(
    loglik_i_sum = loglik_sum,
    fit_loglik = fit_loglik,
    loglik_diff = if (is.finite(fit_loglik)) loglik_sum - fit_loglik else NA_real_,
    max_abs_score_total = if (!is.null(score_total)) max(abs(score_total), na.rm = TRUE) else NA_real_,
    score_total = score_total,
    hess_symmetry_max = hess_symmetry_max,
    hess_neg_min_eigen = eig_neg_min,
    hess_compare_to_fit = hess_compare,
    agh_derivative_grid = if (methInteg == "AGH") if (isTRUE(agh_fixed_grid)) "fixed_grid_after_refresh" else "full_refresh_each_eval" else NA_character_,
    finite_loglik_i = all(is.finite(loglik_i)),
    finite_score_ind = if (!is.null(score_ind)) all(is.finite(score_ind)) else NA,
    finite_hess_ind = if (!is.null(hess_ind)) all(is.finite(hess_ind)) else NA
  )

  out <- list(
    loglik_i = loglik_i,
    score_ind = score_ind,
    hess_ind = hess_ind,
    hess_total_loglik = hess_total_loglik,
    hess_total_negloglik = hess_total_negloglik,
    sandwich = sandwich,
    checks = checks,
    settings = list(
      eps_score = eps_score,
      eps_hess = eps_hess,
      par_index = par_index,
      par_names = par_names[par_index],
      agh_fixed_grid = agh_fixed_grid,
      deriv_update_gh0 = deriv_update,
      fit_id0 = fit_id0_c,
      dims = ci$dims
    ),
    ll_i_vec = ll_i_vec
  )

  if (isTRUE(attach_to_fit)) {
    fit$loglik_i <- out$loglik_i
    fit$score_ind <- out$score_ind
    fit$hess_ind <- out$hess_ind
    fit$hess_total_loglik_from_ind <- out$hess_total_loglik
    fit$hess_total_negloglik_from_ind <- out$hess_total_negloglik
    fit$sandwich <- out$sandwich
    fit$individual_derivative_checks <- out$checks
    fit$individual_derivative_settings <- out$settings
    return(fit)
  }

  out
}

growth_tt_get_agh_grid <- function() {
  info <- native_call(
    "get_agh_tt_grid_info",
    ready0 = integer(1),
    dim0 = integer(1),
    nquad0 = integer(1),
    nsub0 = integer(1),
    fit_id0 = integer(1)
  )

  if (info$ready0 != 1L) {
    stop("Griglia AGH non pronta in Fortran.")
  }

  dim0 <- info$dim0
  nquad0 <- info$nquad0
  nsub0 <- info$nsub0

  out <- native_call(
    "get_agh_tt_grid_copy",
    z_out = double(dim0 * nquad0 * nsub0),
    modes_out = double(dim0 * nsub0),
    logconst_out = double(nsub0),
    ready0 = integer(1),
    dim0 = integer(1),
    nquad0 = integer(1),
    nsub0 = integer(1),
    fit_id0 = integer(1)
  )

  if (out$ready0 != 1L) {
    stop("Copia griglia AGH fallita.")
  }

  list(
    ready = TRUE,
    type = "AGH",
    dim = out$dim0,
    nquad = out$nquad0,
    nsub = out$nsub0,
    fit_id0_c = out$fit_id0,
    z = array(out$z_out, dim = c(out$dim0, out$nquad0, out$nsub0)),
    modes = array(out$modes_out, dim = c(out$dim0, out$nsub0)),
    logconst = out$logconst_out
  )
}

growth_tt_get_agh_full_equiv_grid <- function() {
  info <- native_call(
    "get_growth_agh_tt_fe_grid_info",
    ready0 = integer(1),
    dim_eta0 = integer(1),
    dim_z0 = integer(1),
    nvar0 = integer(1),
    nquad_eta0 = integer(1),
    nquad_re0 = integer(1),
    nsub0 = integer(1),
    fit_id0 = integer(1)
  )

  if (info$ready0 != 1L) {
    stop("Griglia AGH full-equivalent non pronta in Fortran.")
  }

  out <- native_call(
    "get_growth_agh_tt_fe_grid_copy",
    zeta_eta_out = double(info$dim_eta0 * info$nquad_eta0 * info$nsub0),
    alpha_out = double(info$nvar0 * info$nquad_eta0 * info$nsub0),
    sdev_out = double(info$nvar0 * info$nsub0),
    logconst_out = double(info$nsub0),
    modes_out = double(info$dim_z0 * info$nsub0),
    ready0 = integer(1),
    dim_eta0 = integer(1),
    dim_z0 = integer(1),
    nvar0 = integer(1),
    nquad_eta0 = integer(1),
    nquad_re0 = integer(1),
    nsub0 = integer(1),
    fit_id0 = integer(1)
  )

  if (out$ready0 != 1L) {
    stop("Copia griglia AGH full-equivalent fallita.")
  }

  list(
    ready = TRUE,
    type = "AGH_FULL_EQUIV",
    dim_eta = out$dim_eta0,
    dim_z = out$dim_z0,
    nvar = out$nvar0,
    nquad_eta = out$nquad_eta0,
    nquad_re = out$nquad_re0,
    nsub = out$nsub0,
    fit_id0_c = out$fit_id0,
    zeta_eta = array(
      out$zeta_eta_out,
      dim = c(out$dim_eta0, out$nquad_eta0, out$nsub0)
    ),
    alpha = array(
      out$alpha_out,
      dim = c(out$nvar0, out$nquad_eta0, out$nsub0)
    ),
    sdev = array(out$sdev_out, dim = c(out$nvar0, out$nsub0)),
    logconst = out$logconst_out,
    modes = array(out$modes_out, dim = c(out$dim_z0, out$nsub0))
  )
}

growth_tt_set_agh_full_equiv_grid <- function(grid) {
  if (is.null(grid) || !isTRUE(grid$ready)) {
    stop("grid AGH full-equivalent mancante o non pronta.")
  }

  out <- native_call(
    "set_growth_agh_tt_fe_grid_copy",
    zeta_eta_in = as.double(grid$zeta_eta),
    alpha_in = as.double(grid$alpha),
    sdev_in = as.double(grid$sdev),
    logconst_in = as.double(grid$logconst),
    modes_in = as.double(grid$modes),
    dim_eta0 = as.integer(grid$dim_eta),
    dim_z0 = as.integer(grid$dim_z),
    nvar0 = as.integer(grid$nvar),
    nquad_eta0 = as.integer(grid$nquad_eta),
    nquad_re0 = as.integer(grid$nquad_re),
    nsub0 = as.integer(grid$nsub),
    fit_id0 = as.integer(grid$fit_id0_c),
    ready0 = integer(1)
  )

  if (out$ready0 != 1L) {
    stop("Ripristino griglia AGH full-equivalent fallito.")
  }

  invisible(TRUE)
}

growth_tt_set_agh_grid <- function(grid) {
  if (is.null(grid) || !isTRUE(grid$ready)) {
    stop("grid AGH mancante o non pronta.")
  }

  if (identical(grid$type, "AGH_FULL_EQUIV")) {
    return(growth_tt_set_agh_full_equiv_grid(grid))
  }

  out <- native_call(
    "set_agh_tt_grid_copy",
    z_in = as.double(grid$z),
    modes_in = as.double(grid$modes),
    logconst_in = as.double(grid$logconst),
    dim0 = as.integer(grid$dim),
    nquad0 = as.integer(grid$nquad),
    nsub0 = as.integer(grid$nsub),
    fit_id0 = as.integer(grid$fit_id0_c),
    ready0 = integer(1)
  )

  if (out$ready0 != 1L) {
    stop("Ripristino griglia AGH fallito.")
  }

  invisible(TRUE)
}
add_agh_re_wald <- function(
  fit,
  hess_eps = 1e-3,
  jac_eps = 1e-5,
  progress_every = 100L,
  checkpoint_file = "agh_re_full_hessian_checkpoint.rds"
) {

  cat("[1/6] controllo configurazione\n")

  ci <- fit$c_inputs
  st <- fit$structure
  go <- ci$growth_options
  dims <- if (!is.null(ci$dims)) ci$dims else fit$dims

  if (is.null(ci)) stop("fit$c_inputs mancante.")

  meth <- if (!is.null(ci$methInteg)) ci$methInteg else st$methInteg

  if (!identical(meth, "AGH_RE_FULL_EQUIV")) {
    stop("questa funzione e' specifica per AGH_RE_FULL_EQUIV.")
  }

  if (!identical(go$identification, "baseline_std")) {
    stop("questa funzione e' specifica per identification = baseline_std.")
  }

  if (!identical(go$growth_shape, "rho")) {
    stop("questa funzione e' specifica per growth_shape = rho.")
  }

  nvar <- as.integer(dims$nvar)
  nsub <- as.integer(dims$nsub)
  n <- as.integer(dims$n)
  Q <- as.integer(dims$Q)
  p_fix <- as.integer(dims$p_fix)
  p_ran <- as.integer(dims$p_ran)

  namvar <- ci$namvar
  if (is.null(namvar)) namvar <- rownames(fit$decode$lambda)

  factor_names <- st$growth$factor_names
  if (is.null(factor_names)) factor_names <- colnames(fit$decode$lambda)

  fixed_names <- st$growth$fixed_names
  if (!is.list(fixed_names)) {
    fixed_names <- rep(list(fixed_names), Q)
  }

  u_names <- rownames(fit$decode$Lu)

  if (is.null(u_names)) {
    stop("nomi dei random effects non trovati.")
  }

  cm1 <- as.integer(ci$cm10_c)

  loadmask <- matrix(
    as.integer(ci$loadmask0_c),
    nrow = nvar,
    ncol = Q
  )

  intercept_mask <- matrix(
    as.logical(go$intercept_mask),
    nrow = nvar,
    ncol = n
  )

  theta_time <- as.integer(ci$theta_time0_c) == 1L
  delta_free <- as.integer(ci$delta_free0_c) == 1L
  monotone_tau <- as.integer(ci$monotone_tau0_c) == 1L
  re_item <- as.integer(ci$re_item0_c) == 1L

  par_hat <- as.double(fit$par_opt)
  names(par_hat) <- names(fit$par_opt)

  npm <- length(par_hat)

  if (npm != ci$npm0_c) {
    stop("lunghezza par_opt diversa da npm0.")
  }

  fit_id_hess <- as.integer(
    sample.int(.Machine$integer.max, 1L)
  )

  set_options <- function() {
    .growth_set_fortran_options(
      identification = go$identification,
      threshold_identification = go$threshold_identification,
      growth_shape = go$growth_shape,
      intercept_mask = go$intercept_mask,
      gamma_dep = go$gamma_dep,
      x_slope_col = go$x_slope_col,
      nvar = nvar,
      n = n,
      Q = Q,
      p_fix = p_fix,
      rho_fixed = go$rho_fixed,
      marker_index = go$marker_index
    )
  }

  nll_raw <- function(par) {

    set_options()

    val <- loglik_growth_re_tt_agh_full_equiv(
      b0 = par,
      Y0 = ci$Y0_c,
      nvar0 = ci$nvar0_c,
      nsub0 = ci$nsub0_c,
      nobs0 = ci$nobs0_c,
      n0 = ci$n0_c,
      nmes0 = ci$nmes0_c,
      c0 = ci$c0_c,
      cm10 = ci$cm10_c,
      ntre0 = ci$ntre0_c,
      nquad_u0 = ci$nquad_u0_c,
      nodes_u0 = ci$nodes_u0_c,
      weights_u0 = ci$weights_u0_c,
      nquad_zeta0 = ci$nquad_zeta0_c,
      nodes_zeta0 = ci$nodes_zeta0_c,
      weights_zeta0 = ci$weights_zeta0_c,
      nquad_re0 = ci$nquad_re0_c,
      nodes_re0 = ci$nodes_re0_c,
      weights_re0 = ci$weights_re0_c,
      nitem0 = ci$nitem0_c,
      npm0 = ci$npm0_c,
      re_item0 = ci$re_item0_c,
      theta_time0 = ci$theta_time0_c,
      delta_free0 = ci$delta_free0_c,
      monotone_tau0 = ci$monotone_tau0_c,
      update_gh0 = 1L,
      fit_id0 = fit_id_hess,
      Q0 = ci$Q0_c,
      loadmask0 = ci$loadmask0_c,
      p_fix0 = ci$p_fix0_c,
      p_ran0 = ci$p_ran0_c,
      X0 = ci$X0_c,
      Z0 = ci$Z0_c,
      delta_mask0 = ci$delta_mask0_c,
      delta_x_col0 = ci$delta_x_col0_c,
      delta_z_col0 = ci$delta_z_col0_c,
      subject_specific_design0 = ci$subject_specific_design0_c,
      prepared_design0 = TRUE
    )

    val <- as.double(val)

    if (!is.finite(val) || abs(val) >= 1e8) {
      stop("likelihood non valida durante il calcolo della Hessiana.")
    }

    val
  }

  cat("[2/6] preparo Hessiana diretta\n")

  h <- hess_eps * (abs(par_hat) + 1)
  names(h) <- names(par_hat)

  p <- npm

  total_evals <- as.integer(
    1 + 2 * p + 4 * choose(p, 2)
  )

  H <- matrix(
    NA_real_,
    p,
    p,
    dimnames = list(names(par_hat), names(par_hat))
  )

  f0 <- NULL

  save_checkpoint <- function() {

    if (is.null(checkpoint_file)) return(invisible(NULL))

    saveRDS(
      list(
        par = par_hat,
        step = h,
        f0 = f0,
        H = H,
        hess_eps = hess_eps
      ),
      checkpoint_file
    )

    invisible(NULL)
  }

  resumed <- FALSE

  if (!is.null(checkpoint_file) &&
      file.exists(checkpoint_file)) {

    cp <- readRDS(checkpoint_file)

    compatible <-
      length(cp$par) == length(par_hat) &&
      max(abs(cp$par - par_hat)) < 1e-12 &&
      length(cp$step) == length(h) &&
      max(abs(cp$step - h)) < 1e-12

    if (compatible) {

      H <- cp$H
      f0 <- cp$f0
      resumed <- TRUE

      cat(
        "      checkpoint trovato: riprendo il calcolo precedente\n"
      )

    } else {

      stop(
        "il checkpoint esistente appartiene a parametri ",
        "o step diversi. Rinominalo o cancellalo."
      )
    }
  }

  diag_done <- sum(is.finite(diag(H)))
  off_done <- sum(is.finite(H[lower.tri(H)]))

  completed_evals <-
    if (!is.null(f0)) 1L else 0L

  completed_evals <-
    completed_evals +
    2L * diag_done +
    4L * off_done

  eval_count <- completed_evals

  next_report <- (
    floor(eval_count / progress_every) + 1L
  ) * progress_every

  report_progress <- function(force = FALSE) {

    if (force || eval_count >= next_report) {

      cat(
        sprintf(
          "      Hessiana: %d/%d valutazioni (%.1f%%)\n",
          eval_count,
          total_evals,
          100 * eval_count / total_evals
        )
      )

      if (!force) {
        next_report <<- (
          floor(eval_count / progress_every) + 1L
        ) * progress_every
      }
    }

    invisible(NULL)
  }

  nll <- function(par) {

    val <- nll_raw(par)

    eval_count <<- eval_count + 1L

    report_progress()

    val
  }

  if (is.null(f0)) {

    f0 <- nll(par_hat)

    save_checkpoint()
  }

  ll_diff <- f0 - fit$negloglik

  cat(
    sprintf(
      "      likelihood: %.8f | fit: %.8f | diff: %.3e\n",
      f0,
      fit$negloglik,
      ll_diff
    )
  )

  if (abs(ll_diff) > 1e-3) {
    stop(
      "la likelihood ricostruita non coincide con fit$negloglik."
    )
  }

  cat(
    sprintf(
      "[3/6] Hessiana: %d parametri, %d valutazioni totali\n",
      p,
      total_evals
    )
  )

  if (resumed) {
    report_progress(force = TRUE)
  }

  for (j in seq_len(p)) {

    if (is.finite(H[j, j])) next

    xp <- par_hat
    xm <- par_hat

    xp[j] <- xp[j] + h[j]
    xm[j] <- xm[j] - h[j]

    fp <- nll(xp)
    fm <- nll(xm)

    H[j, j] <- (
      fp - 2 * f0 + fm
    ) / h[j]^2

    save_checkpoint()
  }

  pairs_since_save <- 0L

  if (p > 1L) {

    for (j in seq_len(p - 1L)) {

      for (k in (j + 1L):p) {

        if (is.finite(H[k, j])) next

        xpp <- par_hat
        xpm <- par_hat
        xmp <- par_hat
        xmm <- par_hat

        xpp[j] <- xpp[j] + h[j]
        xpp[k] <- xpp[k] + h[k]

        xpm[j] <- xpm[j] + h[j]
        xpm[k] <- xpm[k] - h[k]

        xmp[j] <- xmp[j] - h[j]
        xmp[k] <- xmp[k] + h[k]

        xmm[j] <- xmm[j] - h[j]
        xmm[k] <- xmm[k] - h[k]

        fpp <- nll(xpp)
        fpm <- nll(xpm)
        fmp <- nll(xmp)
        fmm <- nll(xmm)

        hij <- (
          fpp - fpm - fmp + fmm
        ) / (4 * h[j] * h[k])

        H[j, k] <- hij
        H[k, j] <- hij

        pairs_since_save <- pairs_since_save + 1L

        if (pairs_since_save >= 20L) {
          save_checkpoint()
          pairs_since_save <- 0L
        }
      }

      save_checkpoint()
    }
  }

  H <- (H + t(H)) / 2

  save_checkpoint()

  report_progress(force = TRUE)

  cat("[4/6] controllo Hessiana e inversione\n")

  eig <- eigen(
    H,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  pd <- all(is.finite(eig)) && min(eig) > 0

  cond <- if (pd) {
    max(eig) / min(eig)
  } else {
    Inf
  }

  fit2 <- fit

  fit2$hessian <- H

  fit2$hessian_diagnostics <- list(
    positive_definite = pd,
    min_eigenvalue = min(eig),
    max_eigenvalue = max(eig),
    condition_number = cond,
    hess_eps = hess_eps,
    hessian_step = h,
    total_likelihood_evaluations = total_evals,
    negloglik_recomputed = f0,
    negloglik_fit = fit$negloglik,
    negloglik_difference = ll_diff
  )

  cat(
    sprintf(
      "      min eigenvalue: %.6e\n",
      min(eig)
    )
  )

  cat(
    sprintf(
      "      max eigenvalue: %.6e\n",
      max(eig)
    )
  )

  cat(
    sprintf(
      "      positive definite: %s\n",
      pd
    )
  )

  cat(
    sprintf(
      "      condition number: %.6e\n",
      cond
    )
  )

  if (!pd) {

    warning(
      "Hessiana non positiva definita. ",
      "Mi fermo prima di vcov e Wald."
    )

    return(fit2)
  }

  Vopt <- tryCatch(
    solve(H),
    error = function(e) NULL
  )

  if (is.null(Vopt)) {
    warning("impossibile invertire la Hessiana.")
    return(fit2)
  }

  Vopt <- (Vopt + t(Vopt)) / 2

  dimnames(Vopt) <- list(
    names(par_hat),
    names(par_hat)
  )

  dopt <- diag(Vopt)

  se_opt <- setNames(
    ifelse(
      is.finite(dopt) & dopt >= 0,
      sqrt(dopt),
      NA_real_
    ),
    names(par_hat)
  )

  cat("[5/6] delta method sui parametri naturali\n")

  decode_model <- function(par) {

    par <- as.double(par)
    names(par) <- names(par_hat)

    pos <- 0L

    td <- .growth_decode_thresholds(
      par = par,
      pos = pos,
      cm1 = cm1,
      namvar = namvar,
      threshold_identification = go$threshold_identification,
      monotone_tau = monotone_tau
    )

    tau <- td$tau
    pos <- td$pos

    lambda <- matrix(
      0,
      nrow = nvar,
      ncol = Q,
      dimnames = list(namvar, factor_names)
    )

    for (q in seq_len(Q)) {
      for (k in which(loadmask[, q] == 1L)) {
        pos <- pos + 1L
        lambda[k, q] <- par[pos]
      }
    }

    gamma <- setNames(
      numeric(nvar),
      namvar
    )

    gamma_time <- matrix(
      0,
      nrow = nvar,
      ncol = n,
      dimnames = list(
        namvar,
        paste0("t", seq_len(n))
      )
    )

    for (k in seq_len(nvar)) {
      pos <- pos + 1L
      gamma[k] <- par[pos]
    }

    gamma_time[,] <- gamma

    if (n > 1L) {
      for (tt in 2:n) {
        for (k in which(intercept_mask[, tt])) {
          pos <- pos + 1L
          gamma_time[k, tt] <- par[pos]
        }
      }
    }

    theta <- matrix(
      1,
      nrow = nvar,
      ncol = n,
      dimnames = list(
        namvar,
        paste0("t", seq_len(n))
      )
    )

    if (!theta_time) {

      for (k in seq_len(nvar)) {

        if (cm1[k] > 1L) {

          pos <- pos + 1L
          theta[k, ] <- exp(par[pos])
        }
      }

    } else {

      for (k in seq_len(nvar)) {

        if (cm1[k] > 1L) {

          for (tt in seq_len(n)) {

            pos <- pos + 1L
            theta[k, tt] <- exp(par[pos])
          }
        }
      }
    }

    beta <- matrix(
      0,
      nrow = Q,
      ncol = p_fix
    )

    for (r in seq_len(p_fix)) {

      for (q in seq_len(Q)) {

        if (r == 1L) next

        pos <- pos + 1L
        beta[q, r] <- par[pos]
      }
    }

    rho <- setNames(
      numeric(Q),
      factor_names
    )

    if (delta_free) {

      for (q in seq_len(Q)) {

        pos <- pos + 1L
        rho[q] <- par[pos]
      }
    }

    cv <- .growth_decode_covariance(
      par = par,
      pos = pos,
      Q = Q,
      p_ran = p_ran,
      n = n,
      u_names = u_names,
      identification = "baseline_std"
    )

    pos <- cv$pos

    sigmak <- NULL

    if (re_item) {

      sigmak <- setNames(
        numeric(nvar),
        namvar
      )

      for (k in seq_len(nvar)) {

        pos <- pos + 1L
        sigmak[k] <- exp(par[pos])
      }
    }

    if (pos != npm) {
      stop(
        "decode_model: pos = ",
        pos,
        ", npm = ",
        npm
      )
    }

    baseline_var <-
      diag(cv$Sigma_u)[seq_len(Q)] +
      cv$zeta_var[1L, ]

    names(baseline_var) <- factor_names

    list(
      tau = tau,
      lambda = lambda,
      gamma = gamma,
      gamma_time = gamma_time,
      theta = theta,
      beta = beta,
      rho = rho,
      Lu = cv$Lu,
      Sigma_u = cv$Sigma_u,
      zeta_var = cv$zeta_var,
      sigmak = sigmak,
      baseline_var = baseline_var
    )
  }

  natural_fun <- function(par) {

    d <- decode_model(par)

    out <- numeric(0)

    for (k in seq_len(nvar)) {

      th <- d$tau[[k]]
      out[names(th)] <- as.numeric(th)
    }

    for (q in seq_len(Q)) {

      for (k in which(loadmask[, q] == 1L)) {

        out[
          sprintf(
            "lambda[%s,%s]",
            namvar[k],
            factor_names[q]
          )
        ] <- d$lambda[k, q]
      }
    }

    for (k in seq_len(nvar)) {

      out[
        sprintf("gamma[%s]", namvar[k])
      ] <- d$gamma[k]

      for (tt in seq_len(n)) {

        out[
          sprintf(
            "gamma[%s,t%d]",
            namvar[k],
            tt
          )
        ] <- d$gamma_time[k, tt]
      }
    }

    for (k in seq_len(nvar)) {

      for (tt in seq_len(n)) {

        out[
          sprintf(
            "theta[%s,t%d]",
            namvar[k],
            tt
          )
        ] <- d$theta[k, tt]
      }
    }

    for (r in seq_len(p_fix)) {

      for (q in seq_len(Q)) {

        nm0 <- fixed_names[[q]][r]

        if (r == 1L) {

          out[
            sprintf(
              "alpha0[%s]",
              factor_names[q]
            )
          ] <- d$beta[q, r]

        } else if (r == 2L) {

          out[
            sprintf(
              "alpha1[%s]",
              factor_names[q]
            )
          ] <- d$beta[q, r]

          out[
            sprintf(
              "muS[%s]",
              factor_names[q]
            )
          ] <- d$beta[q, r]

        } else {

          out[
            sprintf(
              "beta[%s,%s]",
              nm0,
              factor_names[q]
            )
          ] <- d$beta[q, r]
        }
      }
    }

    for (q in seq_len(Q)) {

      out[
        sprintf(
          "rho[%s]",
          factor_names[q]
        )
      ] <- d$rho[q]
    }

    for (cc in seq_len(nrow(d$Lu))) {

      for (rr in cc:nrow(d$Lu)) {

        out[
          sprintf(
            "Lu[%s,%s]",
            rownames(d$Lu)[rr],
            colnames(d$Lu)[cc]
          )
        ] <- d$Lu[rr, cc]
      }
    }

    for (r in seq_len(nrow(d$Sigma_u))) {

      out[
        sprintf(
          "Var_u[%s]",
          rownames(d$Sigma_u)[r]
        )
      ] <- d$Sigma_u[r, r]
    }

    if (nrow(d$Sigma_u) >= 2L) {

      for (cc in seq_len(nrow(d$Sigma_u) - 1L)) {

        for (rr in (cc + 1L):nrow(d$Sigma_u)) {

          out[
            sprintf(
              "Cov_u[%s,%s]",
              rownames(d$Sigma_u)[rr],
              colnames(d$Sigma_u)[cc]
            )
          ] <- d$Sigma_u[rr, cc]
        }
      }
    }

    for (tt in seq_len(n)) {

      for (q in seq_len(Q)) {

        out[
          sprintf(
            "zeta_var[t%d,%s]",
            tt,
            factor_names[q]
          )
        ] <- d$zeta_var[tt, q]
      }
    }

    if (!is.null(d$sigmak)) {

      for (k in seq_len(nvar)) {

        out[
          sprintf(
            "sigmak[%s]",
            namvar[k]
          )
        ] <- d$sigmak[k]
      }
    }

    for (q in seq_len(Q)) {

      out[
        sprintf(
          "baseline_var[%s]",
          factor_names[q]
        )
      ] <- d$baseline_var[q]
    }

    out
  }

  nat0 <- natural_fun(par_hat)

  common <- intersect(
    names(nat0),
    names(fit$par_natural)
  )

  if (!length(common)) {
    stop("nessun parametro naturale comune trovato.")
  }

  max_nat_diff <- max(
    abs(
      nat0[common] -
      fit$par_natural[common]
    ),
    na.rm = TRUE
  )

  cat(
    sprintf(
      "      controllo decoder: max diff = %.3e\n",
      max_nat_diff
    )
  )

  if (!is.finite(max_nat_diff) ||
      max_nat_diff > 1e-7) {

    stop(
      "il decoder post-hoc non coincide con fit$par_natural."
    )
  }

  numeric_jacobian <- function(par, f, eps) {

    f0j <- f(par)

    J <- matrix(
      NA_real_,
      nrow = length(f0j),
      ncol = length(par),
      dimnames = list(
        names(f0j),
        names(par)
      )
    )

    for (j in seq_along(par)) {

      hj <- eps * (abs(par[j]) + 1)

      pp <- par
      pm <- par

      pp[j] <- pp[j] + hj
      pm[j] <- pm[j] - hj

      J[, j] <- (
        f(pp) - f(pm)
      ) / (2 * hj)
    }

    J
  }

  J <- numeric_jacobian(
    par_hat,
    natural_fun,
    jac_eps
  )

  Vnat <- J %*% Vopt %*% t(J)
  Vnat <- (Vnat + t(Vnat)) / 2

  dnat <- diag(Vnat)

  se_nat <- setNames(
    ifelse(
      is.finite(dnat) & dnat >= 0,
      sqrt(dnat),
      NA_real_
    ),
    rownames(Vnat)
  )

  se_full <- setNames(
    rep(
      NA_real_,
      length(fit$par_natural)
    ),
    names(fit$par_natural)
  )

  se_full[names(se_nat)] <- se_nat

  fixed_nat_names <- character(0)

  for (k in seq_len(nvar)) {

    m <- cm1[k]

    if (m == 1L) {

      fixed_nat_names <- c(
        fixed_nat_names,
        sprintf("tau[%s]_c1", namvar[k])
      )

    } else {

      j0 <- (m + 1L) %/% 2L

      fixed_nat_names <- c(
        fixed_nat_names,
        sprintf(
          "tau[%s]_c%d",
          namvar[k],
          j0
        ),
        sprintf(
          "tau[%s]_c%d",
          namvar[k],
          j0 + 1L
        )
      )
    }

    if (m == 1L) {

      fixed_nat_names <- c(
        fixed_nat_names,
        sprintf(
          "theta[%s,t%d]",
          namvar[k],
          seq_len(n)
        )
      )
    }
  }

  fixed_nat_names <- c(
    fixed_nat_names,
    sprintf("alpha0[%s]", factor_names),
    sprintf("baseline_var[%s]", factor_names)
  )

  fixed_nat_names <- intersect(
    unique(fixed_nat_names),
    names(se_full)
  )

  se_full[fixed_nat_names] <- 0

  make_wald <- function(est, se) {

    zcrit <- qnorm(0.975)

    z <- ifelse(
      is.finite(se) & se > 0,
      est / se,
      NA_real_
    )

    p_raw <- ifelse(
      is.finite(z),
      2 * pnorm(
        abs(z),
        lower.tail = FALSE
      ),
      NA_real_
    )

    data.frame(
      name = names(est),
      est = as.numeric(est),
      se = as.numeric(se),
      z = as.numeric(z),
      p_raw = as.numeric(p_raw),
      ci_low = ifelse(
        is.finite(se) & se > 0,
        est - zcrit * se,
        ifelse(se == 0, est, NA_real_)
      ),
      ci_high = ifelse(
        is.finite(se) & se > 0,
        est + zcrit * se,
        ifelse(se == 0, est, NA_real_)
      ),
      row.names = NULL,
      check.names = FALSE
    )
  }

  wald_all <- make_wald(
    fit$par_natural,
    se_full
  )

  boundary <- grepl(
    paste0(
      "^(Var_u\\[|",
      "zeta_var\\[|",
      "sigmak\\[|",
      "theta\\[|",
      "baseline_var\\[|",
      "Lu\\[)"
    ),
    wald_all$name
  )

  wald_all$wald_zero_ok <-
    is.finite(wald_all$se) &
    wald_all$se > 0 &
    !boundary

  wald_all$p <- ifelse(
    wald_all$wald_zero_ok,
    wald_all$p_raw,
    NA_real_
  )

  wald_nat <- wald_all[
    is.finite(wald_all$se) &
      wald_all$se > 0,
    ,
    drop = FALSE
  ]

  fit2$vcov_opt <- Vopt
  fit2$se_opt <- se_opt
  fit2$se <- se_full

  fit2$vcov_model_natural <- Vnat
  fit2$jacobian_model_natural <- J

  fit2$wald_all <- wald_all
  fit2$wald_nat <- wald_nat

  fit2$hessian_diagnostics$max_natural_decoder_difference <-
    max_nat_diff

  cat("[6/6] completato\n")

  fit2
}

growth_tt_get_agh_re_full_equiv_grid <- function() {

  sym_info <- native_symbol(
    "get_growth_re_agh_full_equiv_grid_info"
  )

  sym_copy <- native_symbol(
    "get_growth_re_agh_full_equiv_grid_copy"
  )

  info <- native_call(
    sym_info,
    ready0 = integer(1L),
    fit_id0 = integer(1L),
    dim_w0 = integer(1L),
    dim_z0 = integer(1L),
    nsub0 = integer(1L),
    nquad_u0 = integer(1L),
    nquad_zeta0 = integer(1L)
  )

  if (info$ready0 != 1L) {
    stop("cache AGH_RE_FULL_EQUIV non pronta.")
  }

  dz <- as.integer(info$dim_z0)
  ns <- as.integer(info$nsub0)

  if (dz < 1L || ns < 1L) {
    stop("dimensioni cache AGH_RE_FULL_EQUIV non valide.")
  }

  out <- native_call(
    sym_copy,
    mode_out = double(dz * ns),
    B_out = double(dz * dz * ns),
    logdet_w_out = double(ns),
    ready0 = integer(1L),
    fit_id0 = integer(1L),
    dim_w0 = integer(1L),
    dim_z0 = integer(1L),
    nsub0 = integer(1L),
    nquad_u0 = integer(1L),
    nquad_zeta0 = integer(1L)
  )

  if (out$ready0 != 1L) {
    stop("copia cache AGH_RE_FULL_EQUIV fallita.")
  }

  list(
    ready = TRUE,
    type = "AGH_RE_FULL_EQUIV",
    fit_id0_c = as.integer(out$fit_id0),
    dim_w = as.integer(out$dim_w0),
    dim_z = as.integer(out$dim_z0),
    nsub = as.integer(out$nsub0),
    nquad_u = as.integer(out$nquad_u0),
    nquad_zeta = as.integer(out$nquad_zeta0),
    mode = array(
      out$mode_out,
      dim = c(
        out$dim_z0,
        out$nsub0
      )
    ),
    B = array(
      out$B_out,
      dim = c(
        out$dim_z0,
        out$dim_z0,
        out$nsub0
      )
    ),
    logdet_w = out$logdet_w_out
  )
}

growth_tt_set_agh_re_full_equiv_score_grid <- function(grid) {

  if (is.null(grid)) {
    stop("grid mancante.")
  }

  if (!identical(grid$type, "AGH_RE_FULL_EQUIV")) {
    stop("grid non AGH_RE_FULL_EQUIV.")
  }

  sym <- native_symbol(
    "set_growth_re_agh_full_equiv_score_grid"
  )

  z <- native_call(
    sym,
    mode_in = as.double(grid$mode),
    B_in = as.double(grid$B),
    logdet_w_in = as.double(grid$logdet_w),
    dim_w0 = as.integer(grid$dim_w),
    dim_z0 = as.integer(grid$dim_z),
    nsub0 = as.integer(grid$nsub),
    nquad_u0 = as.integer(grid$nquad_u),
    nquad_zeta0 = as.integer(grid$nquad_zeta),
    fit_id0 = as.integer(grid$fit_id0_c),
    ready0 = integer(1L)
  )

  if (z$ready0 != 1L) {
    stop("caricamento griglia nella cache score fallito.")
  }

  invisible(TRUE)
}

growth_add_analytic_hessian <- function(
  fit,
  hess_eps = 1e-4,
  grad_eps = 1e-5,
  verify_score = TRUE,
  score_tol = 1e-4,
  nll_tol = 1e-6,
  asym_tol = 1e-3,
  verbose = TRUE
) {
  if (is.null(fit$c_inputs)) stop("fit$c_inputs mancante.")
  if (is.null(fit$agh_grid)) stop("fit$agh_grid mancante.")
  if (is.null(fit$par_opt)) stop("fit$par_opt mancante.")
  if (is.null(fit$negloglik)) stop("fit$negloglik mancante.")

  ci <- fit$c_inputs
  grid <- fit$agh_grid
  go <- ci$growth_options
  dims <- if (!is.null(ci$dims)) ci$dims else fit$dims
  meth <- ci$methInteg

  if (is.null(go)) stop("c_inputs$growth_options mancante.")
  if (is.null(dims)) stop("dimensioni del modello mancanti.")

  if (!meth %in% c("AGH_FULL_EQUIV", "AGH_RE_FULL_EQUIV")) {
    stop("la funzione richiede AGH_FULL_EQUIV o AGH_RE_FULL_EQUIV.")
  }

  required_functions <- ".growth_set_fortran_options"
  if (identical(meth, "AGH_FULL_EQUIV")) {
    required_functions <- c(
      required_functions,
      "growth_tt_set_agh_full_equiv_grid"
    )
  } else {
    required_functions <- c(
      required_functions,
      "growth_tt_set_agh_re_full_equiv_score_grid"
    )
  }

  missing_functions <- required_functions[
    !vapply(required_functions, exists, logical(1L), mode = "function", inherits = TRUE)
  ]

  if (length(missing_functions)) {
    stop(
      "funzioni mancanti. eseguire prima source('growthFun.R'): ",
      paste(missing_functions, collapse = ", ")
    )
  }

  score_name <- if (identical(meth, "AGH_FULL_EQUIV")) {
    "loglik_growth_tt_agh_full_equiv_score_q1"
  } else {
    "loglik_growth_re_tt_agh_full_equiv_score_q1"
  }

  if (!native_has(score_name)) {
    stop("routine ", score_name, " non caricata.")
  }

  if (!go$identification %in% c("marker", "mixed_marker", "baseline_std")) {
    stop("identification deve essere 'marker', 'mixed_marker' o 'baseline_std'.")
  }

  if (!identical(go$growth_shape, "rho")) {
    stop("growth_shape deve essere 'rho'.")
  }

  if (as.integer(dims$Q) != 1L) {
    stop("lo score analitico richiede Q=1.")
  }

  if (as.integer(dims$p_ran) != 2L) {
    stop("lo score analitico richiede p_ran=2.")
  }

  if (!isTRUE(grid$ready)) stop("griglia AGH non pronta.")
  if (!identical(grid$type, meth)) {
    stop("tipo di griglia AGH incompatibile.")
  }

  if (as.integer(grid$fit_id0_c) != as.integer(fit$fit_id0_c)) {
    stop("fit_id della griglia diverso dal fit.")
  }

  if (as.integer(grid$fit_id0_c) != as.integer(ci$fit_id0_c)) {
    stop("fit_id della griglia diverso da c_inputs.")
  }

  if (as.integer(grid$nsub) != as.integer(ci$nsub0_c)) {
    stop("nsub della griglia incompatibile.")
  }

  if (identical(meth, "AGH_FULL_EQUIV")) {
    if (as.integer(grid$nquad_eta) != as.integer(ci$nquad_eta0_c)) {
      stop("nquad_eta della griglia incompatibile.")
    }
    if (as.integer(grid$nquad_re) != as.integer(ci$nquad_re0_c)) {
      stop("nquad_re della griglia incompatibile.")
    }
    if (as.integer(grid$dim_eta) != as.integer(dims$dim_eta)) {
      stop("dim_eta della griglia incompatibile.")
    }
    if (as.integer(grid$nvar) != as.integer(dims$nvar)) {
      stop("nvar della griglia incompatibile.")
    }
  } else {
    if (as.integer(grid$nquad_u) != as.integer(ci$nquad_u0_c)) {
      stop("nquad_u della griglia incompatibile.")
    }
    if (as.integer(grid$nquad_zeta) != as.integer(ci$nquad_zeta0_c)) {
      stop("nquad_zeta della griglia incompatibile.")
    }
  }

  par_hat <- fit$par_opt
  nm <- names(par_hat)
  npm <- length(par_hat)

  if (is.null(nm)) stop("fit$par_opt deve avere i nomi dei parametri.")
  if (npm != as.integer(ci$npm0_c)) stop("length(par_opt) != npm0_c.")

  marker_item <- NULL
  marker_lambda <- NA_real_

  if (go$identification %in% c("marker", "mixed_marker")) {
    marker_index <- as.integer(go$marker_index)

    if (
      length(marker_index) != 1L ||
      is.na(marker_index) ||
      marker_index < 1L ||
      marker_index > as.integer(dims$nvar)
    ) {
      stop("marker_index non valido.")
    }

    marker_item <- ci$namvar[marker_index]

    if (length(marker_item) != 1L || is.na(marker_item)) {
      stop("marker item non valido.")
    }

    if (is.null(fit$decode$lambda)) stop("fit$decode$lambda mancante.")

    lambda_mat <- as.matrix(fit$decode$lambda)
    marker_lambda <- as.numeric(lambda_mat[marker_index, 1L])

    if (!is.finite(marker_lambda) || abs(marker_lambda - 1) > 1e-10) {
      stop("il loading marker non risulta fissato a 1.")
    }

    marker_loading_idx <- grep("^lambda.*\\[", nm)

    if (length(marker_loading_idx)) {
      marker_in_par <- grepl(marker_item, nm[marker_loading_idx], fixed = TRUE)
      if (any(marker_in_par)) {
        stop("il loading marker fissato compare erroneamente in par_opt.")
      }
    }
  }

  rho_idx <- grep("^rho\\[", nm)
  rho_fixed <- go$rho_fixed

  if (is.null(rho_fixed)) {
    if (as.integer(ci$delta_free0_c) != 1L) {
      stop("rho libero ma delta_free0_c != 1.")
    }
    if (length(rho_idx) != 1L) {
      stop("rho libero non compare esattamente una volta in par_opt.")
    }
    rho_type <- "libero"
  } else {
    rho_fixed <- as.numeric(rho_fixed)
    if (length(rho_fixed) != 1L || !is.finite(rho_fixed)) {
      stop("rho_fixed non valido.")
    }
    if (as.integer(ci$delta_free0_c) != 0L) {
      stop("rho fissato ma delta_free0_c != 0.")
    }
    if (length(rho_idx) != 0L) {
      stop("rho fissato ma compare in par_opt.")
    }
    rho_type <- paste0("fisso a ", rho_fixed)
  }

  subject_specific_design <- if (!is.null(ci$subject_specific_design0_c)) {
    as.integer(ci$subject_specific_design0_c)
  } else {
    as.integer(as.integer(ci$p_fix0_c) > 0L)
  }

  if (
    length(subject_specific_design) != 1L ||
    is.na(subject_specific_design) ||
    !subject_specific_design %in% c(0L, 1L)
  ) {
    stop("subject_specific_design0_c non valido.")
  }

  p_fix0_score <- if (subject_specific_design == 1L) {
    abs(as.integer(ci$p_fix0_c))
  } else {
    -abs(as.integer(ci$p_fix0_c))
  }

  design_nsub <- if (subject_specific_design == 1L) {
    as.integer(ci$nsub0_c)
  } else {
    1L
  }

  expected_X <- design_nsub * as.integer(ci$n0_c) * as.integer(ci$Q0_c) * abs(p_fix0_score)
  expected_Z <- design_nsub * as.integer(ci$n0_c) * as.integer(ci$Q0_c) * as.integer(ci$p_ran0_c)
  expected_delta <- design_nsub * as.integer(ci$n0_c) * as.integer(ci$Q0_c)

  if (length(ci$X0_c) != expected_X) stop("lunghezza X0_c incompatibile.")
  if (length(ci$Z0_c) != expected_Z) stop("lunghezza Z0_c incompatibile.")
  if (length(ci$delta_mask0_c) != expected_delta) stop("lunghezza delta_mask0_c incompatibile.")

  if (identical(meth, "AGH_FULL_EQUIV")) {
    growth_tt_set_agh_full_equiv_grid(grid)
  } else {
    growth_tt_set_agh_re_full_equiv_score_grid(grid)
  }

  sym_score <- native_symbol(score_name)

  set_options <- function() {
    .growth_set_fortran_options(
      identification = go$identification,
      threshold_identification = go$threshold_identification,
      growth_shape = go$growth_shape,
      intercept_mask = go$intercept_mask,
      gamma_dep = go$gamma_dep,
      x_slope_col = go$x_slope_col,
      nvar = as.integer(dims$nvar),
      n = as.integer(dims$n),
      Q = as.integer(dims$Q),
      p_fix = as.integer(dims$p_fix),
      rho_fixed = go$rho_fixed,
      marker_index = go$marker_index
    )
  }

  score_eval <- function(par) {
    set_options()

    if (identical(meth, "AGH_FULL_EQUIV")) {
      z <- native_call(
        sym_score,
        b0 = as.double(par),
        Y0 = as.integer(ci$Y0_c),
        nvar0 = as.integer(ci$nvar0_c),
        nsub0 = as.integer(ci$nsub0_c),
        nobs0 = as.integer(ci$nobs0_c),
        n0 = as.integer(ci$n0_c),
        nmes0 = as.integer(ci$nmes0_c),
        c0 = as.integer(ci$c0_c),
        cm10 = as.integer(ci$cm10_c),
        ntre0 = as.integer(ci$ntre0_c),
        nquad_eta0 = as.integer(ci$nquad_eta0_c),
        nodes_eta0 = as.double(ci$nodes_eta0_c),
        weights_eta0 = as.double(ci$weights_eta0_c),
        nquad_re0 = as.integer(ci$nquad_re0_c),
        nodes_re0 = as.double(ci$nodes_re0_c),
        weights_re0 = as.double(ci$weights_re0_c),
        nitem0 = as.integer(ci$nitem0_c),
        npm0 = as.integer(ci$npm0_c),
        re_item0 = as.integer(ci$re_item0_c),
        theta_time0 = as.integer(ci$theta_time0_c),
        delta_free0 = as.integer(ci$delta_free0_c),
        monotone_tau0 = as.integer(ci$monotone_tau0_c),
        update_gh0 = 0L,
        fit_id0 = as.integer(grid$fit_id0_c),
        Q0 = as.integer(ci$Q0_c),
        loadmask0 = as.integer(ci$loadmask0_c),
        p_fix0 = as.integer(p_fix0_score),
        p_ran0 = as.integer(ci$p_ran0_c),
        X0 = as.double(ci$X0_c),
        Z0 = as.double(ci$Z0_c),
        delta_mask0 = as.double(ci$delta_mask0_c),
        delta_x_col0 = as.integer(ci$delta_x_col0_c),
        delta_z_col0 = as.integer(ci$delta_z_col0_c),
        vrais = double(1L),
        score0 = double(npm),
        ok0 = integer(1L)
      )
    } else {
      z <- native_call(
        sym_score,
        b0 = as.double(par),
        Y0 = as.integer(ci$Y0_c),
        nvar0 = as.integer(ci$nvar0_c),
        nsub0 = as.integer(ci$nsub0_c),
        nobs0 = as.integer(ci$nobs0_c),
        n0 = as.integer(ci$n0_c),
        nmes0 = as.integer(ci$nmes0_c),
        c0 = as.integer(ci$c0_c),
        cm10 = as.integer(ci$cm10_c),
        ntre0 = as.integer(ci$ntre0_c),
        nquad_u0 = as.integer(ci$nquad_u0_c),
        nodes_u0 = as.double(ci$nodes_u0_c),
        weights_u0 = as.double(ci$weights_u0_c),
        nquad_zeta0 = as.integer(ci$nquad_zeta0_c),
        nodes_zeta0 = as.double(ci$nodes_zeta0_c),
        weights_zeta0 = as.double(ci$weights_zeta0_c),
        nquad_re0 = as.integer(ci$nquad_re0_c),
        nodes_re0 = as.double(ci$nodes_re0_c),
        weights_re0 = as.double(ci$weights_re0_c),
        nitem0 = as.integer(ci$nitem0_c),
        npm0 = as.integer(ci$npm0_c),
        re_item0 = as.integer(ci$re_item0_c),
        theta_time0 = as.integer(ci$theta_time0_c),
        delta_free0 = as.integer(ci$delta_free0_c),
        monotone_tau0 = as.integer(ci$monotone_tau0_c),
        update_gh0 = 0L,
        fit_id0 = as.integer(grid$fit_id0_c),
        Q0 = as.integer(ci$Q0_c),
        loadmask0 = as.integer(ci$loadmask0_c),
        p_fix0 = as.integer(p_fix0_score),
        p_ran0 = as.integer(ci$p_ran0_c),
        X0 = as.double(ci$X0_c),
        Z0 = as.double(ci$Z0_c),
        delta_mask0 = as.double(ci$delta_mask0_c),
        delta_x_col0 = as.integer(ci$delta_x_col0_c),
        delta_z_col0 = as.integer(ci$delta_z_col0_c),
        vrais = double(1L),
        score0 = double(npm),
        ok0 = integer(1L)
      )
    }

    if (length(z$ok0) != 1L || z$ok0 != 1L) {
      stop("score analitico Fortran fallito, ok0=", z$ok0)
    }

    if (!is.finite(z$vrais) || abs(z$vrais) >= 1e8) {
      stop("NLL restituita dallo score non valida.")
    }

    if (length(z$score0) != npm || any(!is.finite(z$score0))) {
      stop("score analitico non valido.")
    }

    list(
      nll = as.numeric(z$vrais),
      score = setNames(as.numeric(z$score0), nm)
    )
  }

  e0 <- score_eval(par_hat)
  nll_fit <- as.numeric(fit$negloglik)
  nll_abs_diff <- abs(e0$nll - nll_fit)
  nll_rel_diff <- nll_abs_diff / max(1, abs(nll_fit))

  if (!is.finite(nll_rel_diff) || nll_rel_diff > nll_tol) {
    stop(
      "NLL score e NLL del fit non coincidono. diff relativa = ",
      format(nll_rel_diff, scientific = TRUE)
    )
  }

  grad_max_abs_diff <- NA_real_
  grad_rel_diff <- NA_real_

  if (isTRUE(verify_score)) {
    hg <- grad_eps * (abs(par_hat) + 1)
    grad_numeric <- setNames(numeric(npm), nm)

    for (j in seq_len(npm)) {
      pp <- par_hat
      pm <- par_hat
      pp[j] <- pp[j] + hg[j]
      pm[j] <- pm[j] - hg[j]
      grad_numeric[j] <- (score_eval(pp)$nll - score_eval(pm)$nll) / (2 * hg[j])
    }

    grad_diff <- e0$score - grad_numeric
    grad_max_abs_diff <- max(abs(grad_diff))
    grad_scale <- max(1, max(abs(e0$score)), max(abs(grad_numeric)))
    grad_rel_diff <- grad_max_abs_diff / grad_scale

    if (!is.finite(grad_rel_diff) || grad_rel_diff > score_tol) {
      worst <- order(abs(grad_diff), decreasing = TRUE)
      worst <- worst[seq_len(min(10L, length(worst)))]

      print(data.frame(
        parameter = nm[worst],
        analytic = e0$score[worst],
        numeric = grad_numeric[worst],
        abs_diff = abs(grad_diff[worst]),
        row.names = NULL,
        check.names = FALSE
      ))

      stop(
        "score analitico incompatibile con il gradiente numerico. diff relativa = ",
        format(grad_rel_diff, scientific = TRUE)
      )
    }
  }

  h <- hess_eps * (abs(par_hat) + 1)
  Hraw <- matrix(
    NA_real_, npm, npm,
    dimnames = list(nm, nm)
  )

  t0 <- proc.time()

  for (j in seq_len(npm)) {
    if (isTRUE(verbose)) {
      cat(sprintf("Hessiana: parametro %d/%d\r", j, npm))
      flush.console()
    }

    pp <- par_hat
    pm <- par_hat
    pp[j] <- pp[j] + h[j]
    pm[j] <- pm[j] - h[j]

    sp <- score_eval(pp)$score
    sm <- score_eval(pm)$score
    Hraw[, j] <- (sp - sm) / (2 * h[j])
  }

  if (isTRUE(verbose)) cat("\n")

  tt <- proc.time() - t0
  elapsed <- unname(tt["elapsed"])

  raw_asym_abs <- max(abs(Hraw - t(Hraw)))
  raw_scale <- max(1, max(abs(Hraw)))
  raw_asym_rel <- raw_asym_abs / raw_scale

  if (!is.finite(raw_asym_rel) || raw_asym_rel > asym_tol) {
    stop(
      "Hessiana raw troppo asimmetrica. asimmetria relativa = ",
      format(raw_asym_rel, scientific = TRUE)
    )
  }

  H <- (Hraw + t(Hraw)) / 2
  eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  positive_definite <- all(is.finite(eig)) && min(eig) > 0

  V <- NULL
  se_opt <- setNames(rep(NA_real_, npm), nm)

  if (positive_definite) {
    V <- solve(H)
    V <- (V + t(V)) / 2
    dimnames(V) <- list(nm, nm)

    dV <- diag(V)
    if (any(!is.finite(dV)) || any(dV < 0)) {
      stop("diagonale della matrice di covarianza non valida.")
    }

    se_opt <- setNames(sqrt(dV), nm)
  }

  diagnostics <- list(
    method = meth,
    identification = go$identification,
    marker_item = marker_item,
    marker_lambda = marker_lambda,
    rho = rho_type,
    npar = npm,
    subject_specific_design = subject_specific_design,
    p_fix0_saved = as.integer(ci$p_fix0_c),
    p_fix0_score = p_fix0_score,
    nll_score = e0$nll,
    nll_fit = nll_fit,
    nll_absolute_difference = nll_abs_diff,
    nll_relative_difference = nll_rel_diff,
    verify_score = isTRUE(verify_score),
    score_max_abs_difference = grad_max_abs_diff,
    score_relative_difference = grad_rel_diff,
    score_max_abs_at_optimum = max(abs(e0$score)),
    hessian_raw_max_asymmetry = raw_asym_abs,
    hessian_raw_relative_asymmetry = raw_asym_rel,
    min_eigenvalue = min(eig),
    max_eigenvalue = max(eig),
    positive_definite = positive_definite,
    elapsed = elapsed
  )

  fit$hessian <- H
  fit$hessian_raw <- Hraw
  fit$hessian_method <- "analytic_fixed"
  fit$hessian_step <- setNames(h, nm)
  fit$score_at_optimum <- e0$score
  fit$vcov_opt <- V
  fit$se_opt <- se_opt
  fit$hessian_eigenvalues <- eig
  fit$hessian_positive_definite <- positive_definite
  fit$analytic_fixed_diagnostics <- diagnostics

  if (!is.null(fit$timing)) {
    fit$timing$hess_elapsed <- elapsed
    fit$timing$hess_user <- unname(tt["user.self"])
    fit$timing$hess_system <- unname(tt["sys.self"])

    if (!is.null(fit$timing$opt_elapsed)) {
      fit$timing$total_elapsed <- fit$timing$opt_elapsed + elapsed
    }
  }

  if (isTRUE(verbose)) {
    cat("positive definite:", positive_definite, "\n")
    cat("min eigenvalue:", format(min(eig), scientific = TRUE), "\n")
    cat("elapsed:", sprintf("%.3f", elapsed), "sec\n")
  }

  fit
}
add_sandwich_postfit <- function(fit,
                                 eps_score = 1e-5,
                                 progress = TRUE) {

  if (is.null(fit$par_opt)) {
    stop("fit$par_opt mancante.")
  }

  if (is.null(fit$hessian)) {
    stop("fit$hessian mancante. il modello deve essere stato stimato con hessian.")
  }

  if (is.null(fit$c_inputs)) {
    stop("fit$c_inputs mancante.")
  }

  meth <- if (!is.null(fit$structure$methInteg)) {
    fit$structure$methInteg
  } else {
    fit$c_inputs$methInteg
  }

  if (meth %in% c("AGH_RE_NESTED", "AGH_RE_FULL_EQUIV")) {
    stop(
      "growth_tt_individual_derivatives non supporta ancora ",
      meth,
      "."
    )
  }

  der <- growth_tt_individual_derivatives(
    fit = fit,
    data = NULL,
    namvar = NULL,
    compute_score = TRUE,
    compute_hessian = FALSE,
    eps_score = eps_score,
    agh_fixed_grid = TRUE,
    progress = progress,
    attach_to_fit = FALSE
  )

  S <- der$score_ind

  if (is.null(S) || any(!is.finite(S))) {
    stop("score individuali non validi.")
  }

  A <- (fit$hessian + t(fit$hessian)) / 2
  B <- crossprod(S)

  if (!all(dim(A) == dim(B))) {
    stop("dimensioni incompatibili tra hessian e score.")
  }

  Ainv <- tryCatch(
    solve(A),
    error = function(e) NULL
  )

  if (is.null(Ainv)) {
    stop("impossibile invertire la hessian.")
  }

  V_sand <- Ainv %*% B %*% Ainv
  V_sand <- (V_sand + t(V_sand)) / 2

  nm <- names(fit$par_opt)

  dimnames(V_sand) <- list(nm, nm)

  se_sand <- sqrt(pmax(diag(V_sand), 0))
  names(se_sand) <- nm

  fit$vcov_opt_model <- fit$vcov_opt
  fit$se_opt_model <- fit$se_opt

  fit$vcov_opt_sandwich <- V_sand
  fit$se_opt_sandwich <- se_sand

  fit$score_ind <- S

  fit$sandwich <- list(
    bread = A,
    meat = B,
    vcov = V_sand,
    se = se_sand
  )

  fit$sandwich_checks <- list(
    loglik_individual_sum = sum(der$loglik_i),
    loglik_fit = fit$loglik,
    loglik_difference = sum(der$loglik_i) - fit$loglik,
    max_abs_total_score = max(abs(colSums(S))),
    finite_scores = all(is.finite(S)),
    min_eigen_bread = min(
      eigen(
        A,
        symmetric = TRUE,
        only.values = TRUE
      )$values
    )
  )

  fit$vcov_opt <- V_sand
  fit$se_opt <- se_sand

  fit
}
