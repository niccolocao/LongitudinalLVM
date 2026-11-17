.growth_default_design <- function(nsub, n, Q,
                                   subject_specific_design = TRUE) {
  nsub <- as.integer(nsub)
  n <- as.integer(n)
  Q <- as.integer(Q)
  subject_specific_design <- isTRUE(subject_specific_design)
  design_nsub <- if (subject_specific_design) nsub else 1L

  X <- array(0, dim = c(design_nsub, n, Q, 2L))
  Z <- array(0, dim = c(design_nsub, n, Q, 2L))
  delta_mask <- array(0, dim = c(design_nsub, n, Q))
  time <- seq.int(0, n - 1L)
  time_mat <- matrix(rep(time, each = design_nsub),
                     nrow = design_nsub, ncol = n)

  for (q in seq_len(Q)) {
    X[, , q, 1L] <- 1
    X[, , q, 2L] <- time_mat
    Z[, , q, 1L] <- 1
    Z[, , q, 2L] <- time_mat
    delta_mask[, n, q] <- 1
  }

  list(
    p_fix = 2L,
    p_ran = 2L,
    X = X,
    Z = Z,
    delta_mask = delta_mask,
    delta_x_col = 2L,
    delta_z_col = 2L,
    subject_specific_design = subject_specific_design,
    design_nsub = design_nsub
  )
}

.growth_lowlevel_design <- function(nsub0, n0, Q0,
                                    p_fix0 = NULL, p_ran0 = NULL,
                                    X0 = NULL, Z0 = NULL,
                                    delta_mask0 = NULL,
                                    delta_x_col0 = NULL,
                                    delta_z_col0 = NULL,
                                    subject_specific_design0 = NULL) {
  nsub0 <- as.integer(nsub0)
  n0 <- as.integer(n0)
  Q0 <- as.integer(Q0)

  if (!is.null(subject_specific_design0)) {
    if (length(subject_specific_design0) != 1L || is.na(subject_specific_design0)) {
      stop("subject_specific_design0 deve essere TRUE/FALSE oppure NULL.")
    }
    subject_specific_design0 <- isTRUE(subject_specific_design0)
  }

  if (!is.null(p_fix0)) {
    p_fix0 <- as.integer(p_fix0)
    if (p_fix0 < 0L) {
      if (is.null(subject_specific_design0)) subject_specific_design0 <- FALSE
      p_fix0 <- abs(p_fix0)
    }
  }

  if (is.null(X0) || is.null(Z0)) {
    default_subject_specific <- if (is.null(subject_specific_design0)) TRUE else subject_specific_design0
    d <- .growth_default_design(
      nsub0, n0, Q0,
      subject_specific_design = default_subject_specific
    )
    if (is.null(X0)) X0 <- d$X
    if (is.null(Z0)) Z0 <- d$Z
    if (is.null(p_fix0)) p_fix0 <- d$p_fix
    if (is.null(p_ran0)) p_ran0 <- d$p_ran
    if (is.null(delta_mask0)) delta_mask0 <- d$delta_mask
    if (is.null(delta_x_col0)) delta_x_col0 <- d$delta_x_col
    if (is.null(delta_z_col0)) delta_z_col0 <- d$delta_z_col
    if (is.null(subject_specific_design0)) {
      subject_specific_design0 <- d$subject_specific_design
    }
  }

  infer_p <- function(A, p, label) {
    if (!is.null(p)) return(abs(as.integer(p)))
    d <- dim(A)
    if (!is.null(d) && length(d) == 4L &&
        d[1L] %in% c(1L, nsub0) && d[2L] == n0 && d[3L] == Q0) {
      return(as.integer(d[4L]))
    }
    if (Q0 == 1L && !is.null(d) && length(d) == 3L &&
        d[1L] %in% c(1L, nsub0) && d[2L] == n0) {
      return(as.integer(d[3L]))
    }
    if (!is.null(d) && length(d) == 3L && d[1L] == n0 && d[2L] == Q0) {
      return(as.integer(d[3L]))
    }
    if (is.matrix(A) && nrow(A) %in% c(n0, nsub0 * n0)) {
      return(as.integer(ncol(A)))
    }
    stop(label, ": specificare p oppure passare un array/matrice con dimensioni esplicite.")
  }

  p_fix0 <- infer_p(X0, p_fix0, "X0")
  p_ran0 <- infer_p(Z0, p_ran0, "Z0")
  if (p_fix0 < 1L || p_ran0 < 1L) stop("p_fix0 e p_ran0 devono essere positivi.")

  infer_design_nsub <- function(A, p, label, delta = FALSE) {
    d <- dim(A)
    mult <- if (delta) 1L else p
    common_len <- n0 * Q0 * mult
    full_len <- nsub0 * common_len

    if (!is.null(d)) {
      if (!delta && length(d) == 4L && d[2L] == n0 && d[3L] == Q0 && d[4L] == p &&
          d[1L] %in% c(1L, nsub0)) return(as.integer(d[1L]))
      if (!delta && Q0 == 1L && length(d) == 3L && d[2L] == n0 && d[3L] == p &&
          d[1L] %in% c(1L, nsub0)) return(as.integer(d[1L]))
      if (!delta && length(d) == 3L && identical(as.integer(d), c(n0, Q0, p))) return(1L)
      if (!delta && is.matrix(A) && ncol(A) == p) {
        if (nrow(A) == n0) return(1L)
        if (nrow(A) == nsub0 * n0) return(nsub0)
      }
      if (delta && length(d) == 3L && d[2L] == n0 && d[3L] == Q0 &&
          d[1L] %in% c(1L, nsub0)) return(as.integer(d[1L]))
      if (delta && is.matrix(A) && Q0 == 1L) {
        if (all(dim(A) == c(1L, n0))) return(1L)
        if (all(dim(A) == c(nsub0, n0))) return(nsub0)
      }
    }

    len <- length(A)
    if (nsub0 > 1L && len == common_len) return(1L)
    if (len == full_len) return(nsub0)
    if (nsub0 == 1L && len == common_len) return(1L)
    stop(label, ": lunghezza o dimensioni incompatibili.")
  }

  ns_x <- infer_design_nsub(X0, p_fix0, "X0")
  ns_z <- infer_design_nsub(Z0, p_ran0, "Z0")
  ns_d <- if (is.null(delta_mask0)) integer(0) else {
    infer_design_nsub(delta_mask0, 1L, "delta_mask0", delta = TRUE)
  }

  if (is.null(subject_specific_design0)) {
    subject_specific_design0 <- any(c(ns_x, ns_z, ns_d) == nsub0) && nsub0 > 1L
  }
  design_nsub <- if (subject_specific_design0) nsub0 else 1L

  as_design_array <- function(A, p, label) {
    src_nsub <- infer_design_nsub(A, p, label)
    d <- dim(A)

    if (!is.null(d) && length(d) == 4L) {
      out <- array(as.double(A), dim = c(src_nsub, n0, Q0, p), dimnames = dimnames(A))
    } else if (!is.null(d) && Q0 == 1L && length(d) == 3L &&
               d[1L] == src_nsub && d[2L] == n0 && d[3L] == p) {
      out <- array(as.double(A), dim = c(src_nsub, n0, 1L, p))
    } else if (!is.null(d) && length(d) == 3L &&
               identical(as.integer(d), c(n0, Q0, p))) {
      out <- array(0, dim = c(1L, n0, Q0, p))
      out[1L, , , ] <- A
    } else if (is.matrix(A) && nrow(A) == n0 && ncol(A) == p) {
      if (Q0 != 1L) stop(label, " matrix n x p richiede Q0=1 oppure un array n x Q x p.")
      out <- array(0, dim = c(1L, n0, 1L, p))
      out[1L, , 1L, ] <- A
    } else if (is.matrix(A) && nrow(A) == nsub0 * n0 && ncol(A) == p) {
      if (Q0 != 1L) stop(label, " matrix nsub*n x p richiede Q0=1.")
      out <- array(0, dim = c(nsub0, n0, 1L, p))
      rr <- 0L
      for (i in seq_len(nsub0)) for (t in seq_len(n0)) {
        rr <- rr + 1L
        out[i, t, 1L, ] <- A[rr, ]
      }
    } else {
      out <- array(as.double(A), dim = c(src_nsub, n0, Q0, p))
    }

    if (src_nsub == design_nsub) return(out)
    if (src_nsub == 1L && design_nsub == nsub0) {
      return(out[rep(1L, nsub0), , , , drop = FALSE])
    }
    if (src_nsub == nsub0 && design_nsub == 1L) {
      ref <- out[rep(1L, nsub0), , , , drop = FALSE]
      if (!isTRUE(all.equal(unname(out), unname(ref), tolerance = 0))) {
        stop(label, " varia tra soggetti ma subject_specific_design0=FALSE.")
      }
      return(out[1L, , , , drop = FALSE])
    }
    stop(label, ": impossibile armonizzare la dimensione soggetto.")
  }

  as_delta_array <- function(A) {
    if (is.null(A)) return(array(0, dim = c(design_nsub, n0, Q0)))
    src_nsub <- infer_design_nsub(A, 1L, "delta_mask0", delta = TRUE)
    d <- dim(A)
    if (!is.null(d) && length(d) == 3L) {
      out <- array(as.double(A), dim = c(src_nsub, n0, Q0))
    } else if (is.matrix(A) && Q0 == 1L && all(dim(A) == c(src_nsub, n0))) {
      out <- array(as.double(A), dim = c(src_nsub, n0, 1L))
    } else {
      out <- array(as.double(A), dim = c(src_nsub, n0, Q0))
    }

    if (src_nsub == design_nsub) return(out)
    if (src_nsub == 1L && design_nsub == nsub0) {
      return(out[rep(1L, nsub0), , , drop = FALSE])
    }
    if (src_nsub == nsub0 && design_nsub == 1L) {
      ref <- out[rep(1L, nsub0), , , drop = FALSE]
      if (!isTRUE(all.equal(unname(out), unname(ref), tolerance = 0))) {
        stop("delta_mask0 varia tra soggetti ma subject_specific_design0=FALSE.")
      }
      return(out[1L, , , drop = FALSE])
    }
    stop("delta_mask0: impossibile armonizzare la dimensione soggetto.")
  }

  X0 <- as_design_array(X0, p_fix0, "X0")
  Z0 <- as_design_array(Z0, p_ran0, "Z0")
  delta_mask0 <- as_delta_array(delta_mask0)

  if (is.null(delta_x_col0)) delta_x_col0 <- 0L
  if (is.null(delta_z_col0)) delta_z_col0 <- 0L
  delta_x_col0 <- as.integer(delta_x_col0)
  delta_z_col0 <- as.integer(delta_z_col0)
  if (delta_x_col0 < 0L || delta_x_col0 > p_fix0) stop("delta_x_col0 fuori range.")
  if (delta_z_col0 < 0L || delta_z_col0 > p_ran0) stop("delta_z_col0 fuori range.")

  if (any(!is.finite(X0))) stop("X0 contiene valori non finiti.")
  if (any(!is.finite(Z0))) stop("Z0 contiene valori non finiti.")
  if (any(!is.finite(delta_mask0))) stop("delta_mask0 contiene valori non finiti.")

  list(

    p_fix0 = as.integer(if (subject_specific_design0) p_fix0 else -p_fix0),
    p_ran0 = as.integer(p_ran0),
    p_fix = as.integer(p_fix0),
    p_ran = as.integer(p_ran0),
    X0 = as.double(X0),
    Z0 = as.double(Z0),
    delta_mask0 = as.double(delta_mask0),
    delta_x_col0 = delta_x_col0,
    delta_z_col0 = delta_z_col0,
    subject_specific_design0 = as.integer(subject_specific_design0),
    design_nsub = as.integer(design_nsub)
  )
}

.growth_prepare_lowlevel_call <- function(nsub0, n0, Q0,
                                          p_fix0, p_ran0,
                                          X0, Z0, delta_mask0,
                                          delta_x_col0, delta_z_col0,
                                          subject_specific_design0,
                                          prepared_design0 = FALSE) {
  if (!isTRUE(prepared_design0)) {
    return(.growth_lowlevel_design(
      nsub0, n0, Q0, p_fix0, p_ran0, X0, Z0,
      delta_mask0, delta_x_col0, delta_z_col0,
      subject_specific_design0
    ))
  }

  if (is.null(p_fix0) || is.null(p_ran0) ||
      is.null(X0) || is.null(Z0) || is.null(delta_mask0) ||
      is.null(subject_specific_design0)) {
    stop("prepared_design0=TRUE richiede tutti gli input di design gi\u00e0 preparati.")
  }

  nsub0 <- as.integer(nsub0)
  n0 <- as.integer(n0)
  Q0 <- as.integer(Q0)
  p_fix <- abs(as.integer(p_fix0))
  p_ran <- as.integer(p_ran0)
  subject_specific <- isTRUE(as.logical(subject_specific_design0))
  design_nsub <- if (subject_specific) nsub0 else 1L

  if (p_fix < 1L || p_ran < 1L) {
    stop("p_fix0 e p_ran0 devono essere positivi.")
  }
  if (length(X0) != design_nsub * n0 * Q0 * p_fix) {
    stop("X0 preparato ha lunghezza incompatibile.")
  }
  if (length(Z0) != design_nsub * n0 * Q0 * p_ran) {
    stop("Z0 preparato ha lunghezza incompatibile.")
  }
  if (length(delta_mask0) != design_nsub * n0 * Q0) {
    stop("delta_mask0 preparato ha lunghezza incompatibile.")
  }

  list(
    p_fix0 = as.integer(if (subject_specific) p_fix else -p_fix),
    p_ran0 = p_ran,
    p_fix = p_fix,
    p_ran = p_ran,
    X0 = X0,
    Z0 = Z0,
    delta_mask0 = delta_mask0,
    delta_x_col0 = as.integer(delta_x_col0),
    delta_z_col0 = as.integer(delta_z_col0),
    subject_specific_design0 = as.integer(subject_specific),
    design_nsub = as.integer(design_nsub)
  )
}

.growth_build_design <- function(data_sorted, ID, occ, IDs, occ_vals, Q,
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
  nsub <- length(IDs)
  n <- length(occ_vals)
  nrow_long <- nrow(data_sorted)
  ii <- match(data_sorted[[ID]], IDs)
  tt <- match(data_sorted[[occ]], occ_vals)

  if (anyNA(ii) || anyNA(tt)) stop("impossibile indicizzare ID/occasioni per il design growth.")
  key <- paste(ii, tt, sep = ":")
  if (anyDuplicated(key)) stop("il design growth richiede al massimo una riga per coppia ID-occasione.")

  object_is_subject_specific <- function(obj) {
    if (is.null(obj)) return(FALSE)
    if (is.character(obj)) return(TRUE)
    if (is.list(obj) && !is.data.frame(obj)) {
      return(any(vapply(obj, object_is_subject_specific, logical(1L))))
    }
    d <- dim(obj)
    if (!is.null(d)) {
      if (length(d) == 4L && d[1L] == nsub) return(TRUE)
      if (length(d) == 3L && d[1L] == nsub && d[2L] == n) return(TRUE)
      if (is.matrix(obj) && nrow(obj) %in% c(nsub, nsub * n, nrow_long) && nrow(obj) != n) return(TRUE)
      if (is.matrix(obj) && all(dim(obj) == c(nsub, n))) return(TRUE)
      return(FALSE)
    }
    length(obj) == nsub && nsub != n
  }

  if (is.null(subject_specific_design)) {
    subject_specific_design <- if (isTRUE(simple_growth)) {
      any(vapply(
        list(time_scores, fixed_covariates, delta_mask),
        object_is_subject_specific,
        logical(1L)
      ))
    } else {
      object_is_subject_specific(X) || object_is_subject_specific(Z)
    }
  }
  if (length(subject_specific_design) != 1L || is.na(subject_specific_design)) {
    stop("subject_specific_design deve essere TRUE/FALSE oppure NULL.")
  }
  subject_specific_design <- isTRUE(subject_specific_design)
  design_nsub <- if (subject_specific_design) nsub else 1L

  map_long_matrix <- function(M, label) {
    M <- as.matrix(M)
    storage.mode(M) <- "double"
    p <- ncol(M)
    out <- array(NA_real_, dim = c(design_nsub, n, p),
                 dimnames = list(NULL, NULL, colnames(M)))

    if (subject_specific_design) {
      if (nrow(M) == nrow_long) {
        for (r in seq_len(nrow_long)) out[ii[r], tt[r], ] <- M[r, ]
      } else if (nrow(M) == nsub) {
        for (j in seq_len(n)) out[, j, ] <- M
      } else if (nrow(M) == nsub * n) {
        rr <- 0L
        for (i in seq_len(nsub)) for (j in seq_len(n)) {
          rr <- rr + 1L
          out[i, j, ] <- M[rr, ]
        }
      } else if (nrow(M) == n) {
        for (i in seq_len(nsub)) out[i, , ] <- M
      } else if (nrow(M) == 1L) {
        for (i in seq_len(nsub)) for (j in seq_len(n)) out[i, j, ] <- M[1L, ]
      } else {
        stop(label, " deve avere nrow uguale a nrow(data), nsub, n, 1 oppure nsub*n.")
      }
    } else {
      if (nrow(M) == n) {
        out[1L, , ] <- M
      } else if (nrow(M) == 1L) {
        for (j in seq_len(n)) out[1L, j, ] <- M[1L, ]
      } else if (nrow(M) == nrow_long) {
        for (j in seq_len(n)) {
          rows <- which(tt == j)
          if (!length(rows)) stop(label, ": occasione senza valori.")
          for (cc in seq_len(p)) {
            vals <- unique(M[rows, cc])
            vals <- vals[is.finite(vals)]
            if (length(vals) != 1L) {
              stop(label, " varia tra soggetti all'occasione ", j,
                   "; usare subject_specific_design=TRUE.")
            }
            out[1L, j, cc] <- vals
          }
        }
      } else if (nrow(M) == nsub) {
        for (cc in seq_len(p)) {
          vals <- unique(M[, cc])
          vals <- vals[is.finite(vals)]
          if (length(vals) != 1L) {
            stop(label, " varia tra soggetti; usare subject_specific_design=TRUE.")
          }
          out[1L, , cc] <- vals
        }
      } else if (nrow(M) == nsub * n) {
        full <- array(NA_real_, dim = c(nsub, n, p))
        rr <- 0L
        for (i in seq_len(nsub)) for (j in seq_len(n)) {
          rr <- rr + 1L
          full[i, j, ] <- M[rr, ]
        }
        for (j in seq_len(n)) for (cc in seq_len(p)) {
          vals <- unique(full[, j, cc])
          vals <- vals[is.finite(vals)]
          if (length(vals) != 1L) {
            stop(label, " varia tra soggetti; usare subject_specific_design=TRUE.")
          }
          out[1L, j, cc] <- vals
        }
      } else {
        stop(label, " comune deve avere n, 1, nrow(data), nsub oppure nsub*n righe.")
      }
    }

    if (any(!is.finite(out))) {
      stop(label, " contiene valori mancanti o non finiti nel reticolo richiesto.")
    }
    out
  }

  to_grid <- function(obj, label) {
    if (is.character(obj)) {
      if (!all(obj %in% names(data_sorted))) stop("colonne non trovate in ", label, ".")
      return(map_long_matrix(data_sorted[, obj, drop = FALSE], label))
    }

    d <- dim(obj)
    if (!is.null(d) && length(d) == 3L &&
        identical(as.integer(d[1:2]), c(design_nsub, n))) {
      out <- array(as.double(obj), dim = d, dimnames = dimnames(obj))
      if (any(!is.finite(out))) stop(label, " contiene valori non finiti.")
      return(out)
    }
    if (!subject_specific_design && !is.null(d) && length(d) == 2L && nrow(obj) == n) {
      return(map_long_matrix(obj, label))
    }
    if (is.vector(obj) && !is.list(obj)) obj <- matrix(obj, ncol = 1L)
    map_long_matrix(obj, label)
  }

  to_time <- function(obj, label) {
    if (is.null(obj)) {
      return(matrix(rep(seq.int(0, n - 1L), each = design_nsub),
                    nrow = design_nsub, ncol = n))
    }
    if (is.character(obj) && length(obj) == 1L) {
      g <- to_grid(obj, label)
      return(g[, , 1L])
    }
    if (is.numeric(obj) && is.null(dim(obj))) {
      if (length(obj) == n) {
        return(matrix(rep(obj, each = design_nsub), nrow = design_nsub, ncol = n))
      }
      if (subject_specific_design && length(obj) == nsub) {
        return(matrix(rep(obj, n), nrow = nsub, ncol = n))
      }
      if (!subject_specific_design && length(obj) == 1L) {
        return(matrix(rep(obj, n), nrow = 1L, ncol = n))
      }
      g <- to_grid(obj, label)
      return(g[, , 1L])
    }
    if (is.matrix(obj) && all(dim(obj) == c(design_nsub, n))) {
      out <- matrix(as.double(obj), nrow = design_nsub, ncol = n)
      if (any(!is.finite(out))) stop(label, " contiene valori non finiti.")
      return(out)
    }
    g <- to_grid(obj, label)
    if (dim(g)[3L] != 1L) stop(label, " deve descrivere un solo score temporale.")
    g[, , 1L]
  }

  factor_objects <- function(obj, label) {
    if (is.null(obj)) return(rep(list(NULL), Q))
    is_factor_list <- is.list(obj) && !is.data.frame(obj)
    if (Q == 1L && !is_factor_list) return(list(obj))
    if (isTRUE(same_design) && !is_factor_list) return(rep(list(obj), Q))
    if (is_factor_list) {
      if (isTRUE(same_design) && length(obj) == 1L) return(rep(obj, Q))
      if (length(obj) != Q) stop(label, " deve essere una lista di lunghezza Q.")
      return(obj)
    }
    stop(label, " deve essere una lista per fattore quando same_design=FALSE e Q>1.")
  }

  normalize_factor_design <- function(obj, label) {
    d <- dim(obj)
    if (!is.null(d) && length(d) == 4L &&
        identical(as.integer(d[1:3]), c(design_nsub, n, Q))) {
      p <- as.integer(d[4L])
      out <- array(as.double(obj), dim = d, dimnames = dimnames(obj))
      if (any(!is.finite(out))) stop(label, " contiene valori non finiti.")
      grids <- lapply(seq_len(Q), function(q) {
        dn <- dimnames(out); if (is.null(dn)) dn <- vector("list", 4L)
        array(out[, , q, , drop = FALSE], dim = c(design_nsub, n, p),
              dimnames = list(dn[[1L]], dn[[2L]], dn[[4L]]))
      })
      return(list(array = out, grids = grids, p = p))
    }
    if (Q == 1L && !is.null(d) && length(d) == 3L &&
        identical(as.integer(d[1:2]), c(design_nsub, n))) {
      p <- as.integer(d[3L])
      base <- array(as.double(obj), dim = c(design_nsub, n, p), dimnames = dimnames(obj))
      out <- array(base, dim = c(design_nsub, n, 1L, p))
      return(list(array = out, grids = list(base), p = p))
    }
    if (!subject_specific_design && !is.null(d) && length(d) == 3L &&
        identical(as.integer(d[1:2]), c(n, Q))) {
      p <- as.integer(d[3L])
      out <- array(0, dim = c(1L, n, Q, p))
      out[1L, , , ] <- obj
      grids <- lapply(seq_len(Q), function(q) array(out[1L, , q, ], dim = c(1L, n, p)))
      return(list(array = out, grids = grids, p = p))
    }

    objs <- factor_objects(obj, label)
    grids <- lapply(seq_len(Q), function(q) {
      to_grid(objs[[q]], paste0(label, "[[", q, "]]"))
    })
    pp <- vapply(grids, function(g) dim(g)[3L], integer(1L))
    if (length(unique(pp)) != 1L) stop(label, " deve avere lo stesso numero di colonne per ogni fattore.")
    p <- pp[1L]
    out <- array(0, dim = c(design_nsub, n, Q, p))
    for (q in seq_len(Q)) out[, , q, ] <- grids[[q]]
    list(array = out, grids = grids, p = p)
  }

  normalize_time <- function(obj, label) {
    d <- dim(obj)
    if (!is.null(d) && length(d) == 3L &&
        identical(as.integer(d), c(design_nsub, n, Q))) {
      out <- array(as.double(obj), dim = d, dimnames = dimnames(obj))
      if (any(!is.finite(out))) stop(label, " contiene valori non finiti.")
      return(out)
    }
    if (!subject_specific_design && is.matrix(obj) && all(dim(obj) == c(n, Q))) {
      out <- array(0, dim = c(1L, n, Q))
      out[1L, , ] <- obj
      return(out)
    }

    objs <- factor_objects(obj, label)
    out <- array(0, dim = c(design_nsub, n, Q))
    for (q in seq_len(Q)) out[, , q] <- to_time(objs[[q]], paste0(label, "[[", q, "]]"))
    out
  }

  normalize_names <- function(given, grids, prefix) {
    if (is.null(given)) {
      return(lapply(seq_len(Q), function(q) {
        dn <- dimnames(grids[[q]])
        nm <- if (is.null(dn)) NULL else dn[[3L]]
        if (is.null(nm)) nm <- paste0(prefix, seq_len(dim(grids[[q]])[3L]))
        make.unique(as.character(nm))
      }))
    }
    if (Q == 1L && is.character(given)) given <- list(given)
    if (isTRUE(same_design) && is.character(given)) given <- rep(list(given), Q)
    if (is.list(given) && isTRUE(same_design) && length(given) == 1L) given <- rep(given, Q)
    if (!is.list(given) || length(given) != Q) stop(prefix, " names deve essere una lista di lunghezza Q.")
    for (q in seq_len(Q)) {
      if (length(given[[q]]) != dim(grids[[q]])[3L]) stop("numero di nomi incompatibile per il fattore ", q, ".")
      given[[q]] <- make.unique(as.character(given[[q]]))
    }
    given
  }

  if (isTRUE(simple_growth)) {
    time_arr <- normalize_time(time_scores, "time_scores")

    cov_arr <- NULL
    cov_names <- replicate(Q, character(0), simplify = FALSE)
    if (is.null(fixed_covariates) && !is.null(fixed_names)) {
      stop("fixed_names richiede fixed_covariates in modalit\u00e0 simple_growth.")
    }
    if (!is.null(fixed_covariates)) {
      cov <- normalize_factor_design(fixed_covariates, "fixed_covariates")
      cov_arr <- cov$array
      cov_names <- normalize_names(fixed_names, cov$grids, "cov")
      if (is.null(fixed_names) && is.character(fixed_covariates) &&
          (Q == 1L || isTRUE(same_design))) {
        cov_names <- rep(list(as.character(fixed_covariates)), Q)
      }
    }

    if (is.list(time_interactions) && isTRUE(same_design) && length(time_interactions) == 1L) {
      time_interactions <- rep(time_interactions, Q)
    }
    if (is.list(time_interactions) && length(time_interactions) != Q) {
      stop("time_interactions deve essere una lista di lunghezza Q.")
    }
    int_idx <- vector("list", Q)
    for (q in seq_len(Q)) {
      if (is.null(time_interactions)) {
        int_idx[[q]] <- integer(0)
      } else {
        ti <- if (is.list(time_interactions)) time_interactions[[q]] else time_interactions
        if (is.character(ti)) {
          idxq <- match(ti, cov_names[[q]])
          if (anyNA(idxq)) stop("time_interactions contiene covariate non presenti per il fattore ", q, ".")
          int_idx[[q]] <- idxq
        } else {
          idxq <- as.integer(ti)
          if (any(idxq < 1L | idxq > length(cov_names[[q]]))) stop("indice time_interactions fuori range.")
          int_idx[[q]] <- idxq
        }
      }
    }

    p_cov <- if (is.null(cov_arr)) 0L else dim(cov_arr)[4L]
    p_int <- vapply(int_idx, length, integer(1L))
    if (length(unique(p_int)) != 1L) stop("lo stesso numero di interazioni temporali \u00e8 richiesto per tutti i fattori.")
    p_fix <- 2L + p_cov + p_int[1L]
    p_ran <- 2L
    Xout <- array(0, dim = c(design_nsub, n, Q, p_fix))
    Zout <- array(0, dim = c(design_nsub, n, Q, p_ran))
    fixed_names_out <- vector("list", Q)
    random_names_out <- rep(list(c("intercept", "time")), Q)
    if (!is.null(random_names)) {
      rn <- random_names
      if (Q == 1L && is.character(rn)) rn <- list(rn)
      if (isTRUE(same_design) && is.character(rn)) rn <- rep(list(rn), Q)
      if (is.list(rn) && isTRUE(same_design) && length(rn) == 1L) rn <- rep(rn, Q)
      if (!is.list(rn) || length(rn) != Q || any(lengths(rn) != 2L)) {
        stop("random_names deve contenere due nomi per fattore in modalit\u00e0 simple_growth.")
      }
      random_names_out <- lapply(rn, function(x) make.unique(as.character(x)))
    }

    for (q in seq_len(Q)) {
      Xout[, , q, 1L] <- 1
      Xout[, , q, 2L] <- time_arr[, , q]
      Zout[, , q, 1L] <- 1
      Zout[, , q, 2L] <- time_arr[, , q]
      nm <- c("intercept", "time")
      pos <- 2L
      if (p_cov > 0L) {
        Xout[, , q, pos + seq_len(p_cov)] <- cov_arr[, , q, ]
        nm <- c(nm, cov_names[[q]])
        pos <- pos + p_cov
      }
      if (length(int_idx[[q]])) {
        for (jj in int_idx[[q]]) {
          pos <- pos + 1L
          Xout[, , q, pos] <- cov_arr[, , q, jj] * time_arr[, , q]
          nm <- c(nm, paste0(cov_names[[q]][jj], ":time"))
        }
      }
      fixed_names_out[[q]] <- nm
    }

    if (is.null(delta_mask)) {
      dmask <- array(0, dim = c(design_nsub, n, Q))
      dmask[, n, ] <- 1
    } else {
      dmask <- normalize_time(delta_mask, "delta_mask")
    }

    return(list(
      X = Xout,
      Z = Zout,
      delta_mask = dmask,
      delta_x_col = 2L,
      delta_z_col = 2L,
      time_interactions = int_idx,
      interaction_columns = p_fix - p_int[1L] + seq_len(p_int[1L]),
      p_fix = p_fix,
      p_ran = p_ran,
      fixed_names = fixed_names_out,
      random_names = random_names_out,
      time = time_arr,
      simple_growth = TRUE,
      same_design = isTRUE(same_design),
      subject_specific_design = subject_specific_design,
      design_nsub = design_nsub
    ))
  }

  if (is.null(X) || is.null(Z)) stop("con simple_growth=FALSE devi fornire sia X sia Z.")
  xd <- normalize_factor_design(X, "X")
  zd <- normalize_factor_design(Z, "Z")
  fixed_names_out <- normalize_names(fixed_names, xd$grids, "x")
  random_names_out <- normalize_names(random_names, zd$grids, "z")

  list(
    X = xd$array,
    Z = zd$array,
    delta_mask = array(0, dim = c(design_nsub, n, Q)),
    delta_x_col = 0L,
    delta_z_col = 0L,
    p_fix = xd$p,
    p_ran = zd$p,
    fixed_names = fixed_names_out,
    random_names = random_names_out,
    time = NULL,
    simple_growth = FALSE,
    same_design = isTRUE(same_design),
    subject_specific_design = subject_specific_design,
    design_nsub = design_nsub
  )
}
loglik_growth_tt <- function(b0,
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

  native_call("loglik_growth_tt",
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
     vrais = double(1))$vrais
}

loglik_growth_tt_agh <- function(b0,
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

  native_call("loglik_growth_tt_agh",
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
     vrais = double(1))$vrais
}

loglik_growth_tt_agh_full_equiv <- function(b0,
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

  native_call("loglik_growth_tt_agh_full_equiv",
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
     vrais = double(1))$vrais
}

loglik_growth_re_tt_agh_nested <- function(b0,
                                            Y0, nvar0, nsub0, nobs0, n0,
                                            nmes0, c0, cm10, ntre0,
                                            nquad_u0,
                                            nodes_u0, weights_u0,
                                            nquad_zeta0,
                                            nodes_zeta0, weights_zeta0,
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
  loadmask0 <- matrix(as.integer(loadmask0),
                      nrow = as.integer(nvar0),
                      ncol = as.integer(Q0))

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

  native_call("loglik_growth_re_tt_agh_nested",
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
     nquad_u0 = as.integer(nquad_u0),
     nodes_u0 = as.double(nodes_u0),
     weights_u0 = as.double(weights_u0),
     nquad_zeta0 = as.integer(nquad_zeta0),
     nodes_zeta0 = as.double(nodes_zeta0),
     weights_zeta0 = as.double(weights_zeta0),
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
     vrais = double(1))$vrais
}

loglik_growth_re_tt_agh_full_equiv <- function(b0,
                                                Y0, nvar0, nsub0, nobs0, n0,
                                                nmes0, c0, cm10, ntre0,
                                                nquad_u0,
                                                nodes_u0, weights_u0,
                                                nquad_zeta0,
                                                nodes_zeta0, weights_zeta0,
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
  loadmask0 <- matrix(as.integer(loadmask0),
                      nrow = as.integer(nvar0),
                      ncol = as.integer(Q0))

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

  native_call("loglik_growth_re_tt_agh_full_equiv",
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
     nquad_u0 = as.integer(nquad_u0),
     nodes_u0 = as.double(nodes_u0),
     weights_u0 = as.double(weights_u0),
     nquad_zeta0 = as.integer(nquad_zeta0),
     nodes_zeta0 = as.double(nodes_zeta0),
     weights_zeta0 = as.double(weights_zeta0),
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
     vrais = double(1))$vrais
}
.growth_build_design_legacy <- .growth_build_design

.growth_build_design <- function(data_sorted, ID, occ, IDs, occ_vals, Q,
                                 simple_growth = TRUE,
                                 same_design = TRUE,
                                 subject_specific_design = NULL,
                                 time_scores = NULL,
                                 fixed_covariates = NULL,
                                 time_interactions = NULL,
                                 intercept_covariates = NULL,
                                 slope_covariates = NULL,
                                 X = NULL,
                                 Z = NULL,
                                 delta_mask = NULL,
                                 fixed_names = NULL,
                                 random_names = NULL) {
  use_split <- !is.null(intercept_covariates) || !is.null(slope_covariates)

  nsub <- length(IDs)
  n <- length(occ_vals)
  if (nsub < 1L || n < 2L) stop("design growth non valido.")

  ii <- match(data_sorted[[ID]], IDs)
  tt <- match(data_sorted[[occ]], occ_vals)
  if (anyNA(ii) || anyNA(tt)) stop("impossibile indicizzare ID/occasioni per il design growth.")
  key <- paste(ii, tt, sep = ":")
  if (anyDuplicated(key)) {
    stop("il design growth richiede al massimo una riga per coppia ID-occasione.")
  }

  subject_constant_matrix <- function(vars, label, allow_time_varying = FALSE) {
    vars <- unique(as.character(vars))
    vars <- vars[nzchar(vars)]
    if (!length(vars)) {
      return(matrix(numeric(0), nrow = nsub, ncol = 0L,
                    dimnames = list(as.character(IDs), character(0))))
    }

    miss <- setdiff(vars, names(data_sorted))
    if (length(miss)) stop("covariate non trovate: ", paste(miss, collapse = ", "))

    out <- matrix(NA_real_, nrow = nsub, ncol = length(vars),
                  dimnames = list(as.character(IDs), vars))

    for (cc in seq_along(vars)) {
      v <- vars[cc]
      x <- data_sorted[[v]]
      if (!(is.numeric(x) || is.integer(x) || is.logical(x))) {
        stop("covariata ", v, " deve essere numerica/logica.")
      }

      for (i in seq_len(nsub)) {
        z <- as.numeric(x[ii == i])
        z <- z[is.finite(z)]

        if (!length(z)) {
          if (isTRUE(allow_time_varying)) return(NULL)
          stop("covariata ", v, " non osservata per ID ", as.character(IDs[i]), ".")
        }

        uz <- unique(z)
        if (length(uz) != 1L) {
          if (isTRUE(allow_time_varying)) return(NULL)
          stop("covariata ", v, " non costante entro ID ", as.character(IDs[i]), ".")
        }
        out[i, cc] <- uz
      }
    }
    out
  }

  if (!use_split) {
    fixed_covariates0 <- fixed_covariates

    if (is.character(fixed_covariates) && length(fixed_covariates)) {
      subj_fixed <- subject_constant_matrix(
        fixed_covariates,
        "fixed_covariates",
        allow_time_varying = TRUE
      )
      if (!is.null(subj_fixed)) fixed_covariates0 <- subj_fixed
    }

    return(.growth_build_design_legacy(
      data_sorted, ID, occ, IDs, occ_vals, Q,
      simple_growth = simple_growth,
      same_design = same_design,
      subject_specific_design = subject_specific_design,
      time_scores = time_scores,
      fixed_covariates = fixed_covariates0,
      time_interactions = time_interactions,
      X = X, Z = Z, delta_mask = delta_mask,
      fixed_names = fixed_names, random_names = random_names
    ))
  }

  if (!isTRUE(simple_growth)) {
    stop("intercept_covariates/slope_covariates richiedono simple_growth=TRUE.")
  }
  if (!is.null(fixed_covariates) || !is.null(time_interactions)) {
    stop("usa fixed_covariates/time_interactions oppure intercept_covariates/slope_covariates, non entrambi.")
  }
  if (!is.null(X) || !is.null(Z)) {
    stop("con intercept_covariates/slope_covariates non passare X/Z.")
  }
  if (!is.null(fixed_names)) {
    stop("con covariate separate i nomi sono costruiti automaticamente.")
  }
  if (!is.null(intercept_covariates) && !is.character(intercept_covariates)) {
    stop("intercept_covariates deve essere un vettore di nomi di colonne.")
  }
  if (!is.null(slope_covariates) && !is.character(slope_covariates)) {
    stop("slope_covariates deve essere un vettore di nomi di colonne.")
  }

  i_cov <- unique(as.character(intercept_covariates))
  s_cov <- unique(as.character(slope_covariates))
  i_cov <- i_cov[nzchar(i_cov)]
  s_cov <- s_cov[nzchar(s_cov)]
  all_cov <- unique(c(i_cov, s_cov))

  subj <- subject_constant_matrix(
    all_cov,
    "intercept_covariates/slope_covariates",
    allow_time_varying = FALSE
  )

  if (isFALSE(subject_specific_design)) {
    stop("intercept_covariates/slope_covariates richiedono subject_specific_design=TRUE.")
  }
  subject_specific_design <- TRUE
  design_nsub <- nsub

  if (is.null(time_scores)) {
    time_mat <- matrix(rep(seq.int(0, n - 1L), each = nsub), nrow = nsub, ncol = n)
  } else if (is.numeric(time_scores) && is.null(dim(time_scores)) && length(time_scores) == n) {
    time_mat <- matrix(rep(as.double(time_scores), each = nsub), nrow = nsub, ncol = n)
  } else if (is.matrix(time_scores) && all(dim(time_scores) == c(nsub, n))) {
    time_mat <- matrix(as.double(time_scores), nrow = nsub, ncol = n)
  } else {
    stop("con covariate separate time_scores deve essere NULL, un vettore di lunghezza n o una matrice nsub x n.")
  }
  if (any(!is.finite(time_mat))) stop("time_scores contiene valori non finiti.")

  p_fix <- 2L + length(i_cov) + length(s_cov)
  p_ran <- 2L
  Xout <- array(0, dim = c(nsub, n, Q, p_fix))
  Zout <- array(0, dim = c(nsub, n, Q, p_ran))
  fixed_names_out <- rep(list(character(p_fix)), Q)
  random_names_out <- rep(list(c("intercept", "time")), Q)

  if (!is.null(random_names)) {
    rn <- random_names
    if (Q == 1L && is.character(rn)) rn <- list(rn)
    if (isTRUE(same_design) && is.character(rn)) rn <- rep(list(rn), Q)
    if (is.list(rn) && isTRUE(same_design) && length(rn) == 1L) rn <- rep(rn, Q)
    if (!is.list(rn) || length(rn) != Q || any(lengths(rn) != 2L)) {
      stop("random_names deve contenere due nomi per fattore.")
    }
    random_names_out <- lapply(rn, function(x) make.unique(as.character(x)))
  }

  for (q in seq_len(Q)) {
    Xout[, , q, 1L] <- 1
    Xout[, , q, 2L] <- time_mat
    Zout[, , q, 1L] <- 1
    Zout[, , q, 2L] <- time_mat

    pos <- 2L
    if (length(i_cov)) {
      for (v in i_cov) {
        pos <- pos + 1L
        Xout[, , q, pos] <- subj[, v]
      }
    }
    if (length(s_cov)) {
      for (v in s_cov) {
        pos <- pos + 1L
        Xout[, , q, pos] <- subj[, v]
      }
    }
    fixed_names_out[[q]] <- c("intercept", "time", sprintf("I:%s", i_cov), sprintf("S:%s", s_cov))
  }

  if (is.null(delta_mask)) {
    dmask <- array(0, dim = c(nsub, n, Q))
    dmask[, n, ] <- 1
  } else if (is.numeric(delta_mask) && is.null(dim(delta_mask)) && length(delta_mask) == n) {
    dmask <- array(0, dim = c(nsub, n, Q))
    for (q in seq_len(Q)) {
      dmask[, , q] <- matrix(rep(as.double(delta_mask), each = nsub), nrow = nsub, ncol = n)
    }
  } else if (!is.null(dim(delta_mask)) &&
             length(dim(delta_mask)) == 3L &&
             identical(as.integer(dim(delta_mask)), c(nsub, n, Q))) {
    dmask <- array(as.double(delta_mask), dim = c(nsub, n, Q))
  } else if (Q == 1L && is.matrix(delta_mask) && all(dim(delta_mask) == c(nsub, n))) {
    dmask <- array(as.double(delta_mask), dim = c(nsub, n, 1L))
  } else {
    stop("delta_mask incompatibile con il design soggetto-specifico.")
  }
  if (any(!is.finite(dmask))) stop("delta_mask contiene valori non finiti.")

  list(
    X = Xout,
    Z = Zout,
    delta_mask = dmask,
    delta_x_col = 2L,
    delta_z_col = 2L,
    p_fix = p_fix,
    p_ran = p_ran,
    fixed_names = fixed_names_out,
    random_names = random_names_out,
    time = array(rep(time_mat, Q), dim = c(nsub, n, Q)),
    simple_growth = TRUE,
    same_design = isTRUE(same_design),
    subject_specific_design = TRUE,
    design_nsub = nsub,
    subject_covariates = subj
  )
}

.growth_free_threshold_idx <- function(m, threshold_identification) {
  if (m <= 1L) return(integer(0))
  if (threshold_identification == "first_two") return(if (m >= 3L) 3L:m else integer(0))
  j0 <- (m + 1L) %/% 2L
  j1 <- j0 + 1L
  c(if (j0 > 1L) seq.int(j0 - 1L, 1L) else integer(0),
    if (j1 < m) seq.int(j1 + 1L, m) else integer(0))
}

.growth_decode_thresholds <- function(par, pos, cm1, namvar, threshold_identification, monotone_tau) {
  tau <- vector("list", length(cm1))
  names(tau) <- namvar
  for (k in seq_along(cm1)) {
    m <- cm1[k]
    if (m == 1L) {
      tau[[k]] <- setNames(0, sprintf("tau[%s]_c1", namvar[k]))
      next
    }
    if (threshold_identification == "first_two") {
      th <- rep(NA_real_, m)
      th[1:2] <- c(0, 1)
      idx <- .growth_free_threshold_idx(m, threshold_identification)
      if (length(idx)) {
        u <- par[pos + seq_along(idx)]
        pos <- pos + length(idx)
        if (isTRUE(monotone_tau)) th[idx] <- 1 + cumsum(exp(u)) else th[idx] <- u
      }
    } else {
      j0 <- (m + 1L) %/% 2L
      j1 <- j0 + 1L
      th <- rep(NA_real_, m)
      th[j0] <- 0
      th[j1] <- 1
      lower <- if (j0 > 1L) seq.int(j0 - 1L, 1L) else integer(0)
      upper <- if (j1 < m) seq.int(j1 + 1L, m) else integer(0)
      if (length(lower)) {
        u <- par[pos + seq_along(lower)]
        pos <- pos + length(lower)
        if (isTRUE(monotone_tau)) th[lower] <- -cumsum(exp(u)) else th[lower] <- u
      }
      if (length(upper)) {
        u <- par[pos + seq_along(upper)]
        pos <- pos + length(upper)
        if (isTRUE(monotone_tau)) th[upper] <- 1 + cumsum(exp(u)) else th[upper] <- u
      }
    }
    names(th) <- sprintf("tau[%s]_c%d", namvar[k], seq_len(m))
    tau[[k]] <- th
  }
  list(tau = tau, pos = pos)
}

.growth_normalize_intercept_mask <- function(x, namvar, occ_vals) {
  nvar <- length(namvar)
  n <- length(occ_vals)
  out <- matrix(FALSE, nrow = nvar, ncol = n, dimnames = list(namvar, as.character(occ_vals)))
  if (is.null(x)) return(out)
  if (is.matrix(x)) {
    if (!all(dim(x) == c(nvar, n))) stop("intercept_time_free matrix deve avere nvar x n colonne.")
    out[,] <- as.logical(x)
  } else if (is.character(x) && is.null(names(x))) {
    if (!all(x %in% namvar)) stop("intercept_time_free contiene item non presenti.")
    if (n > 1L) out[x, 2:n] <- TRUE
  } else if (is.list(x)) {
    if (is.null(names(x)) || any(!names(x) %in% namvar)) stop("intercept_time_free list deve essere nominata con gli item.")
    for (nm in names(x)) {
      v <- x[[nm]]
      jj <- if (is.character(v)) match(v, as.character(occ_vals)) else as.integer(v)
      if (anyNA(jj) || any(jj < 1L | jj > n)) stop("occasioni non valide in intercept_time_free per ", nm, ".")
      out[nm, jj] <- TRUE
    }
  } else {
    stop("intercept_time_free deve essere NULL, matrice logica, vettore di item o lista nominata.")
  }
  out[, 1L] <- FALSE
  out
}

.growth_select_dep <- function(loadmask, eligible) {
  Q <- ncol(loadmask)
  eligible <- as.integer(eligible)
  if (length(eligible) < Q) return(NULL)
  cmb <- if (Q == 1L) matrix(eligible, nrow = 1L) else combn(eligible, Q)
  if (Q == 1L) cmb <- matrix(eligible, nrow = 1L)
  for (j in seq_len(ncol(cmb))) {
    idx <- cmb[, j]
    A <- t(loadmask[idx, , drop = FALSE])
    if (qr(A)$rank == Q) return(as.integer(idx))
  }
  NULL
}

.growth_make_gamma_dep <- function(loadmask, mask, identification) {
  Q <- ncol(loadmask)
  n <- ncol(mask)
  dep <- matrix(0L, nrow = Q, ncol = n)
  base <- .growth_select_dep(loadmask, seq_len(nrow(loadmask)))
  if (is.null(base)) stop("non esiste un insieme di intercette dipendenti di rango Q.")
  dep[, 1L] <- base
  if (identification == "effects" && n > 1L) {
    for (tt in 2:n) {
      vv <- which(mask[, tt])
      if (!length(vv)) next
      if (length(vv) <= Q) {
        stop("con effects coding, al tempo ", tt, " servono pi\u00f9 di Q intercette variabili; altrimenti il vincolo di somma non consente variazione.")
      }
      dd <- .growth_select_dep(loadmask, vv)
      if (is.null(dd)) stop("intercette variabili al tempo ", tt, " non consentono vincoli di rango Q.")
      dep[, tt] <- dd
    }
  }
  dep
}

.growth_resolve_marker <- function(marker_item, namvar, loadmask, factor_names) {
  Q <- ncol(loadmask)

  if (is.null(marker_item)) {
    stop("con identification='marker' o 'mixed_marker' specificare marker_item.")
  }

  x <- marker_item
  if (!is.null(names(x))) {
    if (!all(factor_names %in% names(x))) {
      stop("marker_item nominato deve contenere un valore per ogni factor_name.")
    }
    x <- x[factor_names]
  }

  if (length(x) == 1L && Q == 1L) {
    x <- rep(x, Q)
  }
  if (length(x) != Q) {
    stop("marker_item deve avere lunghezza Q, una variabile marker per fattore.")
  }

  if (is.character(x)) {
    idx <- match(x, namvar)
    if (anyNA(idx)) {
      stop("marker_item non presente in namvar: ", paste(x[is.na(idx)], collapse = ", "))
    }
  } else if (is.numeric(x) || is.integer(x)) {
    idx <- as.integer(x)
    if (any(!is.finite(x)) || any(idx != x) || any(idx < 1L | idx > length(namvar))) {
      stop("marker_item numerico deve contenere indici interi validi in namvar.")
    }
  } else {
    stop("marker_item deve essere un nome di variabile o un indice intero.")
  }

  if (Q > 1L && anyDuplicated(idx)) {
    stop("con Q > 1 ogni fattore deve avere un marker_item distinto.")
  }

  for (q in seq_len(Q)) {
    if (loadmask[idx[q], q] != 1L) {
      stop(
        "marker_item '", namvar[idx[q]], "' non carica sul fattore '",
        factor_names[q], "'."
      )
    }
  }

  setNames(as.integer(idx), factor_names)
}

.growth_set_fortran_options <- function(identification, threshold_identification, growth_shape,
                                        intercept_mask, gamma_dep, x_slope_col, nvar, n, Q, p_fix,
                                        rho_fixed = NULL, marker_index = NULL, zero_time = rep.int(1L, Q)) {
  if (is.null(zero_time)) zero_time <- rep.int(1L, Q)
  if (!native_has("set_growth_tt_options")) {
    stop("simbolo set_growth_tt_options non caricato: ricompila il Fortran flessibile.")
  }

  ident0 <- switch(
    identification,
    effects = 0L,
    baseline_std = 1L,
    mixed_marker = 2L,
    marker = 3L,
    stop("identification non riconosciuta: ", identification)
  )

  if (identification %in% c("marker", "mixed_marker")) {
    if (is.null(marker_index) || length(marker_index) != Q) {
      stop("marker_index deve avere lunghezza Q con identification='marker' o 'mixed_marker'.")
    }
    marker_index <- as.integer(marker_index)
    if (anyNA(marker_index) || any(marker_index < 1L | marker_index > nvar)) {
      stop("marker_index contiene indici non validi.")
    }
  } else {
    marker_index <- integer(Q)
  }

  if (is.null(rho_fixed)) {
    rho_fixed_flag0 <- 0L
    rho_fixed0 <- numeric(Q)
  } else {
    if (growth_shape != "rho") stop("rho_fixed \u00e8 disponibile solo con growth_shape='rho'.")
    rho_fixed0 <- as.double(rho_fixed)
    if (!length(rho_fixed0) %in% c(1L, Q) || any(!is.finite(rho_fixed0))) {
      stop("rho_fixed deve essere uno scalare finito oppure un vettore finito di lunghezza Q.")
    }
    if (length(rho_fixed0) == 1L) rho_fixed0 <- rep(rho_fixed0, Q)
    rho_fixed_flag0 <- 1L
  }

  ans <- native_call(
    "set_growth_tt_options",
    ident0 = as.integer(ident0),
    threshold_id0 = as.integer(threshold_identification == "central_01"),
    delta_mode0 = as.integer(growth_shape == "rho"),
    nvar0 = as.integer(nvar), n0 = as.integer(n), Q0 = as.integer(Q), pfix0 = as.integer(p_fix),
    gamma_time_free0 = as.integer(intercept_mask),
    gamma_dep0 = as.integer(gamma_dep),
    x_slope_col0 = as.integer(x_slope_col),
    marker_item0 = as.integer(marker_index),
    rho_fixed_flag0 = as.integer(rho_fixed_flag0),
    rho_fixed0 = as.double(rho_fixed0),
    zero_time0 = as.integer(zero_time),
    ok0 = integer(1L)
  )
  if (ans$ok0 != 1L) stop("set_growth_tt_options ha rifiutato la configurazione.")
  invisible(TRUE)
}

.growth_decode_covariance <- function(par, pos, Q, p_ran, n, u_names, identification,
                                      zero_time = rep.int(1L, Q)) {
  dimU <- Q * p_ran
  Lu <- matrix(0, dimU, dimU, dimnames = list(u_names, u_names))
  if (identification %in% c("effects", "marker", "mixed_marker")) {
    for (col in seq_len(dimU)) for (row in col:dimU) {
      pos <- pos + 1L
      Lu[row, col] <- if (row == col) exp(par[pos]) else par[pos]
    }
    Sigma_u <- Lu %*% t(Lu)
    zeta_var <- matrix(NA_real_, n, Q)
    for (tt in seq_len(n)) for (q in seq_len(Q)) {
      pos <- pos + 1L
      zeta_var[tt, q] <- exp(par[pos])
    }
  } else {
    sd <- numeric(dimU)
    for (d in seq_len(dimU)) {
      pos <- pos + 1L
      sd[d] <- if (d <= Q) sqrt(plogis(par[pos])) else exp(par[pos])
    }
    Lc <- matrix(0, dimU, dimU)
    Lc[1, 1] <- 1
    if (dimU > 1L) for (row in 2:dimU) {
      prod <- 1
      for (col in seq_len(row - 1L)) {
        pos <- pos + 1L
        pc <- tanh(par[pos])
        Lc[row, col] <- pc * prod
        prod <- prod * sqrt(max(1e-10, 1 - pc^2))
      }
      Lc[row, row] <- prod
    }
    Lu <- diag(sd, nrow = dimU) %*% Lc
    dimnames(Lu) <- list(u_names, u_names)
    Sigma_u <- Lu %*% t(Lu)
    zeta_var <- matrix(NA_real_, n, Q)
    zeta_var[cbind(zero_time, seq_len(Q))] <- 1 - diag(Sigma_u)[seq_len(Q)]
    for (tt in seq_len(n)) for (q in seq_len(Q)) {
      if (tt == zero_time[q]) next
      pos <- pos + 1L
      zeta_var[tt, q] <- exp(par[pos])
    }
  }
  dimnames(zeta_var) <- list(paste0("t", seq_len(n)), paste0("F", seq_len(Q)))
  list(Lu = Lu, Sigma_u = Sigma_u, zeta_var = zeta_var, pos = pos)
}

#fit an ordinal latent growth model
growth_tt_fit_flexible <- function(data,
                                   namvar,
                                   occ = "occ",
                                   ID = "ID",
                                   loadmask = NULL,
                                   factor_names = NULL,
                                   identification = c("effects", "baseline_std", "marker", "mixed_marker"),
                                   marker_item = NULL,
                                   threshold_identification = c("first_two", "central_01"),
                                   intercept_time_free = NULL,
                                   growth_shape = c("offset", "rho"),
                                   rho_fixed = NULL,
                                   methInteg = c("GH", "QMC", "AGH", "AGH_FULL_EQUIV", "AGH_RE_NESTED", "AGH_RE_FULL_EQUIV"),
                                   nquadgh = 8,
                                   nquadqmc = 2000,
                                   update_gh = 1L,
                                   analytic_score = TRUE,
                                   re_item = TRUE,
                                   theta_time = TRUE,
                                   delta_free = TRUE,
                                   monotone_tau = FALSE,
                                   simple_growth = TRUE,
                                   same_design = TRUE,
                                   subject_specific_design = NULL,
                                   time_scores = NULL,
                                   fixed_covariates = NULL,
                                   time_interactions = NULL,
                                   intercept_covariates = NULL,
                                   slope_covariates = NULL,
                                   X = NULL,
                                   Z = NULL,
                                   delta_mask = NULL,
                                   fixed_names = NULL,
                                   random_names = NULL,
                                   maxiter = 500,
                                   start = NULL,
                                   verbose = FALSE,
                                   compute_hessian = FALSE,
                                    hessian_method = NULL,
                                    hess_eps = 1e-4,
                                   profile = TRUE,
                                   profile_print = TRUE,
                                   parscale = NULL,
                                   debug_par = TRUE,
                                   debug_par_each_eval = TRUE,
                                   debug_par_full = TRUE,
                                   debug_par_max_eval = Inf,
                                   tol_loglik = 1e-7,
                                   tol_par = 1e-3,
                                   Q = NULL,
                                   simple = NULL,
                                   nquadgh_u = nquadgh,
                                   nquadgh_zeta = nquadgh,
                                   nquadgh_re = nquadgh,
                                   rho_pattern = NULL,
                                   two_tier = NULL) {
  identification <- match.arg(identification)
  threshold_identification <- match.arg(threshold_identification)
  if (missing(growth_shape) && (!is.null(rho_pattern) ||
      (is.numeric(time_scores) && anyNA(time_scores)))) growth_shape <- "rho"
  growth_shape <- match.arg(growth_shape, c("offset", "rho"))
  if (!is.null(rho_pattern) && growth_shape != "rho") {
    stop("rho_pattern requires growth_shape = 'rho'.", call. = FALSE)
  }
  if (length(methInteg) == 1L && methInteg %in% c("RE_AGH_NESTED", "RE_AGH_FULL_EQUIV")) {
    methInteg <- sub("^RE_AGH_", "AGH_RE_", methInteg)
  }
  methInteg <- match.arg(methInteg)
  if (is.null(two_tier)) two_tier <- methInteg != "QMC"
  two_tier <- check_flag(two_tier, "two_tier")
  if (methInteg == "QMC" && two_tier) {
    stop("QMC requires two_tier = FALSE.", call. = FALSE)
  }
  if (!two_tier && methInteg %in% c("AGH_RE_NESTED", "AGH_RE_FULL_EQUIV")) {
    stop("The growth coefficient rules require two_tier = TRUE.", call. = FALSE)
  }
  if (!two_tier && methInteg == "AGH_FULL_EQUIV") methInteg <- "AGH"
  joint_rule <- if (two_tier) 0L else if (methInteg == "AGH") 2L else 1L
  maxiter <- check_count(maxiter, "maxiter", 0L)
  if (!is.null(update_gh)) update_gh <- check_count(update_gh, "update_gh", 0L)
  if (missing(simple_growth) && (!is.null(X) || !is.null(Z))) simple_growth <- FALSE
  if (missing(delta_free) && !isTRUE(simple_growth)) delta_free <- FALSE
  if (is.null(hessian_method)) {
  hessian_method <- if (isTRUE(compute_hessian)) {
    "numeric_fixed"
  } else {
    "none"
  }
}

    hessian_method <- match.arg(
      hessian_method,
      c(
        "none",
        "numeric_fixed",
        "analytic_fixed"
      )
    )

compute_hessian <- hessian_method != "none"
  if (joint_rule != 0L && hessian_method == "analytic_fixed") {
    hessian_method <- "numeric_fixed"
  }

  if (length(analytic_score) != 1L || is.na(analytic_score)) {
    stop("analytic_score deve essere TRUE/FALSE.")
  }

  if (length(tol_loglik) != 1L || !is.finite(tol_loglik) || tol_loglik <= 0) {
    stop("tol_loglik deve essere uno scalare finito > 0.")
  }
  if (length(tol_par) != 1L || !is.finite(tol_par) || tol_par <= 0) {
    stop("tol_par deve essere uno scalare finito > 0.")
  }
  if (length(debug_par_max_eval) != 1L || is.na(debug_par_max_eval) || debug_par_max_eval < 0) {
    stop("debug_par_max_eval deve essere uno scalare >= 0 oppure Inf.")
  }

  is_agh <- grepl("^AGH", methInteg)
  is_growth_re_agh <- methInteg %in% c("AGH_RE_NESTED", "AGH_RE_FULL_EQUIV")

  packed <- prepare_long_data(data, namvar, occ, ID)
  nvar <- length(namvar)
  data_report <- list(
    nrow_input = nrow(data),
    nrow_model = nrow(data),
    nsub = length(packed$subjects),
    occasions = packed$times,
    completed_rows = 0L,
    subject_covariates = unique(c(intercept_covariates, slope_covariates))
  )
  if (isTRUE(verbose)) {
    cat("data check: rows =", nrow(data), "subjects =", data_report$nsub,
        "no structural rows added\n")
  }
  loadmask <- loading_mask(namvar, Q, simple, loadmask)
  Q <- ncol(loadmask)
  if (is.null(factor_names)) factor_names <- colnames(loadmask)
  if (length(factor_names) != Q || anyNA(factor_names) || anyDuplicated(factor_names)) {
    stop("factor_names must give distinct names for all factors.")
  }
  colnames(loadmask) <- factor_names
  if (qr(t(loadmask))$rank < Q) stop("The loading mask must have full column rank.")

  if (identification %in% c("marker", "mixed_marker")) {
    marker_index <- .growth_resolve_marker(marker_item, namvar, loadmask, factor_names)
    marker_names <- setNames(namvar[marker_index], factor_names)
  } else {
    if (!is.null(marker_item)) {
      stop("marker_item \u00e8 utilizzabile solo con identification='marker' o 'mixed_marker'.")
    }
    marker_index <- setNames(integer(Q), factor_names)
    marker_names <- NULL
  }

  data_sorted <- data[order(data[[ID]], data[[occ]]), , drop = FALSE]
  occ_vals <- packed$times
  n <- length(occ_vals)
  if (n < 2L) stop("servono almeno due occasioni.")
  IDs <- packed$subjects
  nsub <- length(IDs)

  if (identification == "baseline_std" && !isTRUE(simple_growth)) {
    stop("baseline_std richiede simple_growth=TRUE, perch\u00e9 il vincolo Var(eta_1|x)=1 \u00e8 definito sul random intercept.")
  }
  rho_fixed_vec <- NULL
  if (!is.null(rho_fixed)) {
    if (growth_shape != "rho") stop("rho_fixed \u00e8 disponibile solo con growth_shape='rho'.")
    rho_fixed_vec <- as.double(rho_fixed)
    if (!length(rho_fixed_vec) %in% c(1L, Q) || any(!is.finite(rho_fixed_vec))) {
      stop("rho_fixed deve essere uno scalare finito oppure un vettore finito di lunghezza Q.")
    }
    if (length(rho_fixed_vec) == 1L) rho_fixed_vec <- rep(rho_fixed_vec, Q)
    names(rho_fixed_vec) <- factor_names
  }

  rho_spec <- NULL
  zero_time <- rep.int(1L, Q)
  if (growth_shape == "rho") {
    if (!isTRUE(simple_growth)) stop("The rho time basis requires simple_growth = TRUE.", call. = FALSE)
    if (is.null(rho_pattern) && is.numeric(time_scores) && anyNA(time_scores)) {
      rho_pattern <- time_scores
    }
    if (is.null(rho_pattern) && n != 3L && is.null(rho_fixed_vec) && is.null(delta_mask)) {
      rho_pattern <- if (is.null(time_scores)) c(0, 1, rep(NA_real_, n - 2L)) else time_scores
    }
    if (!is.null(rho_pattern)) {
      if (!is.null(rho_fixed_vec) || !is.null(delta_mask)) {
        stop("Specify fixed and free scores in rho_pattern when using a time-score pattern.", call. = FALSE)
      }
      rho_spec <- growth_rho_pattern(rho_pattern, occ_vals, factor_names)
      zero_time <- rho_spec$zero_time
      time_scores <- rho_spec$initial[, 1L]
      delta_mask <- rep.int(0L, n)
      delta_free <- FALSE
      if (hessian_method == "analytic_fixed") hessian_method <- "numeric_fixed"
    } else {
      if (is.null(time_scores)) time_scores <- c(0, 1, rep(0, n - 2L))
      if (is.null(delta_mask)) delta_mask <- as.integer(seq_len(n) == n)
      delta_free <- is.null(rho_fixed_vec)
      if (n == 2L && isTRUE(delta_free)) {
        stop("Two occasions supply the anchors 0 and 1. Use a fixed rho pattern.", call. = FALSE)
      }
    }
  }

  intercept_mask <- .growth_normalize_intercept_mask(intercept_time_free, namvar, occ_vals)
  gamma_dep <- .growth_make_gamma_dep(loadmask, intercept_mask, identification)

  if (is.null(random_names) && isTRUE(simple_growth) &&
      identification == "effects" && growth_shape == "offset") {
    random_names <- rep(list(c("u0", "u1")), Q)
  }

  design <- .growth_build_design(
    data_sorted = data_sorted, ID = ID, occ = occ, IDs = IDs, occ_vals = occ_vals, Q = Q,
    simple_growth = simple_growth, same_design = same_design,
    subject_specific_design = subject_specific_design,
    time_scores = time_scores, fixed_covariates = fixed_covariates,
    time_interactions = time_interactions,
    intercept_covariates = intercept_covariates,
    slope_covariates = slope_covariates,
    X = X, Z = Z, delta_mask = delta_mask,
    fixed_names = fixed_names, random_names = random_names
  )
  if (!isTRUE(simple_growth) && isTRUE(delta_free)) stop("delta_free richiede simple_growth=TRUE.")
  if (any(!design$delta_mask %in% c(0, 1))) stop("delta_mask deve contenere 0/1.")
  if (!is.null(rho_spec)) {
    basis <- growth_rho_update_design(design$X, design$Z, rho_spec$initial)
    design$X <- basis$X
    design$Z <- basis$Z
    rho_spec$dim_x <- dim(design$X)
    rho_spec$dim_z <- dim(design$Z)
  }
  p_fix <- as.integer(design$p_fix)
  p_ran <- as.integer(design$p_ran)
  x_slope_col <- as.integer(grepl("^S:", design$fixed_names[[1L]]))
  if (length(x_slope_col) != p_fix) stop("mask delle covariate slope incompatibile con p_fix.")
  if (Q > 1L) {
    ok_mask <- vapply(design$fixed_names, function(z) identical(as.integer(grepl("^S:", z)), x_slope_col), logical(1L))
    if (!all(ok_mask)) stop("la posizione delle covariate slope deve essere la stessa per tutti i fattori.")
  }
  if (!is.null(rho_spec) && length(design$interaction_columns)) {
    columns <- design$interaction_columns
    for (q in seq_len(Q)) {
      design$X[, , q, columns] <- design$X[, , q, 2L + design$time_interactions[[q]]]
    }
    x_slope_col[columns] <- 1L
  }
  dim_eta <- n * Q
  dim_joint <- dim_eta + if (isTRUE(re_item)) nvar else 0L
  dim_u <- Q * p_ran
  dim_zeta <- n * Q
  use_analytic_score_opt <- isTRUE(analytic_score) &&
    methInteg %in% c("AGH_FULL_EQUIV", "AGH_RE_FULL_EQUIV") &&
    growth_shape == "rho" && is.null(rho_spec) && n == 3L &&
    identification != "effects" && Q == 1L && p_ran == 2L &&
    maxiter > 0L && !is.null(update_gh) && update_gh > 0L

  c0 <- as.integer(packed$c)
  cm1 <- c0 - 1L
  cm3 <- pmax(cm1 - 2L, 0L)
  ntre0 <- sum(cm3)
  Y0 <- as.integer(packed$y)
  nmes0 <- packed$nmes
  nitem0 <- as.integer(packed$time)
  nobs <- length(Y0)

  empty_rule <- list(nquad = 1L, nodes = 0, weights = 1)
  if (is_growth_re_agh) {
    u_rule <- make_gh_tensor(dim_u, nquadgh_u)
    zeta_rule <- make_gh_tensor(dim_zeta, nquadgh_zeta)
    eta_rule <- zeta_rule
  } else {
    rule_dimension <- if (two_tier) dim_eta else dim_joint
    eta_rule <- if (methInteg == "QMC") make_qmc(rule_dimension, nquadqmc) else make_gh_tensor(rule_dimension, nquadgh)
    u_rule <- list(nquad = 1L, nodes = rep(0, dim_u), weights = 1)
    zeta_rule <- list(nquad = 1L, nodes = rep(0, dim_zeta), weights = 1)
  }
  re_rule <- if (two_tier && isTRUE(re_item)) make_gh_1d(nquadgh_re) else empty_rule

  ref_item <- vapply(seq_len(Q), function(q) tail(which(loadmask[, q] == 1L), 1L), integer(1L))
  dep0 <- gamma_dep[, 1L]
  free_gamma0 <- setdiff(seq_len(nvar), dep0)
  free_marker_gamma <- setdiff(seq_len(nvar), unique(marker_index))

  nm <- character(0)
  for (k in seq_len(nvar)) {
    idx <- .growth_free_threshold_idx(cm1[k], threshold_identification)
    if (length(idx)) nm <- c(nm, sprintf("tau_free[%s]_c%d", namvar[k], idx))
  }
  if (identification == "effects") {
    for (q in seq_len(Q)) for (k in which(loadmask[, q] == 1L)) if (k != ref_item[q]) {
      nm <- c(nm, sprintf("lambda_free[%s,%s]", namvar[k], factor_names[q]))
    }
  } else if (identification %in% c("marker", "mixed_marker")) {
    for (q in seq_len(Q)) for (k in which(loadmask[, q] == 1L)) if (k != marker_index[q]) {
      nm <- c(nm, sprintf("lambda[%s,%s]", namvar[k], factor_names[q]))
    }
  } else {
    for (q in seq_len(Q)) for (k in which(loadmask[, q] == 1L)) {
      nm <- c(nm, sprintf("lambda[%s,%s]", namvar[k], factor_names[q]))
    }
  }
  if (identification == "effects") {
    nm <- c(nm, sprintf("nu_free[%s]", namvar[free_gamma0]))
    if (n > 1L) for (tt in 2:n) {
      vv <- which(intercept_mask[, tt]); dd <- gamma_dep[, tt]
      ff <- setdiff(vv, dd)
      if (length(ff)) nm <- c(nm, sprintf("nu_dev[%s,t%d]", namvar[ff], tt))
    }
  } else if (identification == "marker") {
    nm <- c(nm, sprintf("nu[%s,t1]", namvar[free_marker_gamma]))
    if (n > 1L) for (tt in 2:n) for (k in which(intercept_mask[, tt])) {
      nm <- c(nm, sprintf("nu[%s,t%d]", namvar[k], tt))
    }
  } else {
    nm <- c(nm, sprintf("nu[%s,t1]", namvar))
    if (n > 1L) for (tt in 2:n) for (k in which(intercept_mask[, tt])) {
      nm <- c(nm, sprintf("nu[%s,t%d]", namvar[k], tt))
    }
  }
  if (!isTRUE(theta_time)) {
    for (k in seq_len(nvar)) if (cm1[k] > 1L) nm <- c(nm, sprintf("log_theta[%s]", namvar[k]))
  } else {
    for (k in seq_len(nvar)) if (cm1[k] > 1L) for (tt in seq_len(n)) nm <- c(nm, sprintf("log_theta[%s,t%d]", namvar[k], tt))
  }
  for (r in seq_len(p_fix)) for (q in seq_len(Q)) {
    if (identification %in% c("baseline_std", "mixed_marker") && r == 1L) next
    fn <- design$fixed_names[[q]][r]
    if (r == 1L) {
      nm <- c(nm, sprintf("alpha0[%s]", factor_names[q]))
    } else if (r == 2L) {
      nm <- c(nm, sprintf("alpha1[%s]", factor_names[q]))
    } else {
      nm <- c(nm, sprintf("beta[%s,%s]", fn, factor_names[q]))
    }
  }
  if (!is.null(rho_spec)) {
    free_names <- rho_spec$names[rho_spec$free]
    rho_spec$index <- length(nm) + seq_along(free_names)
    nm <- c(nm, free_names)
  }
  if (isTRUE(delta_free)) nm <- c(nm, sprintf("%s[%s]", if (growth_shape == "rho") "rho" else "delta", factor_names))

  u_names <- character(dim_u)
  for (r in seq_len(p_ran)) for (q in seq_len(Q)) u_names[(r-1L)*Q+q] <- paste0(design$random_names[[q]][r], "_", factor_names[q])
  if (identification %in% c("effects", "marker", "mixed_marker")) {
    for (col in seq_len(dim_u)) for (row in col:dim_u) {
      nm <- c(nm, if (row == col) sprintf("log_Lu[%s,%s]", u_names[row], u_names[col]) else sprintf("Lu[%s,%s]", u_names[row], u_names[col]))
    }
  } else {
    for (d in seq_len(dim_u)) nm <- c(nm, if (d <= Q) sprintf("logit_Var_u[%s]", u_names[d]) else sprintf("log_SD_u[%s]", u_names[d]))
    if (dim_u > 1L) for (row in 2:dim_u) for (col in seq_len(row-1L)) nm <- c(nm, sprintf("atanh_pcor_u[%s,%s]", u_names[row], u_names[col]))
  }
  zeta_free <- matrix(TRUE, n, Q)
  if (identification == "baseline_std") zeta_free[cbind(zero_time, seq_len(Q))] <- FALSE
  for (tt in seq_len(n)) for (q in which(zeta_free[tt, ])) {
    nm <- c(nm, sprintf("log_zeta_var[t%d,%s]", tt, factor_names[q]))
  }
  if (isTRUE(re_item)) nm <- c(nm, sprintf("log_sigmak[%s]", namvar))
  nm_opt <- nm
  npm0 <- length(nm_opt)

  make_start <- function() {
    b <- numeric(npm0); names(b) <- nm_opt; pos <- 0L
    for (k in seq_len(nvar)) {
      idx <- .growth_free_threshold_idx(cm1[k], threshold_identification)
      if (!length(idx)) next
      if (threshold_identification == "central_01") {
        j0 <- (cm1[k] + 1L) %/% 2L
        for (jj in idx) {
          pos <- pos + 1L
          if (isTRUE(monotone_tau)) b[pos] <- 0 else b[pos] <- if (jj < j0) -(j0-jj) else 1 + (jj-(j0+1L)) + 1
        }
      } else {
        nf <- length(idx)
        thr_abs <- seq(1.6, 3.2, length.out = nf)
        if (nf > 0L) {
          if (isTRUE(monotone_tau)) {
            inc <- diff(c(1, thr_abs))
            inc[inc <= 0] <- 1e-6
            vals <- log(inc)
          } else {
            vals <- thr_abs
          }
          for (h in seq_len(nf)) {
            pos <- pos + 1L
            b[pos] <- vals[h]
          }
        }
      }
    }
    nload <- if (identification == "baseline_std") sum(loadmask) else sum(colSums(loadmask)-1L)
    if (nload) { b[pos + seq_len(nload)] <- 1; pos <- pos + nload }
    ng <- if (identification == "effects") {
      length(free_gamma0) + sum(vapply(2:n, function(tt) as.integer(max(0, sum(intercept_mask[,tt])-Q)), integer(1L)))
    } else if (identification == "marker") {
      length(free_marker_gamma) + sum(intercept_mask[, -1L, drop=FALSE])
    } else nvar + sum(intercept_mask[, -1L, drop=FALSE])
    if (ng) { b[pos + seq_len(ng)] <- 0; pos <- pos + ng }
    nth <- if (!isTRUE(theta_time)) sum(cm1>1L) else sum(cm1>1L)*n
    if (nth) { b[pos + seq_len(nth)] <- 0; pos <- pos + nth }
    for (r in seq_len(p_fix)) for (q in seq_len(Q)) {
      if (identification %in% c("baseline_std", "mixed_marker") && r == 1L) next
      pos <- pos + 1L
      b[pos] <- if (r == 2L) 0.5 else 0
    }
    if (!is.null(rho_spec)) {
      b[rho_spec$index] <- rho_spec$initial[rho_spec$free]
      pos <- pos + length(rho_spec$index)
    }
    if (isTRUE(delta_free)) {
      for (q in seq_len(Q)) { pos <- pos + 1L; b[pos] <- if (growth_shape == "rho") 1.5 else 0 }
    }
    if (identification %in% c("effects", "marker", "mixed_marker")) {
      for (col in seq_len(dim_u)) for (row in col:dim_u) {
        pos <- pos + 1L
        b[pos] <- if (row == col) log(if (((row-1L)%/%Q)+1L == 1L) 1 else 0.5) else 0
      }
    } else {
      for (d in seq_len(dim_u)) { pos <- pos + 1L; b[pos] <- if (d <= Q) qlogis(0.6) else log(0.5) }
      if (dim_u > 1L) for (row in 2:dim_u) for (col in seq_len(row-1L)) { pos <- pos + 1L; b[pos] <- 0 }
    }
    n_zeta <- sum(zeta_free)
    b[pos + seq_len(n_zeta)] <- log(if (identification == "baseline_std") 0.4 else 0.5)
    pos <- pos + n_zeta
    if (isTRUE(re_item)) for (k in seq_len(nvar)) {
      pos <- pos + 1L
      b[pos] <- log(if (identification %in% c("effects", "marker", "mixed_marker")) 0.5 else 0.3)
    }
    if (pos != npm0) stop("errore interno nel conteggio dei parametri.")
    b
  }
  b0 <- if (is.null(start)) make_start() else {
    x <- as.double(start)
    if (length(x) != npm0) stop("start length != npm0.")
    if (!is.null(names(start)) && all(nm_opt %in% names(start))) x <- start[nm_opt]
    setNames(as.double(x), nm_opt)
  }

  make_parscale <- function(bref) {
    ps <- rep(NA_real_, length(bref))
    names(ps) <- nm_opt

    clamp_scale <- function(x, smin, smax) {
      pmin(pmax(abs(x), smin), smax)
    }
    set_scale <- function(pattern, smin, smax) {
      idx <- grepl(pattern, nm_opt)
      if (any(idx)) ps[idx] <<- clamp_scale(bref[idx], smin, smax)
    }

    set_scale("^tau_free\\[", 0.5, 3)
    set_scale("^lambda(_free)?\\[", 0.3, 2)
    set_scale("^nu", 0.5, 3)
    set_scale("^log_theta\\[", 0.3, 2)
    set_scale("^alpha0\\[", 0.5, 3)
    set_scale("^alpha1\\[", 0.3, 2)
    set_scale("^beta\\[", 0.3, 3)
    set_scale("^delta\\[", 0.25, 2)
    set_scale("^rho\\[", 0.25, 2)
    set_scale("^log_Lu\\[", 0.3, 2)
    set_scale("^Lu\\[", 0.2, 2)
    set_scale("^logit_Var_u\\[", 0.3, 3)
    set_scale("^log_SD_u\\[", 0.3, 2)
    set_scale("^atanh_pcor_u\\[", 0.2, 2)
    set_scale("^log_zeta_var\\[", 0.3, 2)
    set_scale("^log_sigmak\\[", 0.3, 2)

    idx_na <- is.na(ps)
    if (any(idx_na)) ps[idx_na] <- clamp_scale(bref[idx_na], 0.25, 3)
    if (any(!is.finite(ps)) || any(ps <= 0)) stop("parscale generato non valido.")
    ps
  }

  ps0 <- if (is.null(parscale)) make_parscale(b0) else {
    x <- as.double(parscale)
    if (length(x) != npm0) stop("parscale length != npm0.")
    if (any(!is.finite(x)) || any(x <= 0)) stop("parscale non valido.")
    x
  }

  decode <- function(par) {
    par <- as.double(par); names(par) <- nm_opt; pos <- 0L
    td <- .growth_decode_thresholds(par,pos,cm1,namvar,threshold_identification,monotone_tau)
    tau <- td$tau; pos <- td$pos
    lambda <- matrix(0,nvar,Q,dimnames=list(namvar,factor_names))
    if (identification == "effects") {
      for (q in seq_len(Q)) {
        lambda[ref_item[q],q] <- sum(loadmask[,q])
        for (k in which(loadmask[,q]==1L)) if (k != ref_item[q]) {
          pos <- pos+1L; lambda[k,q] <- par[pos]; lambda[ref_item[q],q] <- lambda[ref_item[q],q]-lambda[k,q]
        }
      }
    } else if (identification %in% c("marker", "mixed_marker")) {
      for (q in seq_len(Q)) {
        lambda[marker_index[q], q] <- 1
        for (k in which(loadmask[, q] == 1L)) if (k != marker_index[q]) {
          pos <- pos + 1L
          lambda[k, q] <- par[pos]
        }
      }
    } else {
      for (q in seq_len(Q)) for (k in which(loadmask[,q]==1L)) {
        pos <- pos+1L
        lambda[k,q]<-par[pos]
      }
    }

    gamma <- numeric(nvar); names(gamma) <- namvar
    gamma_time <- matrix(0,nvar,n,dimnames=list(namvar,paste0("t",seq_len(n))))
    if (identification == "effects") {
      for (k in free_gamma0) { pos<-pos+1L; gamma[k]<-par[pos] }
      A <- t(loadmask[dep0,,drop=FALSE])
      rhs <- -as.vector(t(loadmask[free_gamma0,,drop=FALSE]) %*% gamma[free_gamma0])
      gamma[dep0] <- as.vector(solve(A,rhs))
      gamma_time[,] <- gamma
      if (n>1L) for (tt in 2:n) {
        vv <- which(intercept_mask[,tt]); if(!length(vv)) next
        dd <- gamma_dep[,tt]; ff <- setdiff(vv,dd); dev <- numeric(nvar)
        for(k in ff){pos<-pos+1L;dev[k]<-par[pos]}
        rhs <- if(length(ff)) -as.vector(t(loadmask[ff,,drop=FALSE])%*%dev[ff]) else rep(0,Q)
        dev[dd] <- as.vector(solve(t(loadmask[dd,,drop=FALSE]),rhs))
        gamma_time[,tt] <- gamma+dev
      }
    } else if (identification == "marker") {
      for(k in free_marker_gamma){pos<-pos+1L;gamma[k]<-par[pos]}
      gamma_time[,] <- gamma
      if(n>1L) for(tt in 2:n) for(k in which(intercept_mask[,tt])){pos<-pos+1L;gamma_time[k,tt]<-par[pos]}
    } else {
      for(k in seq_len(nvar)){pos<-pos+1L;gamma[k]<-par[pos]}
      gamma_time[,] <- gamma
      if(n>1L) for(tt in 2:n) for(k in which(intercept_mask[,tt])){pos<-pos+1L;gamma_time[k,tt]<-par[pos]}
    }

    theta <- matrix(1,nvar,n,dimnames=list(namvar,paste0("t",seq_len(n))))
    if (!isTRUE(theta_time)) {
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

    beta <- matrix(0,Q,p_fix,dimnames=list(factor_names,design$fixed_names[[1L]]))
    for(r in seq_len(p_fix)) for(q in seq_len(Q)) {
      if(identification %in% c("baseline_std", "mixed_marker") && r==1L) next
      pos<-pos+1L;beta[q,r]<-par[pos]
    }
    delta <- setNames(rep(0,Q),factor_names)
    if (growth_shape == "rho" && !is.null(rho_fixed_vec)) delta[] <- rho_fixed_vec
    if(isTRUE(delta_free)) for(q in seq_len(Q)){pos<-pos+1L;delta[q]<-par[pos]}

    rho <- NULL
    if (!is.null(rho_spec)) {
      rho <- growth_rho_scores(rho_spec, par)
      pos <- pos + length(rho_spec$index)
    } else if (growth_shape == "rho") {
      rho <- delta
    }

    cv <- .growth_decode_covariance(par, pos, Q, p_ran, n, u_names, identification, zero_time)
    pos <- cv$pos
    sigmak <- NULL
    if(isTRUE(re_item)){sigmak<-setNames(numeric(nvar),namvar);for(k in seq_len(nvar)){pos<-pos+1L;sigmak[k]<-exp(par[pos])}}
    if (pos != npm0) stop("decode: pos = ", pos, " ma npm0 = ", npm0, ".")

    Xe <- design$X; Ze <- design$Z
    if (!is.null(rho_spec)) {
      basis <- growth_rho_update_design(Xe, Ze, rho)
      Xe <- basis$X
      Ze <- basis$Z
    }
    if(dim(Xe)[1L]==1L && nsub>1L) Xe <- Xe[rep(1L,nsub),,,,drop=FALSE]
    if(dim(Ze)[1L]==1L && nsub>1L) Ze <- Ze[rep(1L,nsub),,,,drop=FALSE]
    dm <- design$delta_mask
    if(dim(dm)[1L]==1L && nsub>1L) dm <- dm[rep(1L,nsub),,,drop=FALSE]
    delta_active <- isTRUE(delta_free) || (growth_shape == "rho" && !is.null(rho_fixed_vec))
    if(delta_active) for(q in seq_len(Q)) {
      hit <- dm[,,q] != 0
      if(growth_shape=="rho") {
        tmp <- Xe[,,q,design$delta_x_col]; tmp[hit] <- delta[q]; Xe[,,q,design$delta_x_col] <- tmp
        tmp <- Ze[,,q,design$delta_z_col]; tmp[hit] <- delta[q]; Ze[,,q,design$delta_z_col] <- tmp
      } else {
        Xe[,,q,design$delta_x_col] <- Xe[,,q,design$delta_x_col] + delta[q]*dm[,,q]
        Ze[,,q,design$delta_z_col] <- Ze[,,q,design$delta_z_col] + delta[q]*dm[,,q]
      }
    }
    if (any(x_slope_col == 1L)) {
      for (q in seq_len(Q)) {
        score_eff <- Xe[, , q, design$delta_x_col]
        for (r in which(x_slope_col == 1L)) Xe[, , q, r] <- Xe[, , q, r] * score_eff
      }
    }
    kappa <- array(0,c(nsub,n,Q),dimnames=list(as.character(IDs),paste0("t",seq_len(n)),factor_names))
    for(q in seq_len(Q)) for(r in seq_len(p_fix)) kappa[,,q] <- kappa[,,q] + Xe[,,q,r]*beta[q,r]

    score <- if (isTRUE(simple_growth) && p_fix >= 2L) {
      array(Xe[, , , design$delta_x_col, drop=FALSE], dim=c(nsub,n,Q),
            dimnames=list(as.character(IDs),paste0("t",seq_len(n)),factor_names))
    } else NULL
    list(tau=tau,lambda=lambda,gamma=gamma,gamma_time=gamma_time,theta=theta,beta=beta,
         alpha0=if(p_fix>=1L)beta[,1L] else NULL,
         alpha1=if(p_fix>=2L)beta[,2L] else NULL,
         muS=if(growth_shape=="rho" && p_fix>=2L)beta[,2L] else NULL,
         delta=if(growth_shape=="rho")NULL else delta,rho=rho,
         score=score,kappa=kappa,X_effective=Xe,Z_effective=Ze,Lu=cv$Lu,Sigma_u=cv$Sigma_u,
         zeta_var=cv$zeta_var,sigmak=sigmak,
         baseline_var=diag(cv$Sigma_u)[seq_len(Q)] + cv$zeta_var[cbind(zero_time, seq_len(Q))])
  }

  invisible(decode(b0))

  re_item0 <- as.integer(isTRUE(re_item))
  theta_time0 <- as.integer(isTRUE(theta_time))
  delta_free0 <- as.integer(isTRUE(delta_free))
  monotone_tau0 <- as.integer(isTRUE(monotone_tau))

  .growth_set_fortran_options(identification, threshold_identification, growth_shape,
                              intercept_mask, gamma_dep, x_slope_col, nvar, n, Q, p_fix,
                              rho_fixed = rho_fixed_vec, marker_index = marker_index, zero_time = zero_time)

  Y0_c <- as.integer(Y0); nvar0_c <- as.integer(nvar); nsub0_c <- as.integer(nsub); nobs0_c <- as.integer(nobs)
  n0_c <- as.integer(n); nmes0_c <- as.integer(nmes0); c0_c <- as.integer(c0); cm10_c <- as.integer(cm1)
  ntre0_c <- as.integer(ntre0); nitem0_c <- as.integer(nitem0); npm0_c <- as.integer(npm0)
  Q0_c <- as.integer(Q); loadmask0_c <- as.integer(loadmask)
  fit_id0_c <- next_fit_id()
  p_fix0_c <- as.integer(p_fix); p_ran0_c <- as.integer(p_ran)
  X0_c <- as.double(design$X); Z0_c <- as.double(design$Z); delta_mask0_c <- as.double(design$delta_mask)
  dx_c <- as.integer(design$delta_x_col); dz_c <- as.integer(design$delta_z_col)

  nquad_eta0_c <- as.integer(eta_rule$nquad); nodes_eta0_c <- as.double(eta_rule$nodes); weights_eta0_c <- as.double(eta_rule$weights)
  nquad_u0_c <- as.integer(u_rule$nquad); nodes_u0_c <- as.double(u_rule$nodes); weights_u0_c <- as.double(u_rule$weights)
  nquad_zeta0_c <- as.integer(zeta_rule$nquad); nodes_zeta0_c <- as.double(zeta_rule$nodes); weights_zeta0_c <- as.double(zeta_rule$weights)
  nquad_re0_c <- as.integer(re_rule$nquad); nodes_re0_c <- as.double(re_rule$nodes); weights_re0_c <- as.double(re_rule$weights)
  joint_grid <- NULL

  loglik_call <- function(par, update_grid = 1L) {
    if (!is.null(rho_spec)) {
      native <- growth_rho_native(par, rho_spec, X0_c, Z0_c)
      par <- native$par
      X0_c <- native$X
      Z0_c <- native$Z
      npm0_c <- length(par)
    }
    .growth_set_fortran_options(identification, threshold_identification, growth_shape,
                                intercept_mask, gamma_dep, x_slope_col, nvar, n, Q, p_fix,
                                rho_fixed = rho_fixed_vec, marker_index = marker_index, zero_time = zero_time)
    if (joint_rule != 0L) {
      contributions <- loglik_growth_tt_agh_full_equiv_i_vec(
        par, Y0_c, nvar0_c, nsub0_c, nobs0_c, n0_c, nmes0_c, c0_c, cm10_c, ntre0_c,
        nquad_eta0_c, nodes_eta0_c, weights_eta0_c, nquad_re0_c, nodes_re0_c, weights_re0_c,
        nitem0_c, npm0_c, re_item0, theta_time0, delta_free0, monotone_tau0,
        as.integer(update_grid), fit_id0_c, Q0_c, loadmask0_c, p_fix0_c, p_ran0_c,
        X0_c, Z0_c, delta_mask0_c, dx_c, dz_c, design$subject_specific_design, TRUE,
        joint_rule0 = joint_rule, joint_grid0 = joint_grid
      )
      if (any(!is.finite(contributions) | contributions <= -1e8)) return(1e9)
      if (joint_rule == 2L && update_grid == 1L) joint_grid <<- attr(contributions, "joint_grid")
      -sum(contributions)
    } else if (methInteg == "AGH") {
      loglik_growth_tt_agh(par,Y0_c,nvar0_c,nsub0_c,nobs0_c,n0_c,nmes0_c,c0_c,cm10_c,ntre0_c,
                           nquad_eta0_c,nodes_eta0_c,weights_eta0_c,nquad_re0_c,nodes_re0_c,weights_re0_c,
                           nitem0_c,npm0_c,re_item0,theta_time0,delta_free0,monotone_tau0,
                           as.integer(update_grid),fit_id0_c,Q0_c,loadmask0_c,p_fix0_c,p_ran0_c,X0_c,Z0_c,delta_mask0_c,dx_c,dz_c,
                           design$subject_specific_design,TRUE)
    } else if (methInteg == "AGH_FULL_EQUIV") {
      loglik_growth_tt_agh_full_equiv(par,Y0_c,nvar0_c,nsub0_c,nobs0_c,n0_c,nmes0_c,c0_c,cm10_c,ntre0_c,
                                      nquad_eta0_c,nodes_eta0_c,weights_eta0_c,nquad_re0_c,nodes_re0_c,weights_re0_c,
                                      nitem0_c,npm0_c,re_item0,theta_time0,delta_free0,monotone_tau0,
                                      as.integer(update_grid),fit_id0_c,Q0_c,loadmask0_c,p_fix0_c,p_ran0_c,X0_c,Z0_c,delta_mask0_c,dx_c,dz_c,
                                      design$subject_specific_design,TRUE)
    } else if (methInteg == "AGH_RE_NESTED") {
      loglik_growth_re_tt_agh_nested(par,Y0_c,nvar0_c,nsub0_c,nobs0_c,n0_c,nmes0_c,c0_c,cm10_c,ntre0_c,
                                     nquad_u0_c,nodes_u0_c,weights_u0_c,nquad_zeta0_c,nodes_zeta0_c,weights_zeta0_c,
                                     nquad_re0_c,nodes_re0_c,weights_re0_c,nitem0_c,npm0_c,re_item0,theta_time0,delta_free0,monotone_tau0,
                                     as.integer(update_grid),fit_id0_c,Q0_c,loadmask0_c,p_fix0_c,p_ran0_c,X0_c,Z0_c,delta_mask0_c,dx_c,dz_c,
                                     design$subject_specific_design,TRUE)
    } else if (methInteg == "AGH_RE_FULL_EQUIV") {
      loglik_growth_re_tt_agh_full_equiv(par,Y0_c,nvar0_c,nsub0_c,nobs0_c,n0_c,nmes0_c,c0_c,cm10_c,ntre0_c,
                                         nquad_u0_c,nodes_u0_c,weights_u0_c,nquad_zeta0_c,nodes_zeta0_c,weights_zeta0_c,
                                         nquad_re0_c,nodes_re0_c,weights_re0_c,nitem0_c,npm0_c,re_item0,theta_time0,delta_free0,monotone_tau0,
                                         as.integer(update_grid),fit_id0_c,Q0_c,loadmask0_c,p_fix0_c,p_ran0_c,X0_c,Z0_c,delta_mask0_c,dx_c,dz_c,
                                         design$subject_specific_design,TRUE)
    } else {
      loglik_growth_tt(par,Y0_c,nvar0_c,nsub0_c,nobs0_c,n0_c,nmes0_c,c0_c,cm10_c,ntre0_c,
                       nquad_eta0_c,nodes_eta0_c,weights_eta0_c,nquad_re0_c,nodes_re0_c,weights_re0_c,
                       nitem0_c,npm0_c,re_item0,theta_time0,delta_free0,monotone_tau0,Q0_c,loadmask0_c,
                       p_fix0_c,p_ran0_c,X0_c,Z0_c,delta_mask0_c,dx_c,dz_c,design$subject_specific_design,TRUE)
    }
  }

  score_opt_symbol <- NULL
  score_opt_cache <- new.env(parent = emptyenv())
  score_opt_cache$par <- NULL
  score_opt_cache$nll <- NULL
  score_opt_cache$score <- NULL
  score_opt_cache$ok <- NULL

  if (isTRUE(use_analytic_score_opt)) {
    score_opt_name <- if (identical(methInteg, "AGH_FULL_EQUIV")) {
      "loglik_growth_tt_agh_full_equiv_score_q1"
    } else {
      "loglik_growth_re_tt_agh_full_equiv_score_q1"
    }
    if (!native_has(score_opt_name)) {
      stop("routine ", score_opt_name, " non caricata.")
    }
    score_opt_symbol <- native_symbol(score_opt_name)
  }

  reset_score_opt_cache <- function() {
    score_opt_cache$par <- NULL
    score_opt_cache$nll <- NULL
    score_opt_cache$score <- NULL
    score_opt_cache$ok <- NULL
    invisible(NULL)
  }

  score_opt_eval_raw <- function(par) {
    .growth_set_fortran_options(
      identification, threshold_identification, growth_shape,
      intercept_mask, gamma_dep, x_slope_col, nvar, n, Q, p_fix,
      rho_fixed = rho_fixed_vec, marker_index = marker_index
    )

    if (identical(methInteg, "AGH_FULL_EQUIV")) {
      z <- native_call(
        score_opt_symbol,
        b0 = as.double(par),
        Y0 = Y0_c,
        nvar0 = nvar0_c,
        nsub0 = nsub0_c,
        nobs0 = nobs0_c,
        n0 = n0_c,
        nmes0 = nmes0_c,
        c0 = c0_c,
        cm10 = cm10_c,
        ntre0 = ntre0_c,
        nquad_eta0 = nquad_eta0_c,
        nodes_eta0 = nodes_eta0_c,
        weights_eta0 = weights_eta0_c,
        nquad_re0 = nquad_re0_c,
        nodes_re0 = nodes_re0_c,
        weights_re0 = weights_re0_c,
        nitem0 = nitem0_c,
        npm0 = npm0_c,
        re_item0 = re_item0,
        theta_time0 = theta_time0,
        delta_free0 = delta_free0,
        monotone_tau0 = monotone_tau0,
        update_gh0 = 0L,
        fit_id0 = fit_id0_c,
        Q0 = Q0_c,
        loadmask0 = loadmask0_c,
        p_fix0 = as.integer(if (isTRUE(design$subject_specific_design)) p_fix else -p_fix),
        p_ran0 = p_ran0_c,
        X0 = X0_c,
        Z0 = Z0_c,
        delta_mask0 = delta_mask0_c,
        delta_x_col0 = dx_c,
        delta_z_col0 = dz_c,
        vrais = double(1L),
        score0 = double(npm0),
        ok0 = integer(1L)
      )
    } else {
      z <- native_call(
        score_opt_symbol,
        b0 = as.double(par),
        Y0 = Y0_c,
        nvar0 = nvar0_c,
        nsub0 = nsub0_c,
        nobs0 = nobs0_c,
        n0 = n0_c,
        nmes0 = nmes0_c,
        c0 = c0_c,
        cm10 = cm10_c,
        ntre0 = ntre0_c,
        nquad_u0 = nquad_u0_c,
        nodes_u0 = nodes_u0_c,
        weights_u0 = weights_u0_c,
        nquad_zeta0 = nquad_zeta0_c,
        nodes_zeta0 = nodes_zeta0_c,
        weights_zeta0 = weights_zeta0_c,
        nquad_re0 = nquad_re0_c,
        nodes_re0 = nodes_re0_c,
        weights_re0 = weights_re0_c,
        nitem0 = nitem0_c,
        npm0 = npm0_c,
        re_item0 = re_item0,
        theta_time0 = theta_time0,
        delta_free0 = delta_free0,
        monotone_tau0 = monotone_tau0,
        update_gh0 = 0L,
        fit_id0 = fit_id0_c,
        Q0 = Q0_c,
        loadmask0 = loadmask0_c,
        p_fix0 = as.integer(if (isTRUE(design$subject_specific_design)) p_fix else -p_fix),
        p_ran0 = p_ran0_c,
        X0 = X0_c,
        Z0 = Z0_c,
        delta_mask0 = delta_mask0_c,
        delta_x_col0 = dx_c,
        delta_z_col0 = dz_c,
        vrais = double(1L),
        score0 = double(npm0),
        ok0 = integer(1L)
      )
    }

    if (length(z$ok0) != 1L || z$ok0 != 1L ||
        !is.finite(z$vrais) || z$vrais >= 1e8 ||
        length(z$score0) != npm0 || any(!is.finite(z$score0))) {
      return(list(
        nll = 1e9,
        score = rep(0, npm0),
        ok = FALSE
      ))
    }

    list(
      nll = as.numeric(z$vrais),
      score = as.numeric(z$score0),
      ok = TRUE
    )
  }

  score_opt_eval <- function(par) {
    p <- as.double(par)
    if (!is.null(score_opt_cache$par) && identical(p, score_opt_cache$par)) {
      return(list(
        nll = score_opt_cache$nll,
        score = score_opt_cache$score,
        ok = score_opt_cache$ok
      ))
    }

    debug_eval_counter <<- debug_eval_counter + 1L
    if (isTRUE(debug_par) && isTRUE(debug_par_each_eval) &&
        debug_eval_counter <= debug_par_max_eval) {
      debug_print_growth_par(par, paste0("analytic score eval ", debug_eval_counter))
    }

    t0_score <- proc.time()[["elapsed"]]
    z <- score_opt_eval_raw(par)

    if (isTRUE(profile)) {
      nfn <<- nfn + 1L
      fn_elapsed <<- c(fn_elapsed, proc.time()[["elapsed"]] - t0_score)
    }

    if (!isTRUE(z$ok) && isTRUE(debug_par)) {
      cat("BAD analytic-score objective at eval", debug_eval_counter, "value =", z$nll, "\n")
    }

    score_opt_cache$par <- p
    score_opt_cache$nll <- z$nll
    score_opt_cache$score <- z$score
    score_opt_cache$ok <- z$ok
    z
  }

  fn_fixed_analytic <- function(par) score_opt_eval(par)$nll
  gr_fixed_analytic <- function(par) score_opt_eval(par)$score

  debug_eval_counter <- 0L
  nfn <- 0L
  fn_elapsed <- numeric(0)

  debug_print_growth_par <- function(par, label = "") {
    p <- as.double(par)
    names(p) <- nm_opt
    cat("\n================ debug growth par:", label, "================\n")
    cat("length:", length(p), "npm0:", npm0, "\n")
    cat("n:", n, "Q:", Q, "nvar:", nvar, "dim_eta:", dim_eta, "\n")
    cat("identification:", identification, "threshold identification:", threshold_identification, "growth shape:", growth_shape, "\n")
    if (any(!is.finite(p))) cat("BAD par non finite:", which(!is.finite(p)), "\n")
    if (isTRUE(debug_par_full)) {
      cat("raw par named:\n")
      print(p)
    } else {
      grp <- list(
        thresholds = grepl("^tau_", nm_opt),
        loadings = grepl("^lambda", nm_opt),
        intercepts = grepl("^nu", nm_opt),
        residuals = grepl("^log_theta", nm_opt),
        growth = grepl("^(alpha0|alpha1|beta\\[|rho|delta)", nm_opt),
        random_growth = grepl("^(log_Lu|Lu\\[|logit_Var_u|log_SD_u|atanh_pcor_u)", nm_opt),
        zeta = grepl("^log_zeta_var", nm_opt),
        item_re = grepl("^log_sigmak", nm_opt)
      )
      for (g in names(grp)) {
        ii <- which(grp[[g]])
        if (length(ii)) {
          rr <- range(p[ii], finite = TRUE)
          cat(g, ": n =", length(ii), "range =", paste(rr, collapse = " "), "\n")
        }
      }
    }
    cat("============================================================\n\n")
    flush.console()
  }

  if (isTRUE(debug_par)) debug_print_growth_par(b0, "initial b0 before optim")

  timed_fn <- function(update_gh_eval) {
    force(update_gh_eval)
    function(par) {
      debug_eval_counter <<- debug_eval_counter + 1L
      if (isTRUE(debug_par) && isTRUE(debug_par_each_eval) &&
          debug_eval_counter <= debug_par_max_eval) {
        debug_print_growth_par(par, paste0("eval ", debug_eval_counter, " update_gh=", update_gh_eval))
      }
      t0_fn <- proc.time()[["elapsed"]]
      val <- loglik_call(par, update_grid = update_gh_eval)
      if (isTRUE(debug_par) && (!is.finite(val) || val >= 1e9)) {
        cat("BAD objective at eval", debug_eval_counter, "value =", val, "\n")
        if (!(isTRUE(debug_par_each_eval) && debug_eval_counter <= debug_par_max_eval)) {
          debug_print_growth_par(par, paste0("BAD eval ", debug_eval_counter))
        }
      }
      if (isTRUE(profile)) {
        nfn <<- nfn + 1L
        fn_elapsed <<- c(fn_elapsed, proc.time()[["elapsed"]] - t0_fn)
      }
      val
    }
  }

  eval_ll_growth <- function(par, update_gh_eval) {
    t0_fn <- proc.time()[["elapsed"]]
    val <- loglik_call(par, update_grid = update_gh_eval)
    if (isTRUE(profile)) {
      nfn <<- nfn + 1L
      fn_elapsed <<- c(fn_elapsed, proc.time()[["elapsed"]] - t0_fn)
    }
    val
  }

  agh_trace <- NULL

  if (maxiter == 0L) {
    fn0 <- timed_fn(if (is_agh) 1L else 0L)
    val0 <- fn0(b0)
    best <- b0
    opt <- list(
      par = best,
      value = val0,
      convergence = NA_integer_,
      message = "single likelihood evaluation only"
    )
    t1 <- structure(
      c(0, 0, 0, 0, 0),
      names = c("user.self", "sys.self", "elapsed", "user.child", "sys.child")
    )
  } else if (!is_agh || is.null(update_gh) || update_gh <= 0L) {
    t0_opt <- proc.time()
    opt <- optim(
      par = b0,
      fn = timed_fn(if (is_agh) 1L else 0L),
      method = "BFGS",
      control = list(
        maxit = maxiter,
        trace = if (verbose) 1 else 0,
        parscale = ps0
      )
    )
    t1 <- proc.time() - t0_opt
    best <- opt$par
  } else {
    if (isTRUE(verbose)) cat("\nfast AGH optimization: outer-inner scheme\n")

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
      rel_loglik = numeric(),
      rel_par = numeric(),
      convergence = integer(),
      nfn = integer()
    )

    t0_opt <- proc.time()
    prev_refreshed <- NA_real_
    prev_best <- NULL
    agh_converged <- FALSE
    agh_failed <- FALSE
    agh_fail_message <- NULL

    for (oo in seq_len(agh_outer)) {
      inner_maxit <- oo * agh_step

      if (isTRUE(verbose)) {
        cat("\nAGH outer:", oo, "\n")
        cat("updating AGH grid...\n")
      }

      val_grid <- eval_ll_growth(best, 1L)
      if (!is.finite(val_grid) || val_grid >= 1e8) {
        agh_failed <- TRUE
        agh_fail_message <- paste0("AGH grid update failed at outer ", oo)
        break
      }

      if (isTRUE(verbose)) {
        cat("initial fixed-grid negloglik:", val_grid, "\n")
        cat("optimizing fixed AGH grid for maxit =", inner_maxit, "\n")
        if (isTRUE(use_analytic_score_opt)) cat("using analytic fixed-grid score\n")
      }

      if (isTRUE(use_analytic_score_opt)) {
        if (identical(methInteg, "AGH_RE_FULL_EQUIV")) {
          score_grid <- growth_tt_get_agh_re_full_equiv_grid()
          if (score_grid$fit_id0_c != fit_id0_c) {
            agh_failed <- TRUE
            agh_fail_message <- paste0("fit_id score grid mismatch at outer ", oo)
            break
          }
          growth_tt_set_agh_re_full_equiv_score_grid(score_grid)
        }
        reset_score_opt_cache()

        score_start <- score_opt_eval(best)
        score_start_diff <- abs(score_start$nll - val_grid) / max(1, abs(val_grid))
        if (!isTRUE(score_start$ok) || !is.finite(score_start_diff) || score_start_diff > 1e-6) {
          agh_failed <- TRUE
          agh_fail_message <- paste0(
            "analytic score objective mismatch at outer ", oo,
            ": relative difference = ", format(score_start_diff, scientific = TRUE)
          )
          break
        }

        opt_inner <- tryCatch(
          optim(
            par = best,
            fn = fn_fixed_analytic,
            gr = gr_fixed_analytic,
            method = "BFGS",
            control = list(
              maxit = inner_maxit,
              trace = if (verbose) 1 else 0,
              parscale = ps0
            )
          ),
          error = function(e) e
        )
      } else {
        opt_inner <- tryCatch(
          optim(
            par = best,
            fn = timed_fn(0L),
            method = "BFGS",
            control = list(
              maxit = inner_maxit,
              trace = if (verbose) 1 else 0,
              parscale = ps0
            )
          ),
          error = function(e) e
        )
      }

      if (inherits(opt_inner, "error")) {
        agh_failed <- TRUE
        agh_fail_message <- paste0(
          "AGH fixed-grid optimization failed at outer ", oo, ": ", opt_inner$message
        )
        break
      }

      best <- opt_inner$par
      val_fixed <- eval_ll_growth(best, 0L)
      val_refreshed <- eval_ll_growth(best, 1L)

      if (!is.finite(val_refreshed) || val_refreshed >= 1e8) {
        agh_failed <- TRUE
        agh_fail_message <- paste0("final AGH grid refresh failed at outer ", oo)
        break
      }

      rel_loglik <- NA_real_
      rel_par <- NA_real_
      if (is.finite(prev_refreshed) && !is.null(prev_best)) {
        rel_loglik <- abs(val_refreshed - prev_refreshed) /
          max(1, abs(val_refreshed), abs(prev_refreshed))
        rel_par <- max(
          abs(best - prev_best) /
            pmax(1, abs(best), abs(prev_best))
        )
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
          rel_loglik = rel_loglik,
          rel_par = rel_par,
          convergence = opt_inner$convergence,
          nfn = nfn
        )
      )

      if (isTRUE(verbose)) {
        cat("outer", oo, "inner maxit:", inner_maxit, "\n")
        cat("outer", oo, "fixed-grid negloglik:", val_fixed, "\n")
        cat("outer", oo, "refreshed AGH negloglik:", val_refreshed, "\n")
        cat("outer", oo, "refreshed AGH loglik:", -val_refreshed, "\n")
        if (is.finite(rel_loglik)) cat("outer", oo, "relative likelihood change:", rel_loglik, "\n")
        if (is.finite(rel_par)) cat("outer", oo, "relative parameter change:", rel_par, "\n")
      }

      if (is.finite(rel_loglik) && is.finite(rel_par) &&
          rel_loglik < tol_loglik && rel_par < tol_par) {
        agh_converged <- TRUE
        break
      }

      prev_refreshed <- val_refreshed
      prev_best <- best
    }

    if (isTRUE(agh_failed) || is.null(agh_trace) || nrow(agh_trace) == 0L) {
      stop(if (!is.null(agh_fail_message)) agh_fail_message else
             "AGH fast failed before producing a valid iteration.")
    }

    t1 <- proc.time() - t0_opt
    final_value <- tail(agh_trace$negloglik_refreshed, 1L)
    opt <- list(
      par = best,
      value = final_value,
      convergence = if (isTRUE(agh_converged)) 0L else 1L,
      message = if (isTRUE(agh_converged)) {
        "AGH fast outer-inner optimization converged"
      } else {
        "AGH fast outer-inner optimization reached maximum outer iterations"
      },
      counts = setNames(c(nfn, NA_integer_), c("function", "gradient"))
    )
  }

  if (!valid_objective(opt$value)) {
    stop("The growth likelihood could not be evaluated. Check the starting parameters and the design.")
  }
  par_opt <- setNames(best, nm_opt)

  agh_grid <- NULL

if (is_agh) {

  val_final_refresh <- loglik_call(
    par_opt,
    update_grid = 1L
  )

  if (!is.finite(val_final_refresh) ||
      val_final_refresh >= 1e8) {
    stop("Refresh finale AGH fallito.")
  }

  opt$value <- val_final_refresh

  if (joint_rule == 2L) {
    agh_grid <- joint_grid
  } else if (methInteg == "AGH") {

    agh_grid <- growth_tt_get_agh_grid()

    if (agh_grid$fit_id0_c != fit_id0_c) {
      stop(
        "fit_id della griglia AGH diverso dal fit corrente."
      )
    }

  } else if (methInteg == "AGH_FULL_EQUIV") {

    agh_grid <- growth_tt_get_agh_full_equiv_grid()

    if (agh_grid$fit_id0_c != fit_id0_c) {
      stop(
        "fit_id della griglia AGH full-equivalent diverso dal fit corrente."
      )
    }

  } else if (methInteg == "AGH_RE_FULL_EQUIV") {

    agh_grid <- growth_tt_get_agh_re_full_equiv_grid()

    if (agh_grid$fit_id0_c != fit_id0_c) {
      stop(
        "fit_id della griglia AGH_RE_FULL_EQUIV diverso dal fit corrente."
      )
    }
  }
}

nfn_opt <- nfn
hess <- NULL
hessian_raw <- NULL

t3 <- structure(
  c(0, 0, 0, 0, 0),
  names = c(
    "user.self",
    "sys.self",
    "elapsed",
    "user.child",
    "sys.child"
  )
)

nfn_hess <- 0L

if (hessian_method == "numeric_fixed" &&
    maxiter != 0L) {

  t2 <- proc.time()

  hess <- optimHess(
    par = par_opt,
    fn = timed_fn(0L)
  )

  hess <- (hess + t(hess)) / 2

  t3 <- proc.time() - t2
  nfn_hess <- nfn - nfn_opt
}

if (hessian_method == "analytic_fixed" &&
    maxiter != 0L) {

  if (!identification %in% c("baseline_std", "marker", "mixed_marker")) {
    stop(
      "analytic_fixed \u00e8 disponibile solo per identification='baseline_std', 'marker' o 'mixed_marker'."
    )
  }

  if (!methInteg %in% c("AGH_FULL_EQUIV", "AGH_RE_FULL_EQUIV")) {
    stop(
      "analytic_fixed \u00e8 disponibile solo per AGH_FULL_EQUIV o AGH_RE_FULL_EQUIV."
    )
  }

  expected_grid_type <- methInteg

  if (is.null(agh_grid) ||
      !identical(agh_grid$type, expected_grid_type)) {
    stop(
      "griglia ", expected_grid_type, " salvata mancante o incompatibile."
    )
  }

  score_name <- if (identical(methInteg, "AGH_FULL_EQUIV")) {
    "loglik_growth_tt_agh_full_equiv_score_q1"
  } else {
    "loglik_growth_re_tt_agh_full_equiv_score_q1"
  }

  if (!native_has(score_name)) {
    stop("routine ", score_name, " non caricata.")
  }

  if (identical(methInteg, "AGH_FULL_EQUIV")) {
    growth_tt_set_agh_full_equiv_grid(agh_grid)
  } else {
    growth_tt_set_agh_re_full_equiv_score_grid(agh_grid)
  }

  sym_score <- native_symbol(score_name)

  score_fixed <- function(par) {

    .growth_set_fortran_options(
      identification,
      threshold_identification,
      growth_shape,
      intercept_mask,
      gamma_dep,
      x_slope_col,
      nvar,
      n,
      Q,
      p_fix,
      rho_fixed = rho_fixed_vec,
      marker_index = marker_index
    )

    if (identical(methInteg, "AGH_FULL_EQUIV")) {
      z <- native_call(
        sym_score,
        b0 = as.double(par),
        Y0 = Y0_c,
        nvar0 = nvar0_c,
        nsub0 = nsub0_c,
        nobs0 = nobs0_c,
        n0 = n0_c,
        nmes0 = nmes0_c,
        c0 = c0_c,
        cm10 = cm10_c,
        ntre0 = ntre0_c,
        nquad_eta0 = nquad_eta0_c,
        nodes_eta0 = nodes_eta0_c,
        weights_eta0 = weights_eta0_c,
        nquad_re0 = nquad_re0_c,
        nodes_re0 = nodes_re0_c,
        weights_re0 = weights_re0_c,
        nitem0 = nitem0_c,
        npm0 = npm0_c,
        re_item0 = re_item0,
        theta_time0 = theta_time0,
        delta_free0 = delta_free0,
        monotone_tau0 = monotone_tau0,
        update_gh0 = 0L,
        fit_id0 = fit_id0_c,
        Q0 = Q0_c,
        loadmask0 = loadmask0_c,
        p_fix0 = as.integer(
          if (isTRUE(design$subject_specific_design)) p_fix else -p_fix
        ),
        p_ran0 = p_ran0_c,
        X0 = X0_c,
        Z0 = Z0_c,
        delta_mask0 = delta_mask0_c,
        delta_x_col0 = dx_c,
        delta_z_col0 = dz_c,
        vrais = double(1L),
        score0 = double(npm0),
        ok0 = integer(1L)
      )
    } else {
      z <- native_call(
        sym_score,
        b0 = as.double(par),
        Y0 = Y0_c,
        nvar0 = nvar0_c,
        nsub0 = nsub0_c,
        nobs0 = nobs0_c,
        n0 = n0_c,
        nmes0 = nmes0_c,
        c0 = c0_c,
        cm10 = cm10_c,
        ntre0 = ntre0_c,
        nquad_u0 = nquad_u0_c,
        nodes_u0 = nodes_u0_c,
        weights_u0 = weights_u0_c,
        nquad_zeta0 = nquad_zeta0_c,
        nodes_zeta0 = nodes_zeta0_c,
        weights_zeta0 = weights_zeta0_c,
        nquad_re0 = nquad_re0_c,
        nodes_re0 = nodes_re0_c,
        weights_re0 = weights_re0_c,
        nitem0 = nitem0_c,
        npm0 = npm0_c,
        re_item0 = re_item0,
        theta_time0 = theta_time0,
        delta_free0 = delta_free0,
        monotone_tau0 = monotone_tau0,
        update_gh0 = 0L,
        fit_id0 = fit_id0_c,
        Q0 = Q0_c,
        loadmask0 = loadmask0_c,
        p_fix0 = as.integer(
          if (isTRUE(design$subject_specific_design)) p_fix else -p_fix
        ),
        p_ran0 = p_ran0_c,
        X0 = X0_c,
        Z0 = Z0_c,
        delta_mask0 = delta_mask0_c,
        delta_x_col0 = dx_c,
        delta_z_col0 = dz_c,
        vrais = double(1L),
        score0 = double(npm0),
        ok0 = integer(1L)
      )
    }

    if (length(z$ok0) != 1L || z$ok0 != 1L ||
        !is.finite(z$vrais) || z$vrais >= 1e8 ||
        length(z$score0) != npm0 || any(!is.finite(z$score0))) {
      stop("score analitico fallito.")
    }

    as.numeric(z$score0)
  }

  t2 <- proc.time()

  h <- hess_eps * (
    abs(par_opt) + 1
  )

  hessian_raw <- matrix(
    NA_real_,
    npm0,
    npm0,
    dimnames = list(
      nm_opt,
      nm_opt
    )
  )

  for (j in seq_len(npm0)) {

    pp <- par_opt
    pm <- par_opt

    pp[j] <- pp[j] + h[j]
    pm[j] <- pm[j] - h[j]

    sp <- score_fixed(pp)
    sm <- score_fixed(pm)

    hessian_raw[, j] <- (
      sp - sm
    ) / (
      2 * h[j]
    )
  }

  hess <- (
    hessian_raw +
      t(hessian_raw)
  ) / 2

  t3 <- proc.time() - t2

  nfn_hess <- 2L * npm0
}
  timing <- list(
    opt_elapsed = unname(t1["elapsed"]),
    opt_user = unname(t1["user.self"]),
    opt_system = unname(t1["sys.self"]),
    hess_elapsed = unname(t3["elapsed"]),
    hess_user = unname(t3["user.self"]),
    hess_system = unname(t3["sys.self"]),
    total_elapsed = unname(t1["elapsed"]) + unname(t3["elapsed"])
  )

  profile_info <- NULL
  if (isTRUE(profile)) {
    profile_info <- list(
      npar = length(b0),
      nfn_total = nfn,
      nfn_opt = nfn_opt,
      nfn_hess = nfn_hess,
      fn_mean = if (length(fn_elapsed)) mean(fn_elapsed) else NA_real_,
      fn_median = if (length(fn_elapsed)) median(fn_elapsed) else NA_real_,
      fn_min = if (length(fn_elapsed)) min(fn_elapsed) else NA_real_,
      fn_max = if (length(fn_elapsed)) max(fn_elapsed) else NA_real_,
      fn_sum = if (length(fn_elapsed)) sum(fn_elapsed) else 0,
      parscale_summary = c(
        min = min(ps0), median = median(ps0), mean = mean(ps0), max = max(ps0)
      )
    )

    if (isTRUE(profile_print)) {
      cat("\nprofiling summary\n")
      cat("nsub:", nsub, "\n")
      cat("npar:", profile_info$npar, "\n")
      cat("model: generalized_growth_tt_flexible\n")
      cat("identification:", identification, "\n")
      if (identification %in% c("marker", "mixed_marker")) {
        cat("marker_item:", paste(sprintf("%s=%s", factor_names, marker_names), collapse = ", "), "\n")
      }
      cat("threshold_identification:", threshold_identification, "\n")
      cat("growth_shape:", growth_shape, "\n")
      cat("methInteg:", methInteg, "\n")
      cat("n:", n, "\n")
      cat("Q:", Q, "\n")
      cat("dim_eta:", dim_eta, "\n")
      if (two_tier) {
        cat("nquad_eta:", eta_rule$nquad, "\n")
      } else {
        cat("dim_joint:", dim_joint, "\n")
        cat("nquad_joint:", eta_rule$nquad, "\n")
      }
      if (is_growth_re_agh) {
        cat("dim_u:", dim_u, "\n")
        cat("nquad_u:", u_rule$nquad, "\n")
        cat("dim_zeta:", dim_zeta, "\n")
        cat("nquad_zeta:", zeta_rule$nquad, "\n")
      }
      if (two_tier) cat("nquad_re:", re_rule$nquad, "\n")
      cat("nfn_total:", profile_info$nfn_total, "\n")
      cat("nfn_opt:", profile_info$nfn_opt, "\n")
      cat("nfn_hess:", profile_info$nfn_hess, "\n")
      cat("fn_mean:", sprintf("%.6f", profile_info$fn_mean), "\n")
      cat("fn_median:", sprintf("%.6f", profile_info$fn_median), "\n")
      cat("fn_min:", sprintf("%.6f", profile_info$fn_min), "\n")
      cat("fn_max:", sprintf("%.6f", profile_info$fn_max), "\n")
      cat("fn_sum:", sprintf("%.6f", profile_info$fn_sum), "\n")
      cat("optim_elapsed:", sprintf("%.6f", timing$opt_elapsed), "\n")
      cat("hess_elapsed:", sprintf("%.6f", timing$hess_elapsed), "\n")
      cat("total_elapsed:", sprintf("%.6f", timing$total_elapsed), "\n\n")
      flush.console()
    }
  }

  vcov_opt <- if (!is.null(hess)) tryCatch(solve(hess), error=function(e) NULL) else NULL
  se_opt <- setNames(rep(NA_real_, npm0), nm_opt)
  if (!is.null(vcov_opt)) {
    dd <- diag(vcov_opt)
    se_opt[] <- ifelse(is.finite(dd) & dd >= 0, sqrt(dd), NA_real_)
  }

  dec <- decode(par_opt)
  if (identification %in% c("marker", "mixed_marker")) {
    marker_loading <- vapply(
      seq_len(Q),
      function(q) dec$lambda[marker_index[q], q],
      numeric(1L)
    )
    if (any(!is.finite(marker_loading)) || any(abs(marker_loading - 1) > 1e-12)) {
      stop("errore interno: il loading marker non \u00e8 fissato esattamente a 1.")
    }
  }
  if (identification == "marker") {
    marker_intercept <- dec$gamma[unique(marker_index)]
    if (any(!is.finite(marker_intercept)) || any(abs(marker_intercept) > 1e-12)) {
      stop("errore interno: l'intercetta dell'item marker non \u00e8 fissata esattamente a 0.")
    }
  }

  par_natural_from_decode <- function(d) {
    out <- numeric(0)
    for (k in seq_len(nvar)) {
      th <- d$tau[[k]]
      out[names(th)] <- as.numeric(th)
    }
    for (q in seq_len(Q)) for (k in which(loadmask[, q] == 1L)) {
      out[sprintf("lambda[%s,%s]", namvar[k], factor_names[q])] <- d$lambda[k, q]
    }
    for (k in seq_len(nvar)) {
      out[sprintf("gamma[%s]", namvar[k])] <- d$gamma[k]
      for (tt in seq_len(n)) {
        out[sprintf("gamma[%s,t%d]", namvar[k], tt)] <- d$gamma_time[k, tt]
      }
    }
    for (k in seq_len(nvar)) for (tt in seq_len(n)) {
      out[sprintf("theta[%s,t%d]", namvar[k], tt)] <- d$theta[k, tt]
    }
    for (r in seq_len(p_fix)) for (q in seq_len(Q)) {
      nm0 <- design$fixed_names[[q]][r]
      lab <- if (r == 1L) "alpha0" else if (r == 2L) "alpha1" else paste0("beta[", nm0, ",", factor_names[q], "]")
      if (r <= 2L) {
        out[sprintf("%s[%s]", lab, factor_names[q])] <- d$beta[q, r]
        if (r == 2L && growth_shape == "rho") {
          out[sprintf("muS[%s]", factor_names[q])] <- d$beta[q, r]
        }
      } else {
        out[lab] <- d$beta[q, r]
      }
    }
    if (!is.null(rho_spec)) {
      out[as.vector(rho_spec$names)] <- as.vector(d$rho)
    } else {
      for (q in seq_len(Q)) {
        if (growth_shape == "rho") out[sprintf("rho[%s]", factor_names[q])] <- d$rho[q]
        else out[sprintf("delta[%s]", factor_names[q])] <- d$delta[q]
      }
    }
    for (cc in seq_len(nrow(d$Lu))) for (rr in cc:nrow(d$Lu)) {
      out[sprintf("Lu[%s,%s]", rownames(d$Lu)[rr], colnames(d$Lu)[cc])] <- d$Lu[rr, cc]
    }
    for (r in seq_len(nrow(d$Sigma_u))) {
      out[sprintf("Var_u[%s]", rownames(d$Sigma_u)[r])] <- d$Sigma_u[r, r]
    }
    if (nrow(d$Sigma_u) >= 2L) {
      for (cc in seq_len(nrow(d$Sigma_u) - 1L)) for (rr in (cc + 1L):nrow(d$Sigma_u)) {
        out[sprintf("Cov_u[%s,%s]", rownames(d$Sigma_u)[rr], colnames(d$Sigma_u)[cc])] <- d$Sigma_u[rr, cc]
      }
    }
    for (tt in seq_len(n)) for (q in seq_len(Q)) {
      out[sprintf("zeta_var[t%d,%s]", tt, factor_names[q])] <- d$zeta_var[tt, q]
    }
    if (!is.null(d$sigmak)) for (k in seq_len(nvar)) {
      out[sprintf("sigmak[%s]", namvar[k])] <- d$sigmak[k]
    }
    for (q in seq_len(Q)) out[sprintf("baseline_var[%s]", factor_names[q])] <- d$baseline_var[q]
    for (ii in seq_len(nsub)) for (tt in seq_len(n)) for (q in seq_len(Q)) {
      if (!is.null(d$score)) {
        out[sprintf("growth_score[%s,t%d,%s]", IDs[ii], tt, factor_names[q])] <- d$score[ii, tt, q]
      }
      out[sprintf("kappa[%s,t%d,%s]", IDs[ii], tt, factor_names[q])] <- d$kappa[ii, tt, q]
    }
    out
  }

  par_nat <- par_natural_from_decode(dec)
  par_natural <- par_nat

  fixed_nat_names <- character(0)
  for (k in seq_len(nvar)) {
    m <- cm1[k]
    if (m == 1L) {
      fixed_nat_names <- c(fixed_nat_names, sprintf("tau[%s]_c1", namvar[k]))
    } else if (threshold_identification == "first_two") {
      fixed_nat_names <- c(fixed_nat_names,
                           sprintf("tau[%s]_c1", namvar[k]),
                           sprintf("tau[%s]_c2", namvar[k]))
    } else {
      j0 <- (m + 1L) %/% 2L
      fixed_nat_names <- c(fixed_nat_names,
                           sprintf("tau[%s]_c%d", namvar[k], j0),
                           sprintf("tau[%s]_c%d", namvar[k], j0 + 1L))
    }
    if (m == 1L) fixed_nat_names <- c(fixed_nat_names, sprintf("theta[%s,t%d]", namvar[k], seq_len(n)))
  }
  if (identification %in% c("baseline_std", "mixed_marker")) {
    fixed_nat_names <- c(fixed_nat_names, sprintf("alpha0[%s]", factor_names))
  }
  if (identification == "baseline_std") {
    fixed_nat_names <- c(fixed_nat_names, sprintf("baseline_var[%s]", factor_names))
  }
  if (identification %in% c("marker", "mixed_marker")) {
    fixed_nat_names <- c(
      fixed_nat_names,
      sprintf("lambda[%s,%s]", marker_names, factor_names)
    )
  }
  if (identification == "marker") {
    fixed_nat_names <- c(
      fixed_nat_names,
      sprintf("gamma[%s]", marker_names),
      sprintf("gamma[%s,t1]", marker_names)
    )
  }
  if (!is.null(rho_spec)) {
    fixed_nat_names <- c(fixed_nat_names, rho_spec$names[!rho_spec$free])
  } else if (!isTRUE(delta_free)) {
    fixed_nat_names <- c(
      fixed_nat_names,
      sprintf("%s[%s]", if (growth_shape == "rho") "rho" else "delta", factor_names)
    )
  }
  fixed_nat_names <- intersect(unique(fixed_nat_names), names(par_nat))

  numeric_jacobian_decode <- function(par, f, eps = 1e-5) {
    f0 <- f(par)
    J <- matrix(NA_real_, nrow = length(f0), ncol = length(par),
                dimnames = list(names(f0), names(par)))
    for (j in seq_along(par)) {
      h <- eps * (abs(par[j]) + 1)
      pp <- par; pm <- par
      pp[j] <- pp[j] + h; pm[j] <- pm[j] - h
      J[, j] <- (f(pp) - f(pm)) / (2 * h)
    }
    J
  }

  vcov_nat <- NULL
  se_nat_all <- setNames(rep(NA_real_, length(par_nat)), names(par_nat))
  if (!is.null(vcov_opt)) {
    J_nat <- numeric_jacobian_decode(par_opt, function(p) par_natural_from_decode(decode(p)))
    vcov_nat <- J_nat %*% vcov_opt %*% t(J_nat)
    vcov_nat <- (vcov_nat + t(vcov_nat)) / 2
    dnat <- diag(vcov_nat)
    se_nat_all[] <- ifelse(is.finite(dnat) & dnat >= 0, sqrt(dnat), NA_real_)
  }
  se_nat_all[fixed_nat_names] <- 0

  make_wald <- function(est, se) {
    zcrit <- qnorm(0.975)
    z <- ifelse(is.finite(se) & se > 0, est / se, NA_real_)
    pval <- ifelse(is.finite(z), 2 * pnorm(abs(z), lower.tail = FALSE), NA_real_)
    data.frame(
      name = names(est), est = as.numeric(est), se = as.numeric(se),
      z = as.numeric(z), p = as.numeric(pval),
      ci_low = as.numeric(ifelse(is.finite(se) & se > 0, est - zcrit * se, ifelse(se == 0, est, NA_real_))),
      ci_high = as.numeric(ifelse(is.finite(se) & se > 0, est + zcrit * se, ifelse(se == 0, est, NA_real_))),
      row.names = NULL, check.names = FALSE
    )
  }
  wald_all <- make_wald(par_nat, se_nat_all)
  wald_nat <- wald_all[is.finite(wald_all$se) & wald_all$se > 0, , drop = FALSE]

  c_inputs <- list(Y0_c=Y0_c,nvar0_c=nvar0_c,nsub0_c=nsub0_c,nobs0_c=nobs0_c,n0_c=n0_c,
                   nmes0_c=nmes0_c,c0_c=c0_c,cm10_c=cm10_c,ntre0_c=ntre0_c,
                   nquad_eta0_c=nquad_eta0_c,nodes_eta0_c=nodes_eta0_c,weights_eta0_c=weights_eta0_c,
                   nquad_u0_c=nquad_u0_c,nodes_u0_c=nodes_u0_c,weights_u0_c=weights_u0_c,
                   nquad_zeta0_c=nquad_zeta0_c,nodes_zeta0_c=nodes_zeta0_c,weights_zeta0_c=weights_zeta0_c,
                   nquad_re0_c=nquad_re0_c,nodes_re0_c=nodes_re0_c,weights_re0_c=weights_re0_c,
                   nitem0_c=nitem0_c,npm0_c=npm0_c,re_item0_c=re_item0,theta_time0_c=theta_time0,
                   delta_free0_c=delta_free0,monotone_tau0_c=monotone_tau0,Q0_c=Q0_c,loadmask0_c=loadmask0_c,
                   p_fix0_c=p_fix0_c,p_ran0_c=p_ran0_c,X0_c=X0_c,Z0_c=Z0_c,delta_mask0_c=delta_mask0_c,
                   delta_x_col0_c=dx_c,delta_z_col0_c=dz_c,fit_id0_c=fit_id0_c,
                   subject_specific_design0_c=as.integer(design$subject_specific_design),
                   IDs=IDs,namvar=namvar,occ=occ,ID=ID,methInteg=methInteg,
                   analytic_score=as.integer(use_analytic_score_opt),rho_spec=rho_spec,
                   joint_rule_c=joint_rule,
                   dims=list(nsub=nsub,nvar=nvar,n=n,Q=Q,dim_eta=dim_eta,dim_u=dim_u,dim_zeta=dim_zeta,
                             p_fix=p_fix,p_ran=p_ran,npm0=npm0),
                   growth_options=list(identification=identification,threshold_identification=threshold_identification,
                                       growth_shape=growth_shape,rho_fixed=rho_fixed_vec,
                                       marker_index=marker_index,marker_item=marker_names,
                                       intercept_mask=intercept_mask,gamma_dep=gamma_dep,
                                       x_slope_col=x_slope_col,zero_time=zero_time))

  rho_table <- if (!is.null(rho_spec)) {
    data.frame(occasion = rep(occ_vals, Q), factor = rep(factor_names, each = n),
               score = as.vector(dec$rho), estimated = as.vector(rho_spec$free),
               se = unname(se_nat_all[as.vector(rho_spec$names)]))
  } else NULL

  new_long_fit(list(
    conv=opt$convergence,convergence=opt$convergence,message=opt$message,loglik=-opt$value,negloglik=opt$value,
    BIC=2 * opt$value + length(par_opt) * log(nsub),
    fit_id0_c=fit_id0_c,c_inputs=c_inputs,agh_grid=agh_grid,
    timing=timing,parscale=ps0,
    par_opt=par_opt,par=par_nat,par_natural=par_natural,
    hessian=hess,se_opt=se_opt,vcov_opt=vcov_opt,vcov=vcov_nat,se=se_nat_all,
    wald_nat=wald_nat,wald_all=wald_all,
    dims=list(nsub=nsub,nvar=nvar,n=n,Q=Q,dim_eta=dim_eta,dim_u=dim_u,dim_zeta=dim_zeta,
              p_fix=p_fix,p_ran=p_ran,npm0=npm0),
    structure=list(
      model="generalized_growth_tt",
      model_extended="generalized_growth_tt_flexible",
      identification=identification,
      marker_item=marker_names,
      marker_index=marker_index,
      threshold_identification=threshold_identification,
      intercept_time_free=intercept_mask,
      growth_shape=growth_shape,
      rho_pattern=if (!is.null(rho_spec)) rho_spec$pattern else NULL,
      zero_anchor=occ_vals[zero_time],
      loadmask=loadmask,methInteg=methInteg,
      analytic_score=isTRUE(use_analytic_score_opt),
      re_item=re_item0,two_tier=two_tier,theta_time=theta_time0,hessian_method=hessian_method,
      delta_free=delta_free0,rho_fixed=rho_fixed_vec,monotone_tau=monotone_tau0,
      growth=list(
        type=if(isTRUE(simple_growth)) "simple_random_intercept_slope" else "general_XZ_growth",
        simple_growth=isTRUE(simple_growth),same_design=isTRUE(same_design),
        subject_specific_design=isTRUE(design$subject_specific_design),design_nsub=design$design_nsub,
        p_fix=p_fix,p_ran=p_ran,fixed_names=design$fixed_names,random_names=design$random_names,
        time_scores=design$time,delta_mask=design$delta_mask,dimU=dim_u,factor_names=factor_names,
        intercept_covariates=intercept_covariates,slope_covariates=slope_covariates
      ),
      integration_target=if(!two_tier) "joint" else if(is_growth_re_agh) "growth_random_effects" else "eta",
      nquad=if(!two_tier)nquad_eta0_c else NULL,nquad_eta=if(two_tier)nquad_eta0_c else NULL,
      nquad_u=if(is_growth_re_agh)nquad_u0_c else NULL,
      nquad_zeta=if(is_growth_re_agh)nquad_zeta0_c else NULL,nquad_re=if(two_tier)nquad_re0_c else NULL,
      dimInt=if(!two_tier)dim_joint else NULL,dim_eta_int=dim_eta,dim_u_int=if(is_growth_re_agh)dim_u else NULL,
      dim_zeta_int=if(is_growth_re_agh)dim_zeta else NULL,
      data_report=data_report
    ),
    decode=dec,rho=rho_table,
    profile=profile_info,
    debug=list(
      nm_opt=nm_opt,nfn=nfn,fn_elapsed=fn_elapsed,agh_trace=agh_trace,
      debug=list(nm_opt=nm_opt,nfn=nfn,fn_elapsed=fn_elapsed,agh_trace=agh_trace)
    )
  ), "Ordinal latent growth model", packed)
}

#use the general growth fitter
growth_tt_fit <- growth_tt_fit_flexible
