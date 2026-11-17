check_count <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || x != floor(x) || x > .Machine$integer.max) {
    stop(name, " must be a whole number of at least ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

check_flag <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !x %in% c(FALSE, TRUE)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  as.logical(x)
}

make_gh_1d <- function(nquadgh) {
  n <- check_count(nquadgh, "nquadgh")
  j <- seq_len(n - 1L)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(j, j + 1L)] <- jacobi[cbind(j + 1L, j)] <- sqrt(j)
  rule <- eigen(jacobi, symmetric = TRUE)
  order <- seq.int(n, 1L)
  weights <- rule$vectors[1L, order]^2
  list(nquad = n, nodes = rule$values[order], weights = weights / sum(weights))
}

make_gh_tensor <- function(dim, nquadgh) {
  dim <- check_count(dim, "quadrature dimension")
  n <- check_count(nquadgh, "nquadgh")
  total <- as.double(n)^dim
  if (total * dim > .Machine$integer.max) {
    stop("The tensor rule is too large for the native interface. Reduce nquadgh or use a lower-dimensional rule.",
         call. = FALSE)
  }

  rule <- make_gh_1d(n)
  nodes <- matrix(0, dim, total)
  weights <- rep(1, total)
  for (d in seq_len(dim)) {
    indices <- rep(rep(seq_len(n), each = n^(dim - d)), times = n^(d - 1L))
    nodes[d, ] <- rule$nodes[indices]
    weights <- weights * rule$weights[indices]
  }
  list(nquad = as.integer(total), nodes = as.double(nodes), weights = weights)
}

make_qmc <- function(dim, nquadqmc) {
  dim <- check_count(dim, "quadrature dimension")
  nquadqmc <- check_count(nquadqmc, "nquadqmc")
  if (!requireNamespace("spacefillr", quietly = TRUE)) {
    stop("Install spacefillr to use QMC, or select GH or AGH.", call. = FALSE)
  }
  u <- spacefillr::generate_sobol_owen_set(nquadqmc, dim)
  u <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
  nodes <- matrix(qnorm(u), nrow = nquadqmc, ncol = dim)
  list(nquad = nquadqmc, nodes = as.double(t(nodes)),
       weights = rep(1 / nquadqmc, nquadqmc))
}

make_qmc_1d <- function(nquadqmc_re) {
  make_qmc(1L, nquadqmc_re)
}
