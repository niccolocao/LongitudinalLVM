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
  rule <- fastGHQuad::gaussHermiteData(n)
  list(nquad = n, nodes = rule$x * sqrt(2), weights = rule$w / sqrt(pi))
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
  tensor <- function(values) {
    as.matrix(do.call(expand.grid, rep(list(values), dim)))[, dim:1L, drop = FALSE]
  }
  list(nquad = as.integer(total), nodes = as.double(t(tensor(rule$nodes))),
       weights = apply(tensor(rule$weights), 1L, prod))
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
