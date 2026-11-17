new_long_fit <- function(fit, family, packed) {
  fit$family <- family
  fit$nobs <- packed$informative
  fit$occasions <- packed$times
  fit$categories <- packed$categories
  fit$integration <- describe_integration(fit$dims, fit$structure)
  structure(fit, class = "longProbitLVM_fit")
}

print.longProbitLVM_fit <- function(x, ...) {
  cat(x$family, "\n", sep = "")
  cat("Log likelihood ", format(x$loglik), " | ", length(x$par_opt),
      " parameters | ", x$nobs, " subjects with observed items\n", sep = "")
  if (!is.null(x$message)) cat(x$message, "\n", sep = "")
  invisible(x)
}

coef.longProbitLVM_fit <- function(object, ..., natural = FALSE) {
  if (isTRUE(natural)) {
    if (is.null(object$par_natural)) {
      stop("Natural parameters are available in the fitted object's decode field.", call. = FALSE)
    }
    return(object$par_natural)
  }
  object$par_opt
}

logLik.longProbitLVM_fit <- function(object, ...) {
  structure(object$loglik, df = length(object$par_opt),
            nobs = nobs(object), class = "logLik")
}

nobs.longProbitLVM_fit <- function(object, ...) {
  object$nobs
}

vcov.longProbitLVM_fit <- function(object, ...) {
  object$vcov_opt
}

next_fit_id <- local({
  value <- 0L
  function() {
    value <<- if (value == .Machine$integer.max) 1L else value + 1L
    value
  }
})

valid_objective <- function(x) {
  is.finite(x) && x < 1e8
}

#report quadrature dimensions and node counts
integration_info <- function(fit) {
  if (!inherits(fit, "longProbitLVM_fit") || is.null(fit$integration)) {
    stop("Supply a fitted model with quadrature information.", call. = FALSE)
  }
  fit$integration
}

describe_integration <- function(dimensions, structure) {
  method <- structure$methInteg
  eta <- dimensions$n * dimensions$Q
  item <- if (isTRUE(as.logical(structure$re_item))) dimensions$nvar else 0L
  growth <- if (is.null(structure$dim_u_int)) 0L else structure$dim_u_int
  total <- eta + item + growth
  rule <- function(block, dimension, nodes, copies = 1L) {
    data.frame(block = block, dimension = as.integer(dimension), copies = as.integer(copies),
               nodes_per_rule = as.double(nodes))
  }

  if (!isTRUE(structure$two_tier)) {
    blocks <- rule("Joint", total, structure$nquad)
    mode <- if (method == "AGH") total else 0L
    layout <- "Standard joint"
  } else if (method == "AGH_RE_NESTED") {
    blocks <- rbind(rule("Growth coefficients", growth, structure$nquad_u),
                    rule("Occasion disturbances", eta, structure$nquad_zeta))
    mode <- growth
    layout <- "Two tier, nested growth effects"
  } else if (method == "AGH_RE_FULL_EQUIV") {
    blocks <- rule("Growth coefficients and occasion disturbances", growth + eta,
                   as.double(structure$nquad_u) * structure$nquad_zeta)
    mode <- total
    layout <- "Two tier, full-equivalent growth effects"
  } else {
    blocks <- rule("Latent trajectory", eta, structure$nquad_eta)
    full <- method == "AGH_FULL_EQUIV" || identical(structure$two_tier_type, "full_equiv") && method == "AGH"
    mode <- if (full) total else if (method == "AGH") eta else 0L
    layout <- if (full) "Two tier, full-equivalent" else if (method == "AGH") {
      "Two tier, node-specific"
    } else {
      "Two tier"
    }
  }
  if (item > 0L && isTRUE(structure$two_tier)) {
    blocks <- rbind(blocks, rule("Item intercept", 1L, structure$nquad_re, copies = item))
  }
  list(method = method, layout = layout, joint_dimension = as.integer(total),
       largest_rule_dimension = max(blocks$dimension), largest_mode_dimension = as.integer(mode),
       rules = blocks)
}
