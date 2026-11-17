library(LongitudinalLVM)

set.seed(71)
n <- 30L
dat <- expand.grid(ID = seq_len(n), occ = 1:3)
t <- dat$occ - 1L
eta <- 0.3 + 0.4 * t + rep(rnorm(n, sd = 0.6), 3L) +
  t * rep(rnorm(n, sd = 0.3), 3L) + rnorm(nrow(dat), sd = 0.65)
items <- c("a", "b", "c")
for (item in items) {
  ystar <- eta + rep(rnorm(n, sd = 0.5), 3L) + rnorm(nrow(dat))
  dat[[item]] <- as.integer(cut(ystar, c(-Inf, 0, 1, 1.6, Inf)))
  dat[[item]][sample.int(nrow(dat), 14L)] <- NA_integer_
}
dat[dat$ID == 1L, items] <- NA_integer_
dat[dat$ID == 7L, "b"] <- NA_integer_
dat <- dat[sample.int(nrow(dat)), ]

fit <- growth_tt_fit(dat, items, methInteg = "AGH_FULL_EQUIV", nquadgh = 3L,
                     maxiter = 3L, delta_free = FALSE, monotone_tau = FALSE,
                     debug_par = FALSE, profile = FALSE)
parts <- growth_tt_individual_derivatives(fit, dat, compute_hessian = FALSE)
stopifnot(abs(sum(parts$loglik_i) - as.numeric(logLik(fit))) < 1e-8,
          abs(parts$loglik_i[["1"]]) < 1e-9, nobs(fit) == n - 1L,
          all(is.finite(parts$score_ind)), isTRUE(parts$settings$agh_fixed_grid),
          "Var_u[u0_F1]" %in% names(coef(fit, natural = TRUE)))
start <- append(coef(fit), c("delta[F1]" = 0), after = match("alpha1[F1]", names(coef(fit))))
full <- growth_tt_fit(dat, items, methInteg = "AGH_FULL_EQUIV",
                      nquadgh = 3L, maxiter = 0L, delta_free = TRUE,
                      monotone_tau = FALSE, start = start, debug_par = FALSE, profile = FALSE)
stopifnot(abs(as.numeric(logLik(full) - logLik(fit))) < 1e-5)

cfa <- function(method, two_tier = FALSE, type = "node_specific") {
  cfa_long_fit_md(dat, items, methInteg = method, two_tier = two_tier,
                  two_tier_type = type, nquadgh = 2L, nquadqmc = 32L,
                  maxiter = 0L, profile = FALSE)
}
joint <- cfa("GH")
split <- cfa("GH", TRUE)
stopifnot(abs(as.numeric(logLik(joint) - logLik(split))) < 1e-10)
for (tier in c(FALSE, TRUE)) stopifnot(is.finite(logLik(cfa("AGH", tier))))
stopifnot(is.finite(logLik(cfa("AGH", TRUE, "full_equiv"))))
if (requireNamespace("spacefillr", quietly = TRUE)) {
  stopifnot(is.finite(logLik(cfa("QMC"))))
}
stopifnot(inherits(try(cfa("QMC", TRUE), silent = TRUE), "try-error"))

method <- c("GH", "AGH", "AGH_RE_NESTED", "AGH_RE_FULL_EQUIV")
for (m in method) {
  g <- growth_tt_fit(dat, items, methInteg = m, nquadgh = 2L, maxiter = 0L,
                     debug_par = FALSE, profile = FALSE)
  stopifnot(is.finite(logLik(g)), nobs(g) == n - 1L)
}

wu <- measurement_invariance_wu2016(dat, items, nquadgh = 2L, maxiter = 0L)
stopifnot(nrow(wu$table) == 5L, all(is.finite(wu$table$logLik)),
          all(is.na(wu$table$p_vs_previous)))
restored <- growth_tt_individual_derivatives(fit, compute_hessian = FALSE, compute_score = FALSE)
stopifnot(identical(parts$loglik_i, restored$loglik_i))

growth <- function(method, two_tier = NULL, ...) {
  growth_tt_fit(dat, items, methInteg = method, two_tier = two_tier,
                nquadgh = 2L, nquadqmc = 32L, maxiter = 0L, debug_par = FALSE,
                profile = FALSE, growth_shape = "rho", threshold_identification = "central_01", ...)
}
reference <- fastGHQuad::gaussHermiteData(2L)
gh_grid <- function(v) as.matrix(expand.grid(v, v, v))[, 3:1, drop = FALSE]
growth_split <- growth("GH", TRUE)
growth_joint <- growth("GH", FALSE)
stopifnot(identical(growth_split$c_inputs$nodes_eta0_c,
                    as.double(t(gh_grid(reference$x * sqrt(2))))),
          identical(growth_split$c_inputs$weights_eta0_c,
                    apply(gh_grid(reference$w / sqrt(pi)), 1L, prod)),
          abs(growth_split$loglik - growth_joint$loglik) < 1e-10,
          growth_joint$BIC == -2 * growth_joint$loglik + length(coef(growth_joint)) * log(n))
growth_adaptive <- growth("AGH", FALSE)
stopifnot(abs(growth("AGH_FULL_EQUIV", TRUE)$loglik - growth_adaptive$loglik) < 1e-10,
          identical(integration_info(growth_adaptive)$largest_mode_dimension, 6L))

check_joint <- function(fit, count) {
  rule <- integration_info(fit)$rules
  individual <- growth_tt_individual_derivatives(fit, par_index = 1L, compute_hessian = FALSE)
  stopifnot(!fit$structure$two_tier, rule$dimension == 6L, rule$nodes_per_rule == count,
            abs(sum(individual$loglik_i) - fit$loglik) < 1e-10,
            individual$loglik_i[["1"]] == 0, all(is.finite(individual$score_ind)))
}
check_joint(growth_joint, 2^6)
check_joint(growth_adaptive, 2^6)
if (requireNamespace("spacefillr", quietly = TRUE)) check_joint(growth("QMC"), 32L)
stopifnot(inherits(try(growth("QMC", TRUE), silent = TRUE), "try-error"))
