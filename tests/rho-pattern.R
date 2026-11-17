library(LongitudinalLVM)

set.seed(13)
n <- 30L
visits <- c(0, 3, 6, 9, 12)
dat <- expand.grid(ID = seq_len(n), wave = visits)
dat$age <- rep(seq(-1, 1, length.out = n), length(visits))
items <- letters[1:4]
t <- match(dat$wave, visits)
for (k in items) {
  dat[[k]] <- as.integer(cut(t / 3 + rnorm(nrow(dat)), c(-Inf, 0, 1, 1.6, Inf)))
  dat[[k]][sample.int(nrow(dat), 14L)] <- NA_integer_
}
dat[dat$ID == 1L, items] <- NA_integer_
p <- c(0, 1, NA, NA, NA)
q <- c(0, NA, NA, NA, 1)
base <- list(data = dat, namvar = items, occ = "wave", nquadgh = 2L, maxiter = 3L)
f <- do.call(growth_tt_fit, c(base, list(rho_pattern = p, methInteg = "AGH_FULL_EQUIV")))
g <- do.call(growth_tt_fit, c(base, list(time_scores = q, methInteg = "GH")))
fixed_anchors <- function(x, pattern) {
  stopifnot(identical(x$rho$score[!is.na(pattern)], pattern[!is.na(pattern)]),
            identical(x$rho$estimated, is.na(pattern)),
            sum(grepl("^rho\\[", names(x$par_opt))) == sum(is.na(pattern)))
}
fixed_anchors(f, p)
fixed_anchors(g, q)
stopifnot(is.finite(logLik(f)), is.finite(logLik(g)))
m <- growth_tt_fit(dat, items, occ = "wave", rho_pattern = c(NA, 1, NA, 0, NA),
                   identification = "baseline_std", intercept_covariates = "age",
                   slope_covariates = "age", nquadgh = 2L, maxiter = 0L)
stopifnot(m$decode$baseline_var == 1, m$structure$zero_anchor == 9,
          !("log_zeta_var[t4,F1]" %in% names(m$par_opt)),
          is.finite(logLik(m)))

s <- growth_tt_fit(dat, items, occ = "wave", rho_pattern = cbind(p, q), simple = c(1, 1, 2, 2),
                   nquadgh = 2L, maxiter = 0L, identification = "baseline_std")
stopifnot(is.finite(logLik(s)), all(dim(s$decode$rho) == c(5L, 2L)),
          identical(as.vector(s$decode$rho[!is.na(cbind(p, q))]), c(0, 1, 0, 1)))

part <- growth_tt_individual_derivatives(f, par_index = match("rho[6,F1]", names(f$par_opt)),
                                         compute_hessian = FALSE)
stopifnot(abs(part$checks$loglik_diff) < 1e-8, all(is.finite(part$score_ind)),
          any(abs(part$score_ind) > 0), abs(part$loglik_i[["1"]]) < 1e-9)
fixed <- growth_tt_fit(dat, items, occ = "wave", rho_pattern = c(0, .25, .5, .75, 1),
                       nquadgh = 2L, maxiter = 0L, re_item = FALSE)
stopifnot(!any(fixed$rho$estimated), all(fixed$rho$se == 0))

inter <- growth_tt_fit(dat, items, occ = "wave", rho_pattern = q, fixed_covariates = "age",
                       time_interactions = "age", nquadgh = 2L, maxiter = 0L)
stopifnot(isTRUE(all.equal(inter$decode$X_effective[1, , 1, 4], -c(0, .25, .5, .75, 1))))

for (method in c("GH", "AGH")) {
  joint <- cfa_long_fit_md(dat, items, occ = "wave", methInteg = method,
                          nquadgh = 2L, maxiter = 0L, profile = FALSE)
  rule <- integration_info(joint)
  stopifnot(rule$joint_dimension == 9L, rule$largest_rule_dimension == 9L,
            rule$rules$nodes_per_rule == 2^9,
            rule$largest_mode_dimension == if (method == "AGH") 9L else 0L)

  nested <- cfa_long_fit_md(dat, items, occ = "wave", methInteg = method,
                           two_tier = TRUE, nquadgh = 2L, maxiter = 0L, profile = FALSE)
  tier <- integration_info(nested)
  stopifnot(tier$joint_dimension == 9L, identical(tier$rules$dimension, c(5L, 1L)),
            identical(tier$rules$copies, c(1L, 4L)),
            identical(tier$rules$nodes_per_rule, c(2^5, 2)),
            tier$largest_mode_dimension == if (method == "AGH") 5L else 0L)
}
stopifnot(integration_info(f)$largest_mode_dimension == 9L,
          integration_info(f)$largest_rule_dimension == 5L,
          integration_info(fixed)$joint_dimension == 5L)
