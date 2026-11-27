#cross-sectional CFA checks
library(LongitudinalLVM)
source(system.file("cross-sectional", "cfa_fit.R", package = "LongitudinalLVM"))

set.seed(142)
N <- 2000L
eta <- rnorm(N, .3, sqrt(1.2))
lam <- c(1, .8, 1.2, .9, 1.1)
nu <- c(0, .3, -.3, .2, -.2)
ystar <- sweep(outer(eta, lam), 2, nu, "+") +
  matrix(rnorm(N * 5L), N, 5L)
dat <- as.data.frame(1L * (ystar > 0))
names(dat) <- paste0("y", 1:5)
theta <- c(lam[-1], nu[-1], .3, 1.2)

f0 <- cfa_fit(dat, Q = 1, nquadgh = 50, start = theta, maxiter = 0)
stopifnot(f0$npar == 10L, nrow(f0$patterns) == 32L,
          abs(sum(f0$probabilities) - 1) < 1e-10,
          abs(min(N * f0$probabilities) - 11.02225) < .02,
          max(abs(f0$coefficients - theta)) < 1e-10,
          max(abs(colSums(f0$jacobian))) < 1e-7,
          max(abs(crossprod(f0$pattern_scores, f0$probabilities))) < 1e-7,
          max(abs(f0$jacobian - f0$pattern_scores * f0$probabilities)) < 1e-12,
          all(dim(f0$scores) == c(N, 10L)),
          is.na(f0$convergence))

p_ref <- apply(f0$patterns, 1, function(y) {
  integrate(function(z) {
    vapply(z, function(zz) {
      pp <- pnorm(nu + lam * (.3 + sqrt(1.2) * zz))
      prod(ifelse(y == 1, pp, 1 - pp)) * dnorm(zz)
    }, numeric(1))
  }, -Inf, Inf, rel.tol = 1e-9)$value
})
stopifnot(max(abs(p_ref - f0$probabilities)) < 1e-7)

fit <- cfa_fit(dat, Q = 1, nquadgh = 50, start = theta, maxiter = 500)
stopifnot(fit$convergence == 0L,
          fit$loglik >= f0$loglik - 1e-6,
          max(abs(colMeans(fit$scores))) < 1e-3,
          all(eigen(fit$information, symmetric = TRUE)$values > 0))

set.seed(143)
Phi <- matrix(c(1.2, .35, .35, .9), 2)
eta2 <- sweep(matrix(rnorm(N * 2), N, 2) %*% chol(Phi), 2, c(.2, -.1), "+")
lambda2 <- rbind(c(1, 0), c(.8, 0), c(1.1, 0),
                 c(0, 1), c(0, .9), c(0, 1.2))
nu2 <- c(0, .2, -.2, 0, -.1, .3)
y2 <- sweep(eta2 %*% t(lambda2), 2, nu2, "+") +
  matrix(rnorm(N * 6), N, 6)
dat2 <- as.data.frame(1L * (y2 > 0))
names(dat2) <- paste0("y", 1:6)
theta2 <- c(.8, 1.1, .9, 1.2, .2, -.2, -.1, .3,
            .2, -.1, 1.2, .35, .9)
f2 <- cfa_fit(dat2, simple = rep(1:2, each = 3),
              nquadgh = 25, start = theta2, maxiter = 0)
stopifnot(f2$npar == 13L, nrow(f2$patterns) == 64L,
          abs(sum(f2$probabilities) - 1) < 1e-10,
          max(abs(f2$decode$Phi - Phi)) < 1e-10,
          max(abs(f2$decode$lambda - lambda2)) < 1e-10,
          max(abs(colSums(f2$jacobian))) < 1e-7)

#ordinal items and marginalization
dat3 <- dat
dat3$y2 <- as.integer(cut(ystar[, 2], c(-Inf, 0, 1, Inf),
                         labels = FALSE)) - 1L
dat3[1, 1] <- NA
theta3 <- c(theta, 1)
f3 <- cfa_fit(dat3, Q = 1, nquadgh = 50, start = theta3, maxiter = 0)
stopifnot(f3$npar == 11L, nrow(f3$patterns) == 48L,
          abs(sum(f3$probabilities) - 1) < 1e-10,
          all(dim(f3$scores) == c(N, 11L)))
pat <- f3$patterns[1:2, , drop = FALSE]
pat[1, 1] <- NA
p_marg <- f3$probability_function(f3$coefficients, pat)[1]
idx <- rowSums(sweep(f3$patterns[, -1, drop = FALSE], 2,
                     pat[1, -1], "!=")) == 0
stopifnot(abs(p_marg - sum(f3$probabilities[idx])) < 1e-10)

#resource guards
f_small <- suppressWarnings(cfa_fit(dat, start = theta, maxiter = 0,
                                   max_patterns = 2L))
stopifnot(is.null(f_small$patterns), nrow(f_small$scores) == N)
err <- try(cfa_fit(dat2, simple = rep(1:2, each = 3),
                   nquadgh = 30, max_nodes = 10, maxiter = 0),
           silent = TRUE)
stopifnot(inherits(err, "try-error"))


#adaptive GH at identical natural parameters
fa <- cfa_fit(dat, Q = 1, nquadgh = 20, methInteg = "AGH",
              start = theta, maxiter = 0)
stopifnot(fa$integration$adaptive,
          fa$integration$nodes == 20L,
          abs(sum(fa$probabilities) - 1) < 1e-7,
          max(abs(fa$probabilities - p_ref)) < 1e-7,
          abs(fa$loglik - f0$loglik) < 1e-4,
          max(abs(fa$jacobian - f0$jacobian)) < 1e-6,
          max(abs(fa$scores - f0$scores)) < 1e-5)

fa2 <- cfa_fit(dat2, simple = rep(1:2, each = 3),
               nquadgh = 20, methInteg = "AGH",
               start = theta2, maxiter = 0)
stopifnot(fa2$integration$nodes == 400L,
          max(abs(fa2$probabilities - f2$probabilities)) < 1e-6,
          max(abs(fa2$jacobian - f2$jacobian)) < 1e-5,
          abs(sum(fa2$probabilities) - 1) < 1e-6)

fa3 <- cfa_fit(dat3, Q = 1, nquadgh = 20, methInteg = "AGH",
               start = theta3, maxiter = 0)
stopifnot(max(abs(fa3$probabilities - f3$probabilities)) < 1e-7,
          max(abs(fa3$scores - f3$scores)) < 1e-5)
p_marg_a <- fa3$probability_function(fa3$coefficients, pat)[1]
stopifnot(abs(p_marg_a - p_marg) < 1e-7)
missing_pat <- matrix(NA_real_, 1L, ncol(dat3))
stopifnot(fa3$probability_function(fa3$coefficients, missing_pat) == 1,
          max(abs(fa3$jacobian_function(fa3$coefficients, missing_pat))) == 0)

#AGH optimization and resource guard
fit_a <- cfa_fit(dat, Q = 1, nquadgh = 16, methInteg = "AGH",
                 start = fit$coefficients, maxiter = 100)
stopifnot(fit_a$convergence == 0L,
          abs(fit_a$loglik - fit$loglik) < 1e-3,
          max(abs(fit_a$coefficients - fit$coefficients)) < 5e-3,
          max(abs(colMeans(fit_a$scores))) < 1e-3)
err_a <- try(cfa_fit(dat2, simple = rep(1:2, each = 3),
                     nquadgh = 30, methInteg = "AGH",
                     max_nodes = 10, maxiter = 0), silent = TRUE)
stopifnot(inherits(err_a, "try-error"))
