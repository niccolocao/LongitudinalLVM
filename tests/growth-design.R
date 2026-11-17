library(longProbitLVM)

set.seed(17)
N <- 25L
dat <- expand.grid(ID = seq_len(N), occ = 0:2)
u0 <- rep(rnorm(N, sd = 0.3), 3)
u1 <- rep(rnorm(N, sd = 0.2), 3)
eta <- 0.5 + 0.3 * dat$occ + u0 + u1 * dat$occ + rnorm(nrow(dat), sd = 0.5)
items <- c("a", "b", "c")
for (item in items) {
  dat[[item]] <- cut(eta + rnorm(nrow(dat)), c(-Inf, 0, 1, 2, Inf), ordered_result = TRUE)
  dat[[item]][sample.int(nrow(dat), 11L)] <- NA
}
dat[dat$ID == 1L, items] <- NA
dat <- dat[sample.int(nrow(dat)), ]

for (method in c("AGH_FULL_EQUIV", "AGH_RE_FULL_EQUIV")) {
  f <- growth_tt_fit(dat, items, methInteg = method, identification = "marker", marker_item = "a",
                     nquadgh = 2L, maxiter = 3L, growth_shape = "rho", rho_fixed = 2,
                     debug_par = FALSE, profile = FALSE)
  stopifnot(is.finite(logLik(f)), isTRUE(f$structure$analytic_score))
}

dat2 <- expand.grid(ID = 1:22, occ = 0:4)
score <- dat2$occ
latent <- 0.4 + score / 3 - score^2 / 12 + rnorm(nrow(dat2))
for (item in c("a", "b")) {
  dat2[[item]] <- as.integer(cut(latent + rnorm(nrow(dat2)), c(-Inf, 0, 1, 2, Inf)))
}
dat2[dat2$ID == 1L, "b"] <- NA
s <- 0:4
XZ <- cbind(intercept = 1, time = s, quadratic = s^2)
g <- growth_tt_fit_flexible(dat2, c("a","b"), simple_growth = FALSE,
                             X = XZ, Z = XZ, delta_free = FALSE,
                             methInteg = "GH", nquadgh = 3L, maxiter = 3L,
                             debug_par = FALSE, profile = FALSE)
stopifnot(is.finite(logLik(g)), all(dim(g$decode$Sigma_u) == c(3L, 3L)),
          all(dim(g$decode$beta) == c(1L, 3L)), g$structure$growth$p_ran == 3L)
